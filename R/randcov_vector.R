randcov_vector <- function(randcov_params = NULL, data, newdata, xlev_list = NULL, context = NULL) {
  if (is.null(randcov_params)) {
    randcov_vectors <- NULL
  } else {
    randcov_names <- names(randcov_params)
    randcov_vectors <- lapply(randcov_names, get_randcov_vectors, randcov_params, data, newdata, xlev_list, context)
    randcov_vectors <- Reduce("+", randcov_vectors)
  }
  randcov_vectors
}

#' Precompute the training-side pieces of every random-effect term once
#'
#' \code{\link{get_randcov_vectors}()} re-derives the training-side grouping
#' labels (and, for a slope term, the training-side covariate values) from
#' \code{data} on every call. When the same \code{data}/\code{randcov_params}
#' pair is reused across many calls (e.g. once per row-chunk of a large
#' prediction set), that training-side work is identical every time; this
#' builds it once so a chunked caller can pass the result back in via
#' \code{\link{get_randcov_vectors}()}'s \code{context} argument.
#'
#' @param randcov_params A named vector of random-effect variances, or
#'   \code{NULL}.
#' @param data The training data.
#' @param newdata_all Every row the resulting context will ever be used
#'   against (e.g. a full prediction set before it is split into chunks) --
#'   used only to extend factor levels the training data alone does not
#'   contain, matching \code{\link{get_randcov_vectors}()}'s own per-call
#'   extension.
#' @param xlev_list A named list of factor levels, one entry per random-effect
#'   term, or \code{NULL}.
#'
#' @return A named list (one entry per random-effect term) of training-side
#'   pieces, or \code{NULL} if \code{randcov_params} is \code{NULL}.
#'
#' @noRd
get_randcov_context <- function(randcov_params, data, newdata_all, xlev_list = NULL) {
  if (is.null(randcov_params)) {
    return(NULL)
  }
  randcov_names <- names(randcov_params)
  context <- lapply(randcov_names, function(randcov_name) {
    bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
    # Preserve package lookup without retaining training-side temporaries.
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE, env = asNamespace("SSN2"))
    base_xlev <- xlev_list[[randcov_name]]

    Z_index_data_mf <- model.frame(reform_bar2, data, xlev = base_xlev)
    Z_index_data <- model_matrix_group_labels(reform_bar2, data, xlev = base_xlev)
    Z_index_data_xlev <- .getXlevels(terms(Z_index_data_mf), Z_index_data_mf)

    Z_index_data_xlev_full <- .getXlevels(terms(Z_index_data_mf), rbind(Z_index_data_mf, model.frame(reform_bar2, newdata_all, xlev = base_xlev)))
    if (!identical(Z_index_data_xlev, Z_index_data_xlev_full)) {
      Z_index_data_xlev <- Z_index_data_xlev_full
    }

    if (bar_split[[1]] != "1") {
      reform_bar1 <- reformulate(bar_split[[1]], intercept = FALSE, env = asNamespace("SSN2"))
      Z_val_data <- as.vector(model.matrix(reform_bar1, data))
    } else {
      reform_bar1 <- NULL
      Z_val_data <- NULL
    }

    # the new-data group lookup below matches against this map on every
    # call; precomputing it here lets a chunked caller build it once
    level_index_map <- split(seq_along(Z_index_data), Z_index_data)

    list(
      reform_bar2 = reform_bar2, reform_bar1 = reform_bar1,
      Z_index_data = Z_index_data, Z_index_data_xlev = Z_index_data_xlev,
      Z_val_data = Z_val_data, level_index_map = level_index_map
    )
  })
  names(context) <- randcov_names
  context
}

get_randcov_vectors <- function(randcov_name, randcov_params, data, newdata, xlev_list = NULL, context = NULL) {
  randcov_param <- as.numeric(randcov_params[randcov_name])
  term_context <- context[[randcov_name]]

  if (is.null(term_context)) {
    bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE, env = asNamespace("SSN2"))
    base_xlev <- xlev_list[[randcov_name]]

    Z_index_data_mf <- model.frame(reform_bar2, data, xlev = base_xlev)
    Z_index_data <- model_matrix_group_labels(reform_bar2, data, xlev = base_xlev)
    Z_index_data_xlev <- .getXlevels(terms(Z_index_data_mf), Z_index_data_mf)

    Z_index_data_xlev_full <- .getXlevels(terms(Z_index_data_mf), rbind(Z_index_data_mf, model.frame(reform_bar2, newdata, xlev = base_xlev)))
    if (!identical(Z_index_data_xlev, Z_index_data_xlev_full)) {
      Z_index_data_xlev <- Z_index_data_xlev_full
    }

    if (bar_split[[1]] != "1") {
      reform_bar1 <- reformulate(bar_split[[1]], intercept = FALSE, env = asNamespace("SSN2"))
      Z_val_data <- as.vector(model.matrix(reform_bar1, data))
    } else {
      reform_bar1 <- NULL
      Z_val_data <- NULL
    }
    level_index_map <- split(seq_along(Z_index_data), Z_index_data)
  } else {
    reform_bar2 <- term_context$reform_bar2
    reform_bar1 <- term_context$reform_bar1
    Z_index_data <- term_context$Z_index_data
    Z_index_data_xlev <- term_context$Z_index_data_xlev
    Z_val_data <- term_context$Z_val_data
    level_index_map <- term_context$level_index_map
  }

  Z_index_newdata <- model_matrix_group_labels(reform_bar2, newdata, xlev = Z_index_data_xlev, na_pass = TRUE)

  Z_val_newdata <- if (!is.null(reform_bar1)) get_randcov_slope_val_newdata(reform_bar1, newdata) else NULL

  # look up, for each newdata row, the (typically few) data rows sharing its
  # group label via a hash/list lookup instead of comparing against every
  # data row elementwise -- scales with match count rather than
  # length(Z_index_data) x length(Z_index_newdata), and matches the pattern
  # already used for partition_vector()'s own grouping match
  n_obs <- length(Z_index_data)
  n_new <- length(Z_index_newdata)
  non_na_new <- !is.na(Z_index_newdata)
  matches <- vector("list", n_new)
  matches[non_na_new] <- level_index_map[Z_index_newdata[non_na_new]]
  match_lengths <- lengths(matches)
  obs_idx <- unlist(matches, use.names = FALSE)
  new_idx <- rep(seq_len(n_new), match_lengths)

  x_vals <- rep(randcov_param, length(obs_idx))
  if (!is.null(Z_val_data)) {
    x_vals <- x_vals * Z_val_data[obs_idx] * Z_val_newdata[new_idx]
  }

  Matrix::sparseMatrix(i = new_idx, j = obs_idx, x = x_vals, dims = c(n_new, n_obs))
}

get_randcov_slope_val_newdata <- function(reform_bar1, newdata) {
  slope_val_newdata <- as.vector(model.matrix(reform_bar1, model.frame(reform_bar1, newdata, na.action = na.pass)))
  if (anyNA(slope_val_newdata)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }
  slope_val_newdata
}

randcov_newvar <- function(randcov_params = NULL, newdata_row, context = NULL) {
  if (is.null(randcov_params)) {
    return(0)
  }
  randcov_names <- names(randcov_params)
  vars <- vapply(randcov_names, function(randcov_name) {
    randcov_param <- as.numeric(randcov_params[randcov_name])
    term_context <- context[[randcov_name]]
    if (is.null(term_context)) {
      bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
      if (bar_split[[1]] == "1") {
        return(randcov_param)
      }
      reform_bar1 <- reformulate(bar_split[[1]], intercept = FALSE, env = asNamespace("SSN2"))
    } else {
      reform_bar1 <- term_context$reform_bar1
      if (is.null(reform_bar1)) {
        return(randcov_param)
      }
    }
    slope_val_newdata <- get_randcov_slope_val_newdata(reform_bar1, newdata_row)
    randcov_param * slope_val_newdata^2
  }, numeric(1))
  sum(vars)
}

get_spatial_nugget_var <- function(params_object, diagtol) {
  de_scale <- sum(params_object$tailup[["de"]], params_object$taildown[["de"]], params_object$euclid[["de"]])
  nugget_val <- max(params_object$nugget[["nugget"]], 1e-4 * de_scale, diagtol)
  de_scale + nugget_val
}
