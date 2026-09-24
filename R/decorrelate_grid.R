#' Convert a printable grid \code{data.frame} back into a candidates list
#'
#' The inverse of \code{\link{get_decorrelate_grid_parameters}()}: parses each
#' row's covariance-type/parameter columns back into \code{tailup_initial}/
#' \code{taildown_initial}/\code{euclid_initial}/\code{nugget_initial}
#' (and, if present, \code{randcov_initial}) objects with \code{known =
#' "given"}, so a user-edited or re-supplied grid can be validated and
#' evaluated the same way as a freshly constructed one. A row where every
#' type column reads \code{"no transformation"} is treated as the
#' untransformed IID baseline.
#'
#' @param grid A grid \code{data.frame}, as returned by
#'   \code{\link{ssn_decorrelate_grid}()} or \code{\link{tidy.ssn_decorrelate_grid}()}.
#'
#' @return A named list of candidates (one per grid row), each a list of
#'   \code{*_initial} objects, suitable for
#'   \code{\link{get_decorrelate_grid_candidates}()}.
#'
#' @noRd
get_decorrelate_candidates_from_table <- function(grid) {
  components <- c("tailup", "taildown", "euclid", "nugget")
  types <- paste0(components, "_type")
  if (!NROW(grid) || anyDuplicated(names(grid)) || !all(types %in% names(grid))) {
    stop("grid must have rows and tailup_type, taildown_type, euclid_type, and nugget_type columns.", call. = FALSE)
  }
  constructors <- list(tailup = tailup_initial, taildown = taildown_initial,
                       euclid = euclid_initial, nugget = nugget_initial)
  random_columns <- names(grid)[startsWith(names(grid), "randcov_")]
  candidates <- lapply(seq_len(NROW(grid)), function(i) {
    row <- grid[i, , drop = FALSE]
    baseline <- all(vapply(row[types], function(x) identical(as.character(x), "no transformation"), logical(1)))
    candidate <- list()
    for (component in components) {
      type <- as.character(row[[paste0(component, "_type")]])
      if (baseline) type <- if (component == "nugget") "nugget" else "none"
      if (length(type) != 1L || is.na(type)) stop("grid covariance types must be nonmissing.", call. = FALSE)
      args <- list(type)
      if (type != "none") {
        parameters <- if (component == "nugget") "nugget" else c("de", "range")
        if (component == "euclid") {
          if (type %in% c("matern", "cauchy", "pexponential")) parameters <- c(parameters, "extra")
          args$rotate <- if ("euclid_rotate" %in% names(row)) row$euclid_rotate else 0
          args$scale <- if ("euclid_scale" %in% names(row)) row$euclid_scale else 1
        }
        for (parameter in parameters) {
          column <- paste(component, parameter, sep = "_")
          value <- if (baseline) 1 else row[[column]]
          if (is.null(value) || !is.numeric(value) || length(value) != 1L || !is.finite(value)) {
            stop("grid needs a finite numeric ", column, " for row ", i, ".", call. = FALSE)
          }
          args[[parameter]] <- value
        }
        args$known <- "given"
      }
      candidate[[paste0(component, "_initial")]] <- do.call(constructors[[component]], args)
    }
    if (length(random_columns)) {
      values <- as.list(row[random_columns])
      if (baseline) values[] <- 0
      names(values) <- sub("^randcov_", "", random_columns)
      candidate$randcov_initial <- do.call(randcov_initial, c(values, list(known = "given")))
    }
    candidate
  })
  names(candidates) <- paste0("grid", seq_along(candidates))
  candidates
}

#' Validate a required single-logical grid-construction flag
#'
#' @param value The value to validate.
#' @param name The argument's name, used in the error message.
#'
#' @return \code{NULL}, invisibly, if \code{value} is a single non-missing
#'   logical; otherwise an error.
#'
#' @noRd
check_decorrelate_grid_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

#' Construct an automatic joint SSN covariance candidate grid
#'
#' Builds a heuristic candidate grid of decorrelation/estimation covariance
#' parameters: an OLS residual variance anchors overall scale, active
#' components (tailup/taildown/euclid/nugget/random effects) split that
#' variance across targeted allocation regimes, and each active spatial
#' component is crossed with two stream- or Euclidean-distance-based range
#' choices (and, if \code{anisotropy}, rotate/scale choices). Fixed
#' (\code{is_known}) fields in \code{initial_object}/\code{randcov_initial}
#' are pinned rather than varied. Called by
#' \code{\link{ssn_decorrelate_grid_internal}()} (with \code{add_iid = TRUE})
#' and by \code{\link{ssn_decorrelate}()}'s automatic-grid-construction path
#' (with \code{add_iid = FALSE}, adding the baseline separately via
#' \code{\link{add_decorrelate_iid}()} only when needed).
#'
#' @param formula A model formula.
#' @param ssn.object A fitted-model-ready SSN object.
#' @param initial_object A joint covariance initial-value object (tailup,
#'   taildown, euclid, nugget), with any fixed fields already marked known.
#' @param additive The additive function value column name, or \code{NULL}.
#' @param anisotropy Whether Euclidean anisotropy is active.
#' @param random A one- or two-sided random effect formula, or \code{NULL}.
#' @param randcov_initial A random-effect variance initial-value object, or
#'   \code{NULL}.
#' @param dense_grid Whether to use the denser grid density, or \code{NULL}
#'   to resolve it from the observed sample size.
#' @param add_iid Whether to append an untransformed baseline candidate.
#'
#' @return A named list of candidates, each a list of \code{*_initial}
#'   objects.
#'
#' @noRd
get_decorrelate_grid <- function(formula, ssn.object, initial_object, additive,
                                  anisotropy, random, randcov_initial,
                                  dense_grid, add_iid) {
  observed_index <- get_decorrelate_response_index(formula, ssn.object)
  observed <- ssn.object$obs[observed_index, , drop = FALSE]
  if (is.null(dense_grid)) dense_grid <- NROW(observed) <= 5000L
  check_decorrelate_grid_flag(dense_grid, "dense_grid")
  check_decorrelate_grid_flag(add_iid, "add_iid")
  matrix_object <- get_model_matrix_object(formula, observed, attr(observed, "sf_column"))
  frame <- matrix_object$obdata_model_frame
  y <- as.numeric(model.response(frame))
  offset <- model.offset(frame)
  if (!is.null(offset)) y <- y - offset
  ols <- stats::lm.fit(matrix_object$X, y)
  s2 <- sum(ols$residuals^2) / max(1L, length(y) - ols$rank)
  if (!is.finite(s2) || s2 <= 0) s2 <- 1
  ns2 <- 1.2 * s2

  randcov_names <- get_randcov_names(random)
  if (length(randcov_names)) {
    if (is.null(randcov_initial)) randcov_initial <- randcov_initial()
    given_names <- unlist(lapply(names(randcov_initial$initial), function(name) {
      get_randcov_names(reformulate(name))
    }), use.names = FALSE)
    if (length(given_names) != length(randcov_initial$initial) ||
        any(!given_names %in% randcov_names)) {
      stop("randcov_initial must specify individual terms from random.", call. = FALSE)
    }
    names(randcov_initial$initial) <- names(randcov_initial$is_known) <- given_names
  } else {
    randcov_initial <- NULL
  }
  initial_object$randcov_initial <- randcov_initial
  initial_object <- get_initial_NA_object(
    initial_object, list(anisotropy = anisotropy, randcov_names = randcov_names)
  )
  components <- c("tailup", "taildown", "euclid", "nugget")
  types <- vapply(initial_object[paste0(components, "_initial")], function(x) {
    remove_covtype(class(x)[1])
  }, character(1))
  names(types) <- components
  active <- components[types != "none"]
  variance_names <- c(ifelse(active == "nugget", "nugget", paste0(active, "_de")), randcov_names)
  proportions <- get_decorrelate_variance_grid(length(variance_names), dense_grid)
  variances <- as.data.frame(ns2 * proportions)
  names(variances) <- variance_names
  for (name in setdiff(c("tailup_de", "taildown_de", "euclid_de", "nugget"), names(variances))) {
    variances[[name]] <- 0
  }

  bbox <- sf::st_bbox(observed)
  euclid_max <- sqrt((bbox[["xmax"]] - bbox[["xmin"]])^2 + (bbox[["ymax"]] - bbox[["ymin"]])^2)
  tail_max <- 0
  if (any(types[c("tailup", "taildown")] != "none")) {
    # Match the stream-distance proxy used by local model fitting.
    tail_max <- euclid_max
  }
  if (any(types[c("tailup", "taildown")] != "none") && NROW(observed) <= 5000L) {
    distance_initial <- initial_object
    distance_initial$euclid_initial <- euclid_initial("none")
    backend <- select_square_dist_backend(ssn.object, "obs", types[["tailup"]] == "none", types[["taildown"]] == "none")
    distances <- get_dist_object(ssn.object, distance_initial, additive, FALSE, backend)
    hydro <- distances$hydro_mat[observed_index, observed_index, drop = FALSE]
    mask <- distances$mask_mat[observed_index, observed_index, drop = FALSE]
    tail_max <- max(hydro * mask)
  }
  effective <- list(
    stream = pmax(c(0.25, 0.75) * tail_max, .Machine$double.eps),
    euclid = pmax(c(0.25, 0.75) * euclid_max, .Machine$double.eps)
  )
  tail_range <- function(component, constructor) {
    type <- types[[component]]
    if (type == "none") return(Inf)
    get_range_from_effective(constructor(type, de = 1, range = 1), effective$stream)
  }
  euclid_range <- if (types[["euclid"]] == "none") Inf else effective$euclid
  if (types[["euclid"]] != "none" && !euclid_has_extra(types[["euclid"]])) {
    euclid_range <- get_range_from_effective(
      euclid_params(types[["euclid"]], de = 1, range = 1), effective$euclid
    )
  }
  ranges <- list(
    tailup_range = tail_range("tailup", tailup_params),
    taildown_range = tail_range("taildown", taildown_params),
    euclid_range = euclid_range
  )
  shape_grid <- do.call(expand.grid, ranges)
  if (euclid_has_extra(types[["euclid"]])) {
    extra <- switch(types[["euclid"]], matern = c(1, 4), cauchy = c(0.5, 2), pexponential = c(0.4, 1.6))
    supplied_extra <- initial_object$euclid_initial$initial[["extra"]]
    if (!is.na(supplied_extra)) extra <- supplied_extra
    shape_grid <- merge(shape_grid, data.frame(euclid_extra = extra), by = NULL)
    for (value in extra) {
      rows <- shape_grid$euclid_extra == value
      parameters <- euclid_params(types[["euclid"]], de = 1, range = 1, extra = value)
      shape_grid$euclid_range[rows] <- get_range_from_effective(
        parameters, shape_grid$euclid_range[rows]
      )
    }
  }
  if (anisotropy && types[["euclid"]] != "none") {
    rotate <- if (dense_grid) c(0, 45, 90, 135) * pi / 180 else c(0, pi / 2)
    scale <- if (dense_grid) c(0.25, 0.5, 0.75, 1) else c(0.5, 1)
  } else {
    rotate <- 0
    scale <- 1
  }
  shape_grid <- merge(shape_grid, expand.grid(rotate = rotate, scale = scale), by = NULL)
  cov_grid <- merge(variances, shape_grid, by = NULL)
  cov_grid <- cov_grid_replace_shared(cov_grid, initial_object, list(randcov_names = randcov_names))
  cov_grid$rotate[cov_grid$scale == 1] <- 0
  cov_grid <- unique(cov_grid)
  candidates <- lapply(seq_len(NROW(cov_grid)), function(i) {
    values <- cov_grid[i, , drop = FALSE]
    candidate <- initial_object
    for (component in components) {
      initial_name <- paste0(component, "_initial")
      for (name in names(candidate[[initial_name]]$initial)) {
        column <- if (component == "nugget" || name %in% c("rotate", "scale")) name else paste(component, name, sep = "_")
        candidate[[initial_name]]$initial[[name]] <- values[[column]]
      }
      candidate[[initial_name]]$is_known[] <- TRUE
    }
    if (length(randcov_names)) {
      candidate$randcov_initial$initial <- unlist(values[randcov_names], use.names = TRUE)
      candidate$randcov_initial$is_known[] <- TRUE
    }
    candidate
  })
  names(candidates) <- paste0("grid", seq_along(candidates))
  if (add_iid) candidates <- add_decorrelate_iid(candidates, random)
  candidates
}

#' Build variance-allocation proportions for an automatic covariance grid
#'
#' For one active variance component, returns the trivial single allocation
#' (proportion 1). For two, crosses a set of candidate proportions with their
#' complements. For three or more, avoids a combinatorial cross by varying
#' one dominant component at a time (the rest split its complement evenly),
#' plus an evenly-split row.
#'
#' @param n The number of active variance components (spatial plus random
#'   effects).
#' @param dense_grid Whether to use the denser candidate-proportion set.
#'
#' @return A matrix with \code{n} columns of variance proportions (rows
#'   summing to 1).
#'
#' @noRd
get_decorrelate_variance_grid <- function(n, dense_grid) {
  if (n == 1L) return(matrix(1, nrow = 1L))
  proportions <- if (dense_grid) c(0.05, 0.25, 0.5, 0.75, 0.95) else c(0.5, 0.95)
  if (n == 2L) return(cbind(proportions, 1 - proportions))
  # Vary each dominant component without a Cartesian product of variance shares.
  rows <- lapply(seq_len(n), function(i) {
    values <- matrix(rep((1 - proportions) / (n - 1L), n), ncol = n)
    values[, i] <- proportions
    values
  })
  unique(rbind(rep(1 / n, n), do.call(rbind, rows)))
}

#' Append an untransformed IID baseline candidate to a candidates list
#'
#' Builds a candidate with no spatial covariance, unit nugget variance, and
#' (if \code{random} is given) zero random-effect variance, all known, and
#' appends it under a unique name.
#'
#' @param candidates A named list of candidates to append to.
#' @param random A one- or two-sided random effect formula, or \code{NULL}.
#'
#' @return \code{candidates} with the baseline candidate appended.
#'
#' @noRd
add_decorrelate_iid <- function(candidates, random) {
  iid <- list(
    tailup_initial = tailup_initial("none"),
    taildown_initial = taildown_initial("none"),
    euclid_initial = euclid_initial("none"),
    nugget_initial = nugget_initial("nugget", nugget = 1, known = "given"),
    randcov_initial = NULL
  )
  randcov_names <- get_randcov_names(random)
  if (length(randcov_names)) {
    values <- as.list(stats::setNames(rep(0, length(randcov_names)), randcov_names))
    iid$randcov_initial <- do.call(randcov_initial, c(values, list(known = "given")))
  }
  name <- tail(make.unique(c(names(candidates), "iid")), 1L)
  candidates[[name]] <- iid
  candidates
}
