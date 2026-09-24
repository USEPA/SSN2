#' Build a reusable partition-factor context: grouping formula and
#' training-side group labels/levels
#'
#' Reuses training labels and row indices across prediction chunks. Supply
#' all prediction rows in \code{newdata} so the extended factor levels cover
#' every chunk, then pass the context to \code{\link{partition_vector}()}.
#'
#' @param partition_factor A partition-factor formula, or \code{NULL}.
#' @param data The training data.
#' @param newdata Data to extend factor levels against (e.g. every row the
#'   context will ever be used with, not just one chunk).
#' @param reform_bar2 An already-built grouping formula, or \code{NULL} (built
#'   here).
#' @param xlev Optional factor levels to enforce on the training-side model
#'   frame.
#'
#' @return A list with \code{reform_bar2} and \code{partition_index_data}
#'   (itself a list with \code{reform_bar2_vals}, \code{reform_bar2_xlev}, and
#'   \code{level_index_map} -- the group-label-to-training-row-index lookup
#'   \code{\link{partition_vector}()} otherwise rebuilds on every call), in
#'   the shape \code{\link{partition_vector}()} accepts back.
#'
#' @noRd
get_partition_context <- function(partition_factor, data, newdata, reform_bar2 = NULL, xlev = NULL) {
  if (is.null(reform_bar2)) {
    partition_factor_val <- get_randcov_name(labels(terms(partition_factor)))
    bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE, env = asNamespace("SSN2"))
  }
  p_index_data_mf <- model.frame(reform_bar2, data, xlev = xlev)
  p_index_data_vals <- model_matrix_group_labels(reform_bar2, data, xlev = xlev)
  p_index_data_xlev <- .getXlevels(terms(p_index_data_mf), p_index_data_mf)
  p_index_data_xlev_full <- .getXlevels(terms(p_index_data_mf), rbind(p_index_data_mf, model.frame(reform_bar2, newdata)))
  if (!identical(p_index_data_xlev, p_index_data_xlev_full)) {
    p_index_data_xlev <- p_index_data_xlev_full
  }
  # Cache training-row matches for reuse across chunks.
  level_index_map <- split(seq_along(p_index_data_vals), p_index_data_vals)
  list(
    reform_bar2 = reform_bar2,
    partition_index_data = list(
      reform_bar2_vals = p_index_data_vals, reform_bar2_xlev = p_index_data_xlev,
      level_index_map = level_index_map
    )
  )
}

partition_vector <- function(partition_factor, data, newdata, reform_bar2 = NULL, partition_index_data = NULL, xlev = NULL) {
  if (is.null(partition_factor)) {
    t_partition_index <- NULL
  } else {
    if (is.null(reform_bar2) || is.null(partition_index_data)) {
      built <- get_partition_context(partition_factor, data, newdata, reform_bar2 = reform_bar2, xlev = xlev)
      reform_bar2 <- built$reform_bar2
      if (is.null(partition_index_data)) partition_index_data <- built$partition_index_data
      # partition_index_data <- as.vector(model.matrix(reform_bar2, data))
    }
    partition_index_newdata <- model_matrix_group_labels(reform_bar2, newdata, xlev = partition_index_data$reform_bar2_xlev, na_pass = TRUE)

    # Group lookup avoids a dense prediction-by-observation comparison.
    group_label <- partition_index_data$reform_bar2_vals
    n_obs <- length(group_label)
    n_new <- length(partition_index_newdata)
    # Callers without a cached map still use the same grouping semantics.
    level_index_map <- partition_index_data$level_index_map
    if (is.null(level_index_map)) {
      level_index_map <- split(seq_along(group_label), group_label)
    }
    non_na_new <- !is.na(partition_index_newdata)
    matches <- vector("list", n_new)
    matches[non_na_new] <- level_index_map[partition_index_newdata[non_na_new]]
    match_lengths <- lengths(matches)
    obs_idx <- unlist(matches, use.names = FALSE)
    new_idx <- rep(seq_len(n_new), match_lengths)

    t_partition_index <- Matrix::sparseMatrix(i = new_idx, j = obs_idx, x = rep(1, length(obs_idx)), dims = c(n_new, n_obs))
  }
  t_partition_index
}

get_partition_xlev <- function(partition_factor, obdata) {
  if (is.null(partition_factor)) {
    return(NULL)
  }
  partition_factor_val <- get_randcov_name(labels(terms(partition_factor)))
  bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
  reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE, env = asNamespace("SSN2"))
  mf <- model.frame(reform_bar2, obdata)
  .getXlevels(terms(mf), mf)
}

extend_partition_xlev <- function(partition_factor, partition_xlev, newdata) {
  if (is.null(partition_factor) || is.null(partition_xlev)) {
    return(partition_xlev)
  }
  partition_factor_val <- get_randcov_name(labels(terms(partition_factor)))
  bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
  reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE, env = asNamespace("SSN2"))
  newdata_mf <- model.frame(reform_bar2, newdata, na.action = na.pass)
  newdata_xlev <- .getXlevels(terms(newdata_mf), newdata_mf)
  mapply(
    function(old, variable) union(old, newdata_xlev[[variable]]),
    partition_xlev, names(partition_xlev), SIMPLIFY = FALSE
  )
}
