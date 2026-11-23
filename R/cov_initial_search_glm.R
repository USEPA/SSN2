cov_initial_search_glm <- function(initial_NA_object, ssn.object, data_object, estmethod) {
  cov_grid <- build_cov_initial_grid(initial_NA_object, data_object, is_glm = TRUE)
  # split into list
  cov_grid_splits <- split(cov_grid, seq_len(NROW(cov_grid)))
  # iterate through list
  objvals <- vapply(cov_grid_splits, function(x) eval_grid_glm(x, initial_NA_object, ssn.object, data_object, estmethod), numeric(1))
  # find parameters that yield the minimum -2ll
  min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
  # store this as new NA object
  updated_NA_object <- initial_NA_object
  updated_NA_object$tailup_initial$initial <- c(de = min_params[["tailup_de"]], range = min_params[["tailup_range"]])
  updated_NA_object$taildown_initial$initial <- c(de = min_params[["taildown_de"]], range = min_params[["taildown_range"]])
  updated_NA_object$euclid_initial$initial <- update_euclid_grid_initial(
    updated_NA_object$euclid_initial, min_params
  )
  updated_NA_object$nugget_initial$initial <- c(nugget = min_params[["nugget"]])
  updated_NA_object$dispersion_initial$initial <- c(dispersion = min_params[["dispersion"]])

  if (!is.null(updated_NA_object$randcov_initial)) {
    updated_NA_object$randcov_initial$initial <- min_params[data_object$randcov_names]
  }

  # return best parameters
  best_params <- list(initial_object = updated_NA_object)
}

eval_grid_glm <- function(cov_grid_split, initial_NA_object, ssn.object, data_object, estmethod) {
  cov_grid <- unlist(cov_grid_split)
  # params object
  params_object <- get_params_object_grid_glm(cov_grid, initial_NA_object)

  # gloglik products
  lapll_prods <- laploglik_products(params_object, data_object, estmethod)

  # minus two gloglik
  get_minustwolaploglik(lapll_prods, data_object, estmethod)
}

cov_grid_replace_glm <- function(cov_grid, initial_object, data_object) {
  cov_grid <- cov_grid_replace_shared(cov_grid, initial_object, data_object)

  if (!is.na(initial_object$dispersion_initial$initial[["dispersion"]])) {
    cov_grid[, "dispersion"] <- initial_object$dispersion_initial$initial[["dispersion"]]
  }

  cov_grid
}
