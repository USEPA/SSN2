#' Evaluate Laplace log likelihood
#'
#' @param par Current optimization parameter
#' @param orig2optim_object Optimization object that controls transformations
#' @param data_object Data object
#' @param estmethod Estimation method
#'
#' @noRd
laploglik <- function(par, orig2optim_object, data_object, estmethod) {
  # find covariance parameters on original scale
  cov_orig_val <- optim2orig_glm(orig2optim_object, par)

  # make params object
  params_object <- get_params_object_glm(orig2optim_object$classes, cov_orig_val)

  # find products
  lapll_prods <- laploglik_products(params_object, data_object, estmethod)

  # find minus two lap log lik
  minustwolaploglik <- get_minustwolaploglik(lapll_prods, data_object, estmethod)

  # find second quadrant if anisotropy (spmodel has laploglik_anis)
  if (data_object$anisotropy && !orig2optim_object$is_known[["euclid_rotate_is_known"]]) {
    # change rotate parameter
    params_object_q2 <- params_object
    params_object_q2$euclid[["rotate"]] <- pi - params_object_q2$euclid[["rotate"]]

    # find products
    lapll_prods_q2 <- laploglik_products(params_object_q2, data_object, estmethod)

    # find minus two lap log lik
    minustwolaploglik_q2 <- get_minustwolaploglik(lapll_prods_q2, data_object, estmethod)

    # find minimum rotate
    minustwolaploglik <- min(c(minustwolaploglik, minustwolaploglik_q2))
  }

  # return
  minustwolaploglik
}
