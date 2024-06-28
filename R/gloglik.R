#' Compute log likeihood
#'
#' @param par Optimization parameters
#' @param orig2optim_object An optimization object that performs necessary transformations
#' @param data_object Data object
#' @param estmethod Estimation method
#'
#' @noRd
gloglik <- function(par, orig2optim_object, data_object, estmethod) {
  # find parameters on the original scale
  cov_orig_val <- optim2orig(orig2optim_object, par)

  # crearte a parameter object
  params_object <- get_params_object(orig2optim_object$classes, cov_orig_val)

  # find the -2ll products
  gll_prods <- gloglik_products(params_object, data_object, estmethod)

  # compute the -2ll
  minustwologlik <- get_minustwologlik(gll_prods, data_object, estmethod)

  # perform anisotropy correction when relevant (spmodel has gloglik_anis)
  if (data_object$anisotropy && !orig2optim_object$is_known[["euclid_rotate_is_known"]]) {
    # use rotation parameter from other quadrant
    params_object_q2 <- params_object
    params_object_q2$euclid[["rotate"]] <- pi - params_object_q2$euclid[["rotate"]]

    # find the -2ll products and compute, taking the minimum
    gll_prods_q2 <- gloglik_products(params_object_q2, data_object, estmethod)
    minustwologlik_q2 <- get_minustwologlik(gll_prods_q2, data_object, estmethod)
    minustwologlik <- min(c(minustwologlik, minustwologlik_q2))
  }
  # return -2ll
  minustwologlik
}
