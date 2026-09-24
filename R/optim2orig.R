#' Transform parameters from optim to original scale
#'
#' @param orig2optim_object An object that contains the parameters on the original scale
#' @param par Current parameter values
#'
#' @noRd
optim2orig <- function(orig2optim_object, par) {
  # fill optim parameter vector with known parameters
  fill_optim_par_val <- fill_optim_par(orig2optim_object, par)

  # store all values and perform appropriate inverse transformations
  euclid_type <- remove_covtype(orig2optim_object$classes[["euclid"]])
  fill_orig_val_ssn <- optim2orig_ssn_components(
    fill_optim_par_val$par_ssn, euclid_type, orig2optim_object$range_constrain_value
  )

  # handle random effects
  fill_orig_val_randcov <- optim2orig_randcov_components(fill_optim_par_val$par_randcov)

  # return covariance parameters and random effects
  list(orig_ssn = fill_orig_val_ssn, orig_randcov = fill_orig_val_randcov)
}

optim2orig_glm <- function(orig2optim_object, par) {
  fill_optim_par_val <- fill_optim_par(orig2optim_object, par)

  euclid_type <- remove_covtype(orig2optim_object$classes[["euclid"]])
  fill_orig_val_ssn <- optim2orig_ssn_components(
    fill_optim_par_val$par_ssn, euclid_type, orig2optim_object$range_constrain_value
  )
  dispersion <- exp(fill_optim_par_val$par_ssn[["dispersion_log"]])
  fill_orig_val_ssn <- c(fill_orig_val_ssn, dispersion = dispersion)

  fill_orig_val_randcov <- optim2orig_randcov_components(fill_optim_par_val$par_randcov)

  list(orig_ssn = fill_orig_val_ssn, orig_randcov = fill_orig_val_randcov)
}
