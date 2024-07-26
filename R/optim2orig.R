#' Transform parameters from optim to original scale
#'
#' @param orig2optim_object An object that contains the parameters on the original scale
#' @param par Current parameter values
#'
#' @noRd
optim2orig <- function(orig2optim_object, par) {
  # fill optim parameter vector with known parameters
  fill_optim_par_val <- fill_optim_par(orig2optim_object, par)
  fill_optim_par_val_ssn <- fill_optim_par_val$par_ssn

  # store all values and perform appropriate inverse transformations
  tailup_de <- exp(fill_optim_par_val_ssn[["tailup_de_log"]])
  tailup_range <- exp(fill_optim_par_val_ssn[["tailup_range_log"]])
  taildown_de <- exp(fill_optim_par_val_ssn[["taildown_de_log"]])
  taildown_range <- exp(fill_optim_par_val_ssn[["taildown_range_log"]])
  euclid_de <- exp(fill_optim_par_val_ssn[["euclid_de_log"]])
  euclid_range <- exp(fill_optim_par_val_ssn[["euclid_range_log"]])
  euclid_rotate <- pi * expit(fill_optim_par_val_ssn[["euclid_rotate_logodds"]])
  euclid_scale <- 1 * expit(fill_optim_par_val_ssn[["euclid_scale_logodds"]])
  nugget <- exp(fill_optim_par_val_ssn[["nugget_log"]])

  # create parameter vector on original scale
  fill_orig_val_ssn <- c(
    tailup_de = tailup_de,
    tailup_range = tailup_range,
    taildown_de = taildown_de,
    taildown_range = taildown_range,
    euclid_de = euclid_de,
    euclid_range = euclid_range,
    euclid_rotate = euclid_rotate,
    euclid_scale = euclid_scale,
    nugget = nugget
  )

  # handle random effects
  if (is.null(fill_optim_par_val$par_randcov)) {
    fill_orig_val_randcov <- NULL
  } else {
    fill_orig_val_randcov <- exp(fill_optim_par_val$par_randcov)
    names(fill_orig_val_randcov) <- gsub("_log", "", names(fill_optim_par_val$par_randcov))
  }

  # return covariance parameters and random effects
  list(orig_ssn = fill_orig_val_ssn, orig_randcov = fill_orig_val_randcov)
}

#' Transform parameters from optim to original scale for glms
#'
#' @param orig2optim_object An object that contains the parameters on the original scale
#' @param par Current parameter values
#'
#' @noRd
optim2orig_glm <- function(orig2optim_object, par) {
  fill_optim_par_val <- fill_optim_par(orig2optim_object, par)

  fill_optim_par_val_ssn <- fill_optim_par_val$par_ssn
  tailup_de <- exp(fill_optim_par_val_ssn[["tailup_de_log"]])
  tailup_range <- exp(fill_optim_par_val_ssn[["tailup_range_log"]])
  taildown_de <- exp(fill_optim_par_val_ssn[["taildown_de_log"]])
  taildown_range <- exp(fill_optim_par_val_ssn[["taildown_range_log"]])
  euclid_de <- exp(fill_optim_par_val_ssn[["euclid_de_log"]])
  euclid_range <- exp(fill_optim_par_val_ssn[["euclid_range_log"]])
  euclid_rotate <- pi * expit(fill_optim_par_val_ssn[["euclid_rotate_logodds"]])
  euclid_scale <- 1 * expit(fill_optim_par_val_ssn[["euclid_scale_logodds"]])
  nugget <- exp(fill_optim_par_val_ssn[["nugget_log"]])
  dispersion <- exp(fill_optim_par_val_ssn[["dispersion_log"]])

  fill_orig_val_ssn <- c(
    tailup_de = tailup_de,
    tailup_range = tailup_range,
    taildown_de = taildown_de,
    taildown_range = taildown_range,
    euclid_de = euclid_de,
    euclid_range = euclid_range,
    euclid_rotate = euclid_rotate,
    euclid_scale = euclid_scale,
    nugget = nugget,
    dispersion = dispersion
  )

  if (is.null(fill_optim_par_val$par_randcov)) {
    fill_orig_val_randcov <- NULL
  } else {
    fill_orig_val_randcov <- exp(fill_optim_par_val$par_randcov)
    names(fill_orig_val_randcov) <- gsub("_log", "", names(fill_optim_par_val$par_randcov))
  }
  list(orig_ssn = fill_orig_val_ssn, orig_randcov = fill_orig_val_randcov)
}
