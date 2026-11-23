#' Transform parameters from original to optim scale
#'
#' @param initial_object Initial value object
#' @param data_object Data object
#' @noRd
orig2optim <- function(initial_object, data_object) {
  # find parameters on optim (transformed) scale
  ssn <- orig2optim_ssn_components(initial_object, data_object)
  ssn$value <- clamp_optim_scale(ssn$value, ssn$is_known)
  n_est_ssn <- sum(!ssn$is_known)

  # handle random effects
  rc <- orig2optim_randcov_components(initial_object)

  # return the orig2optim initial value (and other) information
  orig2optim_list <- list(
    value = ssn$value,
    is_known = ssn$is_known,
    n_est_ssn = n_est_ssn,
    randcov_value = rc$value,
    randcov_is_known = rc$is_known,
    n_est_rand = rc$n_est,
    n_est = n_est_ssn + rc$n_est,
    range_constrain_value = data_object$range_constrain_value,
    classes = c(
      tailup = class(initial_object$tailup_initial),
      taildown = class(initial_object$taildown_initial),
      euclid = class(initial_object$euclid_initial),
      nugget = class(initial_object$nugget_initial)
    )
  )
}

orig2optim_glm <- function(initial_object, data_object) {
  ssn <- orig2optim_ssn_components(initial_object, data_object)

  # dispersion (GLM-only component)
  dispersion <- initial_object$dispersion_initial$initial[["dispersion"]]
  dispersion_log <- log(dispersion)
  dispersion_is_known <- initial_object$dispersion_initial$is_known[["dispersion"]]
  ssn$value <- c(ssn$value, dispersion_log = dispersion_log)
  ssn$is_known <- c(ssn$is_known, dispersion_is_known = dispersion_is_known)

  ssn$value <- clamp_optim_scale(ssn$value, ssn$is_known)
  n_est_ssn <- sum(!ssn$is_known)

  # random effects
  rc <- orig2optim_randcov_components(initial_object)

  orig2optim_list <- list(
    value = ssn$value,
    is_known = ssn$is_known,
    n_est_ssn = n_est_ssn,
    randcov_value = rc$value,
    randcov_is_known = rc$is_known,
    n_est_rand = rc$n_est,
    n_est = n_est_ssn + rc$n_est,
    range_constrain_value = data_object$range_constrain_value,
    classes = c(
      tailup = class(initial_object$tailup_initial),
      taildown = class(initial_object$taildown_initial),
      euclid = class(initial_object$euclid_initial),
      nugget = class(initial_object$nugget_initial),
      dispersion = class(initial_object$dispersion_initial)
    )
  )
}
