#' Warn if \code{optim()} did not converge
#'
#' @param convergence The convergence code returned by \code{optim()}
#'   (\code{0} indicates success), or \code{NA} when \code{optim()} was not
#'   called.
#'
#' @return \code{invisible(NULL)}, after emitting a warning if
#'   \code{convergence} is nonzero.
#'
#' @noRd
warn_optim_convergence <- function(convergence) {
  if (is.na(convergence) || convergence == 0) {
    return(invisible())
  }
  warning(
    paste0(
      "optim() did not converge while fitting the model ",
      "(convergence code ", convergence, "). Fitted model estimates may be unreliable."
    ),
    call. = FALSE
  )
  invisible()
}

#' Warn if fitted spatial variance parameters are near a numerical boundary of zero
#'
#' @param params_object A joint covariance parameter object.
#' @param diagtol A numerical-stability tolerance added to the covariance
#'   diagonal.
#'
#' @return \code{invisible(NULL)}, after emitting a warning if the total
#'   spatially dependent plus nugget variance is within 10 times
#'   \code{diagtol} of zero.
#'
#' @noRd
warn_spcov_boundary <- function(params_object, diagtol) {
  if (diagtol <= 0) {
    return(invisible())
  }
  total_de <- params_object$tailup[["de"]] + params_object$taildown[["de"]] +
    params_object$euclid[["de"]] + params_object$nugget[["nugget"]]
  if (total_de <= 10 * diagtol) {
    warning(
      "The fitted spatial variance parameters are near a numerical boundary of zero, so the likelihood value may be unreliable for model comparisons (e.g., AIC(), AICc(), BIC()).",
      call. = FALSE
    )
  }
  invisible()
}

#' Warn if fitted binomial probabilities indicate perfect separation
#'
#' @param fitted_response The fitted response-scale values.
#' @param family The GLM family name; only \code{"binomial"} is checked.
#'
#' @return \code{invisible(NULL)}, after emitting a warning if at least 99%
#'   of fitted probabilities are numerically 0 or 1.
#'
#' @noRd
warn_fitted_saturation <- function(fitted_response, family) {
  if (family != "binomial") {
    return(invisible())
  }
  tol <- 1e-6
  saturated <- fitted_response < tol | fitted_response > 1 - tol
  if (mean(saturated) >= 0.99) {
    warning(
      "Nearly all fitted probabilities are numerically 0 or 1. Perfect separation detected.",
      call. = FALSE
    )
  }
  invisible()
}
