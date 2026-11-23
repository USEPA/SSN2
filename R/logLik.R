#' Extract log-likelihood
#'
#' @description Find the log-likelihood of a fitted model.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The log-likelihood.
#'
#' @name logLik.SSN2
#' @method logLik ssn_lm
#' @export
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' logLik(ssn_mod)
logLik.ssn_lm <- function(object, ...) {
  # only defined for likelihood-based fits; other estmethods don't optimize a
  # likelihood so there is no minus2loglik to convert
  if (object$estmethod %in% c("reml", "ml")) {
    # optim minimizes -2*loglik, so recover loglik by undoing that transform
    minus2loglik <- object$optim$value
    loglik <- -1 / 2 * minus2loglik
    # number of estimated parameters
    # reml profiles out the fixed effects, so its likelihood is a function of
    # only the npar covariance parameters; ml's likelihood also depends on
    # the p fixed effects, so those are counted too
    if (object$estmethod == "ml") {
      n_est_param <- object$npar + object$p
    } else {
      n_est_param <- object$npar
    }
    # lm also returns nall which does help with reml warnings of models
    # with different fixed effect structures
    # still insufficient compared to previous version
    loglik <- structure(loglik, nobs = object$n, df = n_est_param, class = "logLik")
    return(loglik)
  } else {
    stop("logLik is only defined if estmethod is \"ml\" or \"reml\".", call. = FALSE)
  }
}

#' @rdname logLik.SSN2
#' @method logLik ssn_glm
#' @export
# the log-likelihood is extracted from object$optim$value the same way
# regardless of family, so the GLM method reuses the Gaussian implementation
logLik.ssn_glm <- logLik.ssn_lm
