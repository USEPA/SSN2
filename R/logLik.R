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
  if (object$estmethod %in% c("reml", "ml")) {
    minus2loglik <- object$optim$value
    loglik <- -1 / 2 * minus2loglik
    return(loglik)
  } else {
    stop("log likelihood is only defined for the reml or ml estimation")
  }
}

#' @rdname logLik.SSN2
#' @method logLik ssn_glm
#' @export
logLik.ssn_glm <- logLik.ssn_lm
