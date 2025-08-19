#' Glance at a fitted model object
#'
#' @description Returns a row of model
#'   summaries from a fitted model object. Glance returns the same number of columns for all models
#'   and estimation methods.
#'
#' @param x A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A single-row tibble with columns
#'   \itemize{
#'     \item \code{n} The sample size.
#'     \item \code{p} The number of fixed effects.
#'     \item \code{npar} The number of estimated covariance parameters.
#'     \item \code{value} The optimized value of the fitting function
#'     \item \code{AIC} The AIC.
#'     \item \code{AICc} The AICc.
#'     \item \code{BIC} The BIC.
#'     \item \code{logLik} The log-likelihood
#'     \item \code{deviance} The deviance.
#'     \item \code{pseudo.r.squared} The pseudo r-squared
#'   }
#'
#' @name glance.SSN2
#' @method glance ssn_lm
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
#' glance(ssn_mod)
glance.ssn_lm <- function(x, ...) {
  is_likbased <- x$estmethod %in% c("ml", "reml")
  tibble::tibble(
    n = x$n,
    p = x$p,
    npar = x$npar,
    value = x$optim$value,
    AIC = ifelse(is_likbased, AIC(x), NA),
    AICc = ifelse(is_likbased, AICc(x), NA),
    BIC = ifelse(is_likbased, BIC(x), NA),
    logLik = ifelse(is_likbased, logLik(x), NA),
    deviance = ifelse(is_likbased, deviance(x), NA),
    pseudo.r.squared = pseudoR2(x),
    # cv.crit = loocv(x)
  )
}

#' @rdname glance.SSN2
#' @method glance ssn_glm
#' @export
glance.ssn_glm <- glance.ssn_lm
