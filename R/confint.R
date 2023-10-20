#' Confidence intervals for fitted model parameters
#'
#' @description Computes confidence intervals for one or more parameters in a fitted
#'   model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param parm A specification of which parameters are to be given confidence
#'   intervals (a character vector of names). If missing, all parameters are considered.
#' @param level The confidence level required. The default is \code{0.95}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return Gaussian-based confidence intervals (two-sided and equal-tailed) for the
#'   fixed effect coefficients based on the confidence level specified by \code{level}.
#'   For \code{ssn_glm()} objects, confidence intervals are on the link scale.
#'
#' @name confint.SSN2
#' @method confint ssn_lm
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
#' confint(ssn_mod)
#' confint(ssn_mod, level = 0.9)
confint.ssn_lm <- function(object, parm, level = 0.95, ...) {
  # if (type == "fixed") ## may add other confidence intervals later
  alpha <- 1 - level
  # tstar <- qt(1 - alpha / 2, df = object$n - object$p)
  tstar <- qnorm(1 - alpha / 2)
  estimates <- coef(object, type = "fixed")
  variances <- diag(vcov(object, type = "fixed"))
  lower <- estimates - tstar * sqrt(variances)
  upper <- estimates + tstar * sqrt(variances)
  confints <- cbind(lower, upper)
  rownames(confints) <- names(estimates)
  colnames(confints) <- c(paste(alpha / 2 * 100, "%"), paste((1 - alpha / 2) * 100, "%"))
  if (missing(parm)) {
    return(confints)
  } else {
    return(confints[row.names(confints) %in% parm, , drop = FALSE])
  }
}

#' @rdname confint.SSN2
#' @method confint ssn_glm
#' @export
confint.ssn_glm <- confint.ssn_lm
