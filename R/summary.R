#' Summarize a fitted model object
#'
#' @description Summarize a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details \code{summary.ssn()} creates a summary of a fitted model object
#'   intended to be printed using \code{print()}. This summary contains
#'   useful information like the original function call, residuals,
#'   a coefficients table, a pseudo r-squared, and estimated covariance
#'   parameters.
#'
#' @return A list with several fitted model quantities used to create
#'   informative summaries when printing.
#'
#' @name summary.SSN2
#' @method summary ssn_lm
#' @export
#'
#' @seealso [print.SSN2]
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
#' summary(ssn_mod)
summary.ssn_lm <- function(object, ...) {
  summary_coefficients_fixed <- data.frame(
    estimates = coef(object, type = "fixed"),
    Std_Error = sqrt(diag(vcov(object, type = "fixed")))
  )

  summary_coefficients_fixed$z_value <- summary_coefficients_fixed$estimates / summary_coefficients_fixed$Std_Error
  summary_coefficients_fixed$p <- 2 * (1 - pnorm(abs(summary_coefficients_fixed$z_value)))

  params_object <- object$coefficients$params_object
  coefficients <- list(fixed = summary_coefficients_fixed, params_object = params_object)
  summary_list <- list(
    call = object$call,
    terms = object$terms,
    residuals = object$residuals,
    coefficients = coefficients,
    pseudoR2 = object$pseudoR2,
    vcov = object$vcov,
    is_known = object$is_known,
    # fn = object$fn,
    anisotropy = object$anisotropy
  )
  new_summary_list <- structure(summary_list, class = paste0("summary.", class(object)))
  new_summary_list
}

#' @name summary.SSN2
#' @method summary ssn_glm
#' @export
summary.ssn_glm <- summary.ssn_lm
