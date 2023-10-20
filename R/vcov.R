#' Calculate variance-covariance matrix for a fitted model object
#'
#' @description Calculate variance-covariance matrix for a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The variance-covariance matrix of coefficients obtained via \code{coef()}.
#'   Currently, only the variance-covariance matrix of the fixed effects is supported.
#'
#' @name vcov.SSN2
#' @method vcov ssn_lm
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
#' vcov(ssn_mod)
vcov.ssn_lm <- function(object, ...) {
  type <- "fixed"
  if (type == "fixed") {
    return(object$vcov$fixed)
  }
}

#' @param var_correct A logical indicating whether to return the corrected variance-covariance
#'   matrix for models fit using [ssn_glm()] (when \code{family} is different
#'   from \code{"Gaussian"}). The default is \code{TRUE}.
#' @rdname vcov.SSN2
#' @method vcov ssn_glm
#' @export
vcov.ssn_glm <- function(object, var_correct = TRUE, ...) {
  type <- "fixed"
  if (type == "fixed") {
    if (var_correct) {
      return(object$vcov$fixed$corrected)
    } else {
      return(object$vcov$fixed$uncorrected)
    }
  }
}
