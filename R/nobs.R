#' Extract the number of observations
#'
#' @description Find the number of observations used for model fitting. This
#'   exists so that generic R functions that call \code{nobs()} (e.g.,
#'   \code{stats::nobs()}) work on \code{SSN2} fitted model objects.
#'   Could alternatively rename object$n as object$nobs and rely on stats::nobs.default.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The number of observations used for modeling (i.e., \code{object$n}).
#'
#' @name nobs.SSN2
#' @method nobs ssn_lm
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
#' nobs(ssn_mod)
nobs.ssn_lm <- function(object, ...) {
  object$n
}

#' @rdname nobs.SSN2
#' @method nobs ssn_glm
#' @export
# ssn_glm objects store n the same way as ssn_lm objects, so reuse that method
nobs.ssn_glm <- nobs.ssn_lm
