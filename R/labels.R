#' Find labels from object
#'
#' @description Find a suitable set of labels from a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A character vector containing the terms used for the fixed effects
#'   from a fitted model object.
#'
#' @name labels.SSN2
#' @method labels ssn_lm
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
#' labels(ssn_mod)
labels.ssn_lm <- function(object, ...) {
  labels(terms(formula(object)))
}

#' @rdname labels.SSN2
#' @method labels ssn_glm
#' @export
labels.ssn_glm <- labels.ssn_lm
