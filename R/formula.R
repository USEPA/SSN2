#' Model formulae
#'
#' Return formula used by a fitted model object.
#'
#' @param x A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The formula used by a fitted model object.
#'
#' @name formula.SSN2
#' @method formula ssn_lm
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
#' formula(ssn_mod)
formula.ssn_lm <- function(x, ...) {
  formula(x$formula)
}

#' @rdname formula.SSN2
#' @method formula ssn_glm
#' @export
formula.ssn_glm <- formula.ssn_lm
