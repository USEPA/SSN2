#' Compute leverage (hat) values
#'
#' @description Compute the leverage (hat) value for each observation from a fitted
#'   model object.
#'
#' @param model A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Leverage values measure how far an observation's explanatory variables
#'   are relative to the average of the explanatory variables. In other words, observations with high
#'   leverage are typically considered to have an extreme or unusual combination of explanatory
#'   variables. Leverage values are the diagonal of the hat (projection) matrix.
#'   The larger the hat value, the larger the leverage.
#'
#' @return A vector of leverage (hat) values for each observation from the
#'   fitted model object.
#'
#' @name hatvalues.SSN2
#' @method hatvalues ssn_lm
#' @export
#'
#' @seealso [augment.SSN2()] [cooks.distance.SSN2()] [influence.SSN2()] [residuals.SSN2()]
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
#' hatvalues(ssn_mod)
hatvalues.ssn_lm <- function(model, ...) {
  model$hatvalues
}

#' @rdname hatvalues.SSN2
#' @method hatvalues ssn_glm
#' @export
hatvalues.ssn_glm <- hatvalues.ssn_lm
