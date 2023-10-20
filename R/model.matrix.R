#' Extract the model matrix from a fitted model object
#'
#' @description Extract the model matrix (X) from a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The model matrix (of the fixed effects), whose rows represent
#'   observations and whose columns represent explanatory variables corresponding
#'   to each fixed effect.
#'
#' @name model.matrix.SSN2
#' @method model.matrix ssn_lm
#' @export
#'
#' @seealso [stats::model.matrix()]
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
#' model.matrix(ssn_mod)
model.matrix.ssn_lm <- function(object, ...) {
  # model.matrix(formula(object, ...), model.frame(object, ...), ...) too much customization
  model.matrix(object$formula, model.frame(object), contrasts = object$contrasts)
}

#' @rdname model.matrix.SSN2
#' @method model.matrix ssn_glm
#' @export
model.matrix.ssn_glm <- model.matrix.ssn_lm
