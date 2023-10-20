#' Extract the model frame from a fitted model object
#'
#' @description Extract the model frame from a fitted model object.
#'
#' @param formula A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A model frame that contains the variables used by the formula
#'   for the fitted model object.
#'
#' @name model.frame.SSN2
#' @method model.frame ssn_lm
#' @export
#'
#' @seealso [stats::model.frame()]
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
#' model.frame(ssn_mod)
model.frame.ssn_lm <- function(formula, ...) {
  # model.frame(formula(formula, ...), data = formula$data, ...) too much customization
  model.frame(formula(formula),
    data = sf::st_drop_geometry(formula$ssn.object$obs),
    drop.unused.levels = TRUE, na.action = na.omit
  )
}

#' @rdname model.frame.SSN2
#' @method model.frame ssn_glm
#' @export
model.frame.ssn_glm <- model.frame.ssn_lm
