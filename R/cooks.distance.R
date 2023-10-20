#' Compute Cook's distance
#'
#' @description Compute the Cook's distance for each observation from a fitted
#'   model object.
#'
#' @param model A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Cook's distance measures the influence of an observation on a fitted
#'   model object. If an observation is influential, its omission from the data
#'   noticeably impacts parameter estimates. The larger the Cook's distance, the
#'   larger the influence.
#'
#' @return A vector of Cook's distance values for each observation from the
#'   fitted model object.
#'
#' @name cooks.distance.SSN2
#' @method cooks.distance ssn_lm
#' @export
#'
#' @seealso [augment.SSN2()] [hatvalues.SSN2()] [influence.SSN2()] [residuals.SSN2()]
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
#' cooks.distance(ssn_mod)
cooks.distance.ssn_lm <- function(model, ...) {
  model$cooks_distance
}

#' @rdname cooks.distance.SSN2
#' @method cooks.distance ssn_glm
#' @export
cooks.distance.ssn_glm <- cooks.distance.ssn_lm
