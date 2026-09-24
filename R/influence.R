#' Regression diagnostics
#'
#' @description Provides basic quantities which are used in forming
#'   a wide variety of diagnostics for checking the quality of fitted model objects.
#'
#' @param model A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details This function calls [residuals.SSN2()], [hatvalues.SSN2()],
#'   and [cooks.distance.SSN2()] and puts the results into a tibble. It is
#'   primarily used when calling [augment.SSN2()].
#'
#' @return A tibble with residuals (\code{.resid}), leverage values (\code{.hat}),
#'   cook's distance (\code{.cooksd}), and standardized residuals (\code{.std.resid}).
#'
#' @name influence.SSN2
#' @method influence ssn_lm
#' @export
#'
#' @seealso [augment.SSN2()] [cooks.distance.SSN2()] [hatvalues.SSN2()] [residuals.SSN2()]
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
#' influence(ssn_mod)
influence.ssn_lm <- function(model, ...) {
  # standardized residuals are requested via a separate residuals() call
  # (rather than derived here from .resid and .hat) since the standardization
  # accounts for the full model covariance, not just leverage
  tibble::tibble(
    .resid = residuals(model),
    .hat = hatvalues(model),
    .cooksd = cooks.distance(model),
    .std.resid = residuals(model, type = "standardized")
  )
}

#' @rdname influence.SSN2
#' @method influence ssn_glm
#' @export
influence.ssn_glm <- function(model, ...) {
  tibble::tibble(
    # allow ... so type can be passed to residuals()
    .resid = residuals(model, ...),
    .hat = hatvalues(model),
    .cooksd = cooks.distance(model),
    # standardized residuals always use their own default type regardless of
    # any type passed above, since standardization needs a specific residual definition
    .std.resid = residuals(model, type = "standardized")
  )
}
