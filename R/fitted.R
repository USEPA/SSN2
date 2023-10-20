#' Extract model fitted values
#'
#' @description Extract fitted values from fitted model objects. \code{fitted.values}
#'   is an alias.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param type \code{"response"} for fitted values of the response,
#'   \code{"tailup"} for fitted values of the tailup random errors,
#'   \code{"taildown"} for fitted values of the taildown random errors,
#'   \code{"euclid"} for fitted values of the Euclidean random errors,
#'   \code{"nugget"} for fitted values of the nugget random errors,
#'   or \code{"randcov"} for fitted values of the random effects. If from
#'   [ssn_glm()], \code{"link"} for fitted values on the link scale.
#'   The default is \code{"response"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details When \code{type} is \code{"response"}, the fitted values
#'   for each observation are the standard fitted values \eqn{X \hat{\beta}}.
#'   When \code{type} is \code{"tailup"}, \code{"taildown"}, \code{"euclid"},
#'   or \code{"nugget"} the fitted values for each observation
#'   are (generally) the best linear unbiased predictors of the respective random error.
#'   When \code{type} is \code{"randcov"}, the fitted
#'   values for each level of each random effect are (generally) the best linear unbiased
#'   predictors of the corresponding random effect. The fitted values for \code{type}
#'   \code{"tailup"}, \code{"taildown"}, \code{"euclid"},
#'   \code{"nugget"}, and \code{"randcov"} can generally be used to check assumptions
#'   for each component of the fitted model object (e.g., check a Gaussian assumption).
#'
#'   If from [ssn_glm()], when \code{type} is \code{"response"}, the fitted values
#'   for each observation are the standard fitted values on the inverse link
#'   scale: \eqn{g^{-1}}(\eqn{X \hat{\beta} + \nu}), where \eqn{g(.)} is a link function,
#'   \eqn{\beta} are the fixed effects, and \eqn{\nu} are the spatial and random effects.
#'
#' @return The fitted values according to \code{type}.
#'
#' @name fitted.SSN2
#' @method fitted ssn_lm
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
#' fitted(ssn_mod)
#' fitted.values(ssn_mod)
fitted.ssn_lm <- function(object, type = "response", ...) {
  if (type == "response") {
    fitted_val <- object$fitted$response
  } else if (type == "tailup") {
    fitted_val <- object$fitted$tailup
  } else if (type == "taildown") {
    fitted_val <- object$fitted$taildown
  } else if (type == "euclid") {
    fitted_val <- object$fitted$euclid
  } else if (type == "nugget") {
    fitted_val <- object$fitted$nugget
  } else if (type == "randcov") {
    fitted_val <- object$fitted$randcov
  } else {
    stop("Invalid type argument. The type argument must be \"response\", \"tailup\",  \"taildown\",  \"euclid\",  \"nugget\", or \"randcov\".", call. = FALSE)
  }
  fitted_val
}
#' @rdname fitted.SSN2
#' @method fitted.values ssn_lm
#' @export
fitted.values.ssn_lm <- fitted.ssn_lm

#' @rdname fitted.SSN2
#' @method fitted ssn_glm
#' @export
fitted.ssn_glm <- function(object, type = "response", ...) {
  if (type == "response") {
    fitted_val <- object$fitted$response
  } else if (type == "link") {
    fitted_val <- object$fitted$link
  } else if (type == "tailup") {
    fitted_val <- object$fitted$tailup
  } else if (type == "taildown") {
    fitted_val <- object$fitted$taildown
  } else if (type == "euclid") {
    fitted_val <- object$fitted$euclid
  } else if (type == "nugget") {
    fitted_val <- object$fitted$nugget
  } else if (type == "randcov") {
    fitted_val <- object$fitted$randcov
  } else {
    stop("Invalid type argument. The type argument must be \"response\", \"tailup\",  \"taildown\",  \"euclid\",  \"nugget\", or \"randcov\".", call. = FALSE)
  }
  fitted_val
}
#' @rdname fitted.SSN2
#' @method fitted.values ssn_glm
#' @export
fitted.values.ssn_glm <- fitted.ssn_glm
