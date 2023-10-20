#' Extract fitted model residuals
#'
#' @description Extract residuals from a fitted model object.
#'   \code{resid} is an alias.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param type \code{"response"} for response residuals, \code{"pearson"}
#'   for Pearson residuals, or \code{"standardized"} for standardized residuals.
#'   For \code{ssn_lm()} fitted model objects, the default is \code{"response"}.
#'   For \code{ssn_glm()} fitted model objects, deviance residuals are also
#'   available (\code{"deviance"}) and are the default residual type.
#' @param ... Other arguments. Not used (needed for generic consistency).
#' @param model A fitted model object from [ssn_lm()] or [ssn_glm()].
#'
#' @details The response residuals are taken as the response minus the fitted values
#'   for the response: \eqn{y - X \hat{\beta}}. The Pearson residuals are the
#'   response residuals pre-multiplied by their inverse square root.
#'   The standardized residuals are Pearson residuals divided by the square
#'   root of one minus the leverage (hat) value. The standardized residuals are often used to
#'   check model assumptions, as they have mean zero and variance approximately one.
#'
#'   \code{rstandard()} is an alias for \code{residuals(model, type = "standardized")}.
#'
#' @return The residuals as a numeric vector.
#'
#' @name residuals.SSN2
#' @method residuals ssn_lm
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
#' residuals(ssn_mod)
#' resid(ssn_mod)
#' rstandard(ssn_mod)
residuals.ssn_lm <- function(object, type = "response", ...) {
  if (type == "response") {
    return(object$residuals$response)
  } else if (type == "pearson") {
    return(object$residuals$pearson)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("residuals must be response or pearson or standardized")
  }
}
#' @rdname residuals.SSN2
#' @export
resid.ssn_lm <- residuals.ssn_lm

#' @rdname residuals.SSN2
#' @method rstandard ssn_lm
#' @export
rstandard.ssn_lm <- function(model, ...) {
  residuals.ssn_lm(model, type = "standardized")
}



#' @rdname residuals.SSN2
#' @method residuals ssn_glm
#' @export
residuals.ssn_glm <- function(object, type = "deviance", ...) {
  if (type == "response") {
    return(object$residuals$response)
  } else if (type == "deviance") {
    return(object$residuals$deviance)
  } else if (type == "pearson") {
    return(object$residuals$pearson)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("residuals must be deviance or response or pearson or standardized")
  }
}

#' @rdname residuals.SSN2
#' @export
resid.ssn_glm <- residuals.ssn_glm

#' @rdname residuals.SSN2
#' @method rstandard ssn_glm
#' @export
rstandard.ssn_glm <- function(model, ...) {
  residuals.ssn_glm(model, type = "standardized")
}
