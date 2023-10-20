#' Extract fitted model coefficients
#'
#' @description \code{coef} extracts fitted model coefficients from fitted model objects.
#'   \code{coefficients} is an alias for it.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param type \code{"fixed"} for fixed effect coefficients, \code{"tailup"} for
#'   tailup covariance parameter coefficients, \code{"taildown"} for
#'   taildown covariance parameter coefficients, \code{"euclid"} for
#'   Euclidean covariance parameter coefficients, \code{"nugget"} for
#'   nugget covariance parameter coefficients, \code{"dispersion"} for
#'   the dispersion parameter coefficient (\code{ssn_glm()} objects), \code{"randcov"} for random effect
#'   variance coefficients, or \code{"ssn"} for all of the tailup, taildown,
#'   Euclidean, nugget, and dispersion (\code{ssn_glm()} objects) parameter coefficients.
#'   Defaults to \code{"fixed"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A named vector of coefficients.
#'
#' @name coef.SSN2
#' @method coef ssn_lm
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
#' coef(ssn_mod)
#' coef(ssn_mod, type = "tailup")
#' coefficients(ssn_mod)
coef.ssn_lm <- function(object, type = "fixed", ...) {
  if (type == "fixed") {
    return(object$coefficients$fixed)
  } else if (type == "ssn") {
    randcov_index <- which(names(object$coefficients$params_object) == "randcov")
    return(object$coefficients$params_object[-randcov_index])
  } else if (type == "tailup") {
    return(object$coefficients$params_object$tailup)
  } else if (type == "taildown") {
    return(object$coefficients$params_object$taildown)
  } else if (type == "euclid") {
    return(object$coefficients$params_object$euclid)
  } else if (type == "nugget") {
    return(object$coefficients$params_object$nugget)
  } else if (type == "randcov") {
    return(object$coefficients$params_object$randcov)
  } else {
    stop("Invalid type argument. The type argument must be \"fixed\", \"ssn\", \"tailup\",  \"taildown\",  \"euclid\",  \"nugget\", or \"randcov\".", call. = FALSE)
  }
}
#' @rdname coef.SSN2
#' @export
coefficients.ssn_lm <- coef.ssn_lm

#' @rdname coef.SSN2
#' @method coef ssn_glm
#' @export
coef.ssn_glm <- function(object, type = "fixed", ...) {
  if (type == "fixed") {
    return(object$coefficients$fixed)
  } else if (type == "ssn") {
    randcov_index <- which(names(object$coefficients$params_object) == "randcov")
    return(object$coefficients$params_object[-randcov_index])
  } else if (type == "tailup") {
    return(object$coefficients$params_object$tailup)
  } else if (type == "taildown") {
    return(object$coefficients$params_object$taildown)
  } else if (type == "euclid") {
    return(object$coefficients$params_object$euclid)
  } else if (type == "nugget") {
    return(object$coefficients$params_object$nugget)
  } else if (type == "dispersion") {
    return(object$coefficients$params_object$dispersion)
  } else if (type == "randcov") {
    return(object$coefficients$params_object$randcov)
  } else {
    stop("Invalid type argument. The type argument must be \"fixed\", \"ssn\", \"tailup\",  \"taildown\",  \"euclid\",  \"nugget\", \"dispersion\", or \"randcov\".", call. = FALSE)
  }
}

#' @rdname coef.SSN2
#' @export
coefficients.ssn_glm <- coef.ssn_glm
