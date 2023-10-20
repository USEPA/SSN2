#' Fitted model deviance
#'
#' @description Returns the deviance of a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details The deviance is twice the difference in log-likelihoods between the
#'   saturated (perfect-fit) model and the fitted model.
#'
#' @return The deviance.
#'
#' @name deviance.SSN2
#' @method deviance ssn_lm
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
#' ssn_mod <- ssn_glm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   family = "Gamma",
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' deviance(ssn_mod)
deviance.ssn_lm <- function(object, ...) {
  if (object$estmethod %in% c("reml", "ml")) {
    deviance <- object$deviance
    return(deviance)
  } else {
    stop("deviance is only defined for reml or ml estimation")
  }
}

#' @rdname deviance.SSN2
#' @method deviance ssn_glm
#' @export
deviance.ssn_glm <- deviance.ssn_lm
