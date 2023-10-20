#' Compute a pseudo r-squared
#'
#' @description Compute a pseudo r-squared for a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param adjust A logical indicating whether the pseudo r-squared
#'   should be adjusted to account for the number of explanatory variables. The
#'   default is \code{FALSE}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Several pseudo r-squared statistics exist for in the literature.
#'   We define this pseudo r-squared as one minus the ratio of the deviance of a full model
#'   relative to the deviance of a null (intercept only) model. This pseudo r-squared
#'   can be viewed as a generalization of the classical r-squared definition
#'   seen as one minus the ratio of error sums of squares from the full model relative
#'   to the error sums of squares from the null model. If adjusted, the adjustment
#'   is analogous to the the classical r-squared adjustment.
#'
#' @return The pseudo r-squared as a numeric vector.
#'
#' @name pseudoR2.SSN2
#' @method pseudoR2 ssn_lm
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
#' pseudoR2(ssn_mod)
pseudoR2.ssn_lm <- function(object, adjust = FALSE, ...) {
  if (adjust) {
    has_intercept <- "(Intercept)" %in% tidy(object)$term
    pr2 <- object$pseudoR2
    pr2_adj <- 1 - (1 - pr2) * (object$n - 1 * has_intercept) / (object$n - object$p)
    return(pr2_adj)
  } else {
    return(object$pseudoR2)
  }
}

#' @rdname pseudoR2.SSN2
#' @method pseudoR2 ssn_glm
#' @export
pseudoR2.ssn_glm <- pseudoR2.ssn_lm
