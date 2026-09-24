#' Variability component comparison
#'
#' @description Compare the proportion of total variability explained by the fixed effects
#'   and each variance parameter.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details The total variability in the response is decomposed into a
#'   portion explained by the fixed effects and a portion explained by each
#'   variance parameter in the fitted covariance structure:
#'   \itemize{
#'     \item \code{tailup_de}: the tailup (downstream-moving) random error variance.
#'     \item \code{taildown_de}: the taildown (upstream-moving) random error variance.
#'     \item \code{euclid_de}: the Euclidean random error variance.
#'     \item \code{nugget}: the nugget (spatially independent) random error variance.
#'     \item random effects: if \code{object} was fit with a \code{random}
#'       argument, one additional variance parameter per named random effect
#'       term (e.g., a random intercept's grouping variable), representing
#'       the variance attributable to that grouping.
#'   }
#'   The proportion of variability
#'   explained by the fixed effects is the pseudo R-squared returned by
#'   [pseudoR2()]. The remaining
#'   \code{1 - pseudoR2} proportion is then split among \code{tailup_de},
#'   \code{taildown_de}, \code{euclid_de}, \code{nugget}, and any random
#'   effect variances, in proportion to their share of the total variance
#'   (the sum of \code{tailup_de}, \code{taildown_de}, \code{euclid_de},
#'   \code{nugget}, and all random effect variances).
#'
#' @return A tibble that partitions the total variability by the fixed effects
#'   and each variance parameter (see Details). For \code{ssn_glm()} models,
#'   only the variances on the link scale are considered (i.e., the variance
#'   function of the response is omitted).
#'
#' @name varcomp.SSN2
#' @method varcomp ssn_lm
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
#' varcomp(ssn_mod)
varcomp.ssn_lm <- function(object, ...) {
  PR2 <- pseudoR2(object)
  tailup_de <- object$coefficients$params_object$tailup[["de"]]
  taildown_de <- object$coefficients$params_object$taildown[["de"]]
  euclid_de <- object$coefficients$params_object$euclid[["de"]]
  nugget <- object$coefficients$params_object$nugget[["nugget"]]
  randcov <- as.vector(object$coefficients$params_object$randcov)
  total_var <- sum(tailup_de, taildown_de, euclid_de, nugget, randcov)
  varcomp_names <- c("Covariates (PR-sq)", "tailup_de", "taildown_de", "euclid_de", "nugget", c(names(object$coefficients$params_object$randcov)))
  varcomp_values <- c(PR2, (1 - PR2) * c(tailup_de, taildown_de, euclid_de, nugget, randcov) / total_var)
  tibble::tibble(varcomp = varcomp_names, proportion = varcomp_values)
}

#' @rdname varcomp.SSN2
#' @method varcomp ssn_glm
#' @export
# ssn_glm() reuses the ssn_lm() varcomp logic directly (same variance
# decomposition applies on the link-function scale)
varcomp.ssn_glm <- varcomp.ssn_lm
