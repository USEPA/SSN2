#' Variability component comparison
#'
#' @description Compare the proportion of total variability explained by the fixed effects
#'   and each variance parameter.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A tibble that partitions the the total variability by the fixed effects
#'   and each variance parameter. The proportion of variability explained by the
#'   fixed effects is the pseudo R-squared obtained by \code{psuedoR2()}. The
#'   remaining proportion is spread accordingly among each variance parameter:
#'   \code{"tailup_de"}, \code{"taildown_de"}, \code{"euclid_de"}, \code{"nugget"},
#'   and if random effects are used, each named random effect. For \code{ssn_glm()},
#'   models, only the variances on the link scale are considered (i.e., the variance
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
varcomp.ssn_glm <- varcomp.ssn_lm
