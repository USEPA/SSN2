#' Calculate variance-covariance matrix for a fitted model object
#'
#' @description Calculate variance-covariance matrix for a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param type For \code{type = "fixed"} (the default), the variance-covariance matrix
#'   of the fixed effects. If Satterthwaite degrees of freedom were calculated,
#'   \code{type = "cov"} returns the variance-covariance matrix of the covariance
#'   parameters. In SSN2, \code{type = "ssn"} returns the tailup, taildown,
#'   Euclidean, and nugget subset; \code{"tailup"}, \code{"taildown"},
#'   \code{"euclid"}, and \code{"nugget"} return their individual subsets.
#'   \code{type = "randcov"} returns the variance-covariance matrix of just
#'   the random effects (if relevant).
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return The variance-covariance matrix of fixed effect
#'   coefficients obtained via \code{coef(..., type = "fixed")}
#'   or the variance-covariance matrix of  estimated covariance parameters (when available).
#'
#' @name vcov.SSN2
#' @method vcov ssn_lm
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
#' vcov(ssn_mod)
vcov.ssn_lm <- function(object, type = "fixed", ...) {
  # "cov"/"ssn"/"tailup"/"taildown"/"euclid"/"nugget"/"randcov" are only
  # available when Satterthwaite information was cached at fit time (see
  # ssn_lm()'s ddf argument) -- no on-demand computation and no method
  # argument, matching spmodel's vcov.splm(), which purely returns
  # object$vcov$cov/spcov/randcov
  if (identical(type, "fixed")) {
    return(object$vcov$fixed)
  }

  valid_types <- c("cov", "ssn", "tailup", "taildown", "euclid", "nugget", "randcov")
  if (!type %in% valid_types) {
    stop("type must be \"fixed\", \"cov\", \"ssn\", \"tailup\", \"taildown\", \"euclid\", \"nugget\", or \"randcov\".", call. = FALSE)
  }

  vcov_theta <- object$vcov$cov
  if (is.null(vcov_theta)) {
    return(NULL)
  }

  if (identical(type, "cov")) {
    return(vcov_theta)
  }

  cov_names_free <- rownames(vcov_theta)
  spcov_names_free <- intersect(get_spcov_field_names(), cov_names_free)

  if (identical(type, "ssn")) {
    return(vcov_theta[spcov_names_free, spcov_names_free, drop = FALSE])
  }

  if (type %in% c("tailup", "taildown", "euclid", "nugget")) {
    component_names_free <- cov_names_free[startsWith(cov_names_free, type)]
    if (length(component_names_free) == 0) {
      return(NULL)
    }
    return(vcov_theta[component_names_free, component_names_free, drop = FALSE])
  }

  # type == "randcov"
  randcov_names_free <- setdiff(cov_names_free, spcov_names_free)
  if (length(randcov_names_free) == 0) {
    return(NULL)
  }
  vcov_theta[randcov_names_free, randcov_names_free, drop = FALSE]
}

#' @param var_correct A logical indicating whether to return the corrected variance-covariance
#'   matrix for models fit using [ssn_glm()] (when \code{family} is different
#'   from \code{"Gaussian"}). The default is \code{TRUE}.
#' @rdname vcov.SSN2
#' @method vcov ssn_glm
#' @export
vcov.ssn_glm <- function(object, var_correct = TRUE, ...) {
  # type is hard-coded (see vcov.ssn_lm) since only fixed effects are supported
  type <- "fixed"
  if (type == "fixed") {
    # "corrected" adjusts the naive fixed-effect covariance to account for the
    # extra uncertainty from estimating the latent random effects; "uncorrected"
    # is that naive (asymptotic) estimate
    if (var_correct) {
      return(object$vcov$fixed$corrected)
    } else {
      return(object$vcov$fixed$uncorrected)
    }
  }
}
