#' Back-transform SSN spatial covariance parameters from the optimization scale
#'
#' Inverts the log/logit-odds transforms used during optimization: variance
#' and range parameters are exponentiated, and Euclidean rotate/scale/shape
#' parameters (bounded fields) are inverse-logit-transformed back onto their
#' original (0, pi]/(0, 1]/family-specific ranges.
#'
#' @param fill_optim_par_val_ssn A named vector of optimization-scale spatial
#'   covariance parameters.
#' @param euclid_type The Euclidean covariance type, required only when
#'   \code{euclid_extra} is on the logit-odds scale (\code{"matern"} or
#'   \code{"pexponential"}).
#' @param range_constrain_value The shared range bound used by any range
#'   parameter on the logit-odds scale (see
#'   \code{\link{get_range_constrain_setup}()}); ignored/unused when every
#'   range parameter is on the (unconstrained) log scale.
#'
#' @return A named vector of original-scale spatial covariance parameters
#'   (\code{tailup_de}, \code{tailup_range}, \code{taildown_de},
#'   \code{taildown_range}, \code{euclid_de}, \code{euclid_range},
#'   optionally \code{euclid_extra}, \code{euclid_rotate}, \code{euclid_scale},
#'   \code{nugget}).
#'
#' @noRd
optim2orig_ssn_components <- function(fill_optim_par_val_ssn, euclid_type = NULL, range_constrain_value = NULL) {
  euclid_values <- c(
    tailup_de = exp(fill_optim_par_val_ssn[["tailup_de_log"]]),
    tailup_range = optim2orig_range_component(fill_optim_par_val_ssn, range_constrain_value, "tailup"),
    taildown_de = exp(fill_optim_par_val_ssn[["taildown_de_log"]]),
    taildown_range = optim2orig_range_component(fill_optim_par_val_ssn, range_constrain_value, "taildown"),
    euclid_de = exp(fill_optim_par_val_ssn[["euclid_de_log"]]),
    euclid_range = optim2orig_range_component(fill_optim_par_val_ssn, range_constrain_value, "euclid")
  )

  if ("euclid_extra_log" %in% names(fill_optim_par_val_ssn)) {
    euclid_values <- c(euclid_values, euclid_extra = exp(fill_optim_par_val_ssn[["euclid_extra_log"]]))
  } else if ("euclid_extra_logodds" %in% names(fill_optim_par_val_ssn)) {
    if (is.null(euclid_type)) {
      stop("euclid_type is required to inverse-transform euclid_extra.", call. = FALSE)
    }
    if (euclid_type == "matern") {
      extra <- 0.2 + (5 - 0.2) * expit(fill_optim_par_val_ssn[["euclid_extra_logodds"]])
    } else if (euclid_type == "pexponential") {
      extra <- 2 * expit(fill_optim_par_val_ssn[["euclid_extra_logodds"]])
    } else {
      stop("euclid_extra_logodds is only defined for matern or pexponential covariance.", call. = FALSE)
    }
    euclid_values <- c(euclid_values, euclid_extra = extra)
  }

  c(
    euclid_values[1:4],
    euclid_values[-(1:4)],
    euclid_rotate = pi * expit(fill_optim_par_val_ssn[["euclid_rotate_logodds"]]),
    euclid_scale = 1 * expit(fill_optim_par_val_ssn[["euclid_scale_logodds"]]),
    nugget = exp(fill_optim_par_val_ssn[["nugget_log"]])
  )
}

#' Back-transform a single range parameter from the optimization scale
#'
#' Which transform was used (log vs. logit-odds) is determined by which
#' optimization-scale name is present, mirroring how \code{euclid_extra}'s
#' own log/logit-odds choice is already resolved just above -- not by
#' re-deriving the original constrain decision.
#'
#' @param fill_optim_par_val_ssn A named vector of optimization-scale spatial
#'   covariance parameters.
#' @param range_constrain_value The shared range bound (only used on the
#'   logit-odds branch).
#' @param prefix \code{"tailup"}, \code{"taildown"}, or \code{"euclid"}.
#'
#' @return A single original-scale range value.
#'
#' @noRd
optim2orig_range_component <- function(fill_optim_par_val_ssn, range_constrain_value, prefix) {
  log_name <- paste0(prefix, "_range_log")
  if (log_name %in% names(fill_optim_par_val_ssn)) {
    exp(fill_optim_par_val_ssn[[log_name]])
  } else {
    logodds_name <- paste0(prefix, "_range_logodds")
    expit(fill_optim_par_val_ssn[[logodds_name]]) * range_constrain_value
  }
}

#' Back-transform random-effect variance parameters from the optimization scale
#'
#' @param par_randcov A named vector of optimization-scale (log) random-effect
#'   variances, or \code{NULL}.
#'
#' @return A named vector of original-scale random-effect variances, or
#'   \code{NULL} if \code{par_randcov} is \code{NULL}.
#'
#' @noRd
optim2orig_randcov_components <- function(par_randcov) {
  if (is.null(par_randcov)) {
    return(NULL)
  }
  val <- exp(par_randcov)
  names(val) <- gsub("_log", "", names(par_randcov))
  val
}
