#' Transform SSN spatial covariance starting values onto the optimization scale
#'
#' The inverse of \code{\link{optim2orig_ssn_components}()}: variance and
#' range parameters are log-transformed (unconstrained), and Euclidean
#' rotate/scale/shape parameters (each with a bounded original-scale range)
#' are logit-transformed, so the optimizer can search over an unconstrained
#' space.
#'
#' @param initial_object A joint covariance initial-value object with
#'   starting values for every field.
#' @param data_object A model data object with \code{range_constrain_value}
#'   and \code{tailup_range_constrain}/\code{taildown_range_constrain}/
#'   \code{euclid_range_constrain} (see
#'   \code{\link{get_range_constrain_setup}()}). Defaults to \code{NULL}
#'   (every range parameter unconstrained/log-scale), matching this
#'   function's behavior before \code{range_constrain} existed.
#'
#' @return A list with \code{value} (a named vector of optimization-scale
#'   parameters) and \code{is_known} (a matching named logical vector).
#'
#' @noRd
orig2optim_ssn_components <- function(initial_object, data_object = NULL) {
  ## tailup
  tailup_de <- initial_object$tailup_initial$initial[["de"]]
  tailup_de_log <- log(tailup_de)
  tailup_de_is_known <- initial_object$tailup_initial$is_known[["de"]]

  tailup_range <- orig2optim_range_component(
    initial_object$tailup_initial, data_object$tailup_range_constrain,
    data_object$range_constrain_value, "tailup"
  )

  ## taildown
  taildown_de <- initial_object$taildown_initial$initial[["de"]]
  taildown_de_log <- log(taildown_de)
  taildown_de_is_known <- initial_object$taildown_initial$is_known[["de"]]

  taildown_range <- orig2optim_range_component(
    initial_object$taildown_initial, data_object$taildown_range_constrain,
    data_object$range_constrain_value, "taildown"
  )

  ## euclid
  euclid_de <- initial_object$euclid_initial$initial[["de"]]
  euclid_de_log <- log(euclid_de)
  euclid_de_is_known <- initial_object$euclid_initial$is_known[["de"]]

  euclid_range <- orig2optim_range_component(
    initial_object$euclid_initial, data_object$euclid_range_constrain,
    data_object$range_constrain_value, "euclid"
  )

  euclid_type <- remove_covtype(class(initial_object$euclid_initial)[1])
  if (euclid_has_extra(euclid_type)) {
    euclid_extra <- initial_object$euclid_initial$initial[["extra"]]
    euclid_extra_is_known <- initial_object$euclid_initial$is_known[["extra"]]
    if (euclid_type == "matern") {
      euclid_extra_value <- logit((euclid_extra - 0.2) / (5 - 0.2))
      euclid_extra_name <- "euclid_extra_logodds"
    } else if (euclid_type == "cauchy") {
      euclid_extra_value <- log(euclid_extra)
      euclid_extra_name <- "euclid_extra_log"
    } else {
      euclid_extra_value <- logit(euclid_extra / 2)
      euclid_extra_name <- "euclid_extra_logodds"
    }
    euclid_extra_values <- stats::setNames(euclid_extra_value, euclid_extra_name)
    euclid_extra_known <- stats::setNames(euclid_extra_is_known, euclid_extra_name)
  } else {
    euclid_extra_values <- NULL
    euclid_extra_known <- NULL
  }

  euclid_rotate <- initial_object$euclid_initial$initial[["rotate"]]
  euclid_rotate_prop <- euclid_rotate / pi
  euclid_rotate_logodds <- logit(euclid_rotate_prop)
  euclid_rotate_is_known <- initial_object$euclid_initial$is_known[["rotate"]]

  euclid_scale <- initial_object$euclid_initial$initial[["scale"]]
  euclid_scale_logodds <- logit(euclid_scale)
  euclid_scale_is_known <- initial_object$euclid_initial$is_known[["scale"]]

  ## nugget
  nugget <- initial_object$nugget_initial$initial[["nugget"]]
  nugget_log <- log(nugget)
  nugget_is_known <- initial_object$nugget_initial$is_known[["nugget"]]

  list(
    value = c(
      tailup_de_log = tailup_de_log,
      tailup_range$value,
      taildown_de_log = taildown_de_log,
      taildown_range$value,
      euclid_de_log = euclid_de_log,
      euclid_range$value,
      euclid_extra_values,
      euclid_rotate_logodds = euclid_rotate_logodds,
      euclid_scale_logodds = euclid_scale_logodds,
      nugget_log = nugget_log
    ),
    is_known = c(
      tailup_de_is_known = tailup_de_is_known,
      tailup_range$is_known,
      taildown_de_is_known = taildown_de_is_known,
      taildown_range$is_known,
      euclid_de_is_known = euclid_de_is_known,
      euclid_range$is_known,
      euclid_extra_known,
      euclid_rotate_is_known = euclid_rotate_is_known,
      euclid_scale_is_known = euclid_scale_is_known,
      nugget_is_known = nugget_is_known
    )
  )
}

#' Transform a single range parameter onto the optimization scale
#'
#' Log scale (unconstrained) unless \code{range_constrain} requests (and
#' \code{\link{get_range_constrain_setup}()} approved) bounding this specific
#' range parameter, in which case a logit-odds scale bounded by
#' \code{range_constrain_value} is used instead -- the per-component analogue
#' of spmodel's single-range \code{orig2optim_range()}.
#'
#' @param component_initial A single covariance component's initial-value
#'   object (\code{tailup_initial}, \code{taildown_initial}, or
#'   \code{euclid_initial}).
#' @param constrain Whether this specific range parameter is constrained.
#' @param range_constrain_value The shared range bound (ignored when
#'   \code{constrain} is not \code{TRUE}).
#' @param prefix \code{"tailup"}, \code{"taildown"}, or \code{"euclid"} (used
#'   to name the returned value).
#'
#' @return A list with \code{value} and \code{is_known}, each a single named
#'   entry (\code{<prefix>_range_log} or \code{<prefix>_range_logodds}).
#'
#' @noRd
orig2optim_range_component <- function(component_initial, constrain, range_constrain_value, prefix) {
  range <- component_initial$initial[["range"]]
  is_known <- component_initial$is_known[["range"]]
  if (isTRUE(constrain)) {
    value <- logit(range / range_constrain_value)
    name <- paste0(prefix, "_range_logodds")
  } else {
    value <- log(range)
    name <- paste0(prefix, "_range_log")
  }
  list(
    value = stats::setNames(value, name),
    is_known = stats::setNames(is_known, name)
  )
}

#' Transform random-effect variance starting values onto the optimization scale
#'
#' The random-effect analogue of \code{\link{orig2optim_ssn_components}()}:
#' log-transforms each variance (clamped away from extreme values via
#' \code{\link{clamp_optim_scale}()} for numerical stability).
#'
#' @param initial_object A joint covariance initial-value object; only
#'   \code{randcov_initial} is used.
#'
#' @return A list with \code{value} (a named vector of optimization-scale
#'   log-variances), \code{is_known} (a matching named logical vector), and
#'   \code{n_est} (the number of not-known fields). All \code{NULL}/0 if
#'   there is no random effect.
#'
#' @noRd
orig2optim_randcov_components <- function(initial_object) {
  if (is.null(initial_object$randcov_initial)) {
    return(list(value = NULL, is_known = NULL, n_est = 0))
  }

  value <- log(initial_object$randcov_initial$initial)
  names(value) <- paste(names(initial_object$randcov_initial$initial), "log", sep = "_")
  is_known <- initial_object$randcov_initial$is_known
  names(is_known) <- paste(names(initial_object$randcov_initial$is_known), "log", sep = "_")
  n_est <- sum(!is_known)
  value <- clamp_optim_scale(value, is_known)

  list(value = value, is_known = is_known, n_est = n_est)
}

#' Clamp not-known optimization-scale values to a numerically safe range
#'
#' @param value A named vector of optimization-scale values.
#' @param is_known A matching named logical vector; known values are left
#'   untouched.
#'
#' @return \code{value} with not-known entries clamped to \code{[-50, 50]}.
#'
#' @noRd
clamp_optim_scale <- function(value, is_known) {
  value <- ifelse(value > 50 & !is_known, 50, value)
  value <- ifelse(value < -50 & !is_known, -50, value)
  value
}
