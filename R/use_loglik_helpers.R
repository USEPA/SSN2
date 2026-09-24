#' Build a placeholder optimizer-output object when every parameter is already known
#'
#' Used when no covariance parameters need estimating (so \code{optim()} is
#' never called), keeping the fitted model's \code{optim} element shaped
#' consistently with a real optimizer result.
#'
#' @param value The (already known) objective function value.
#'
#' @return A list with \code{method}, \code{control}, \code{value},
#'   \code{counts}, \code{convergence}, \code{message} (all \code{NA} except
#'   \code{value}), and \code{hessian} (\code{NA}).
#'
#' @noRd
known_optim_output_stub <- function(value) {
  list(
    method = NA,
    control = NA, value = value,
    counts = NA, convergence = NA,
    message = NA,
    hessian = NA
  )
}

#' Package \code{optim()}'s output (and its call settings) for storage on a fitted model
#'
#' @param optim_output The raw result of \code{optim()}.
#' @param optim_dotlist The \code{method}/\code{control}/\code{hessian}
#'   arguments passed to \code{optim()}.
#'
#' @return A list with \code{method}, \code{control}, \code{value},
#'   \code{counts}, \code{convergence}, \code{message}, and \code{hessian}
#'   (\code{FALSE} if not requested).
#'
#' @noRd
trim_optim_output <- function(optim_output, optim_dotlist) {
  list(
    method = optim_dotlist$method,
    control = optim_dotlist$control, value = optim_output$value,
    counts = optim_output$counts, convergence = optim_output$convergence,
    message = optim_output$message,
    hessian = if (optim_dotlist$hessian) optim_output$hessian else FALSE
  )
}

#' Resolve an anisotropy rotation angle's reflective ambiguity
#'
#' Anisotropic rotation has a mirror-image ambiguity (rotate and \eqn{\pi -
#' \mathrm{rotate}} can fit equally well); compares the fitted rotation
#' against its reflection by (negative twice) log-likelihood and keeps
#' whichever fits better.
#'
#' @param params_object A joint covariance parameter object with the fitted
#'   rotation.
#' @param data_object A model data object.
#' @param estmethod The estimation method (\code{"reml"} or \code{"ml"}).
#' @param is_glm Whether to use the Laplace (GLM) log-likelihood instead of
#'   the Gaussian one.
#'
#' @return \code{params_object}, with \code{euclid$rotate} replaced by its
#'   reflection if that fits better.
#'
#' @noRd
resolve_anis_rotation <- function(params_object, data_object, estmethod, is_glm) {
  if (is_glm) {
    products_fn <- laploglik_products
    minustwoll_fn <- get_minustwolaploglik
  } else {
    products_fn <- gloglik_products
    minustwoll_fn <- get_minustwologlik
  }

  prods <- products_fn(params_object, data_object, estmethod)
  minustwoll <- minustwoll_fn(prods, data_object, estmethod)

  params_object_q2 <- params_object
  params_object_q2$euclid[["rotate"]] <- pi - params_object_q2$euclid[["rotate"]]
  prods_q2 <- products_fn(params_object_q2, data_object, estmethod)
  minustwoll_q2 <- minustwoll_fn(prods_q2, data_object, estmethod)

  rotate_min <- which.min(c(minustwoll, minustwoll_q2))
  if (rotate_min == 2) {
    params_object <- params_object_q2
  }
  params_object
}

#' Extract the base (non-random-effect) known-parameter flags from an initial object
#'
#' @param initial_object A joint covariance initial-value object.
#'
#' @return A list with \code{tailup}, \code{taildown}, \code{euclid}, and
#'   \code{nugget}, each a named logical vector of which fields are known.
#'
#' @noRd
loglik_is_known_base <- function(initial_object) {
  list(
    tailup = initial_object$tailup_initial$is_known,
    taildown = initial_object$taildown_initial$is_known,
    euclid = initial_object$euclid_initial$is_known,
    nugget = initial_object$nugget_initial$is_known
  )
}
