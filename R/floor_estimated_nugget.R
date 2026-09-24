#' Floor an estimated nugget variance away from (near-)zero
#'
#' Leaves a known (fixed) nugget untouched. Otherwise, floors the estimated
#' nugget at the larger of 0.0001 times the total spatially dependent
#' variance and \code{diagtol}, avoiding numerically unstable near-singular
#' covariance matrices from an estimated nugget collapsing to (near) zero.
#'
#' @param params_object A joint covariance parameter object.
#' @param nugget_is_known A named logical vector/list indicating whether
#'   each nugget field is known (fixed).
#' @param diagtol A numerical-stability tolerance added to the covariance
#'   diagonal.
#'
#' @return \code{params_object} with its nugget floored if estimated.
#'
#' @noRd
floor_estimated_nugget <- function(params_object, nugget_is_known, diagtol) {
  if (isTRUE(nugget_is_known[["nugget"]])) {
    return(params_object)
  }
  de_scale <- sum(
    params_object$tailup[["de"]], params_object$taildown[["de"]],
    params_object$euclid[["de"]]
  )
  params_object$nugget[["nugget"]] <- max(
    params_object$nugget[["nugget"]], 1e-4 * de_scale, diagtol
  )
  params_object
}
