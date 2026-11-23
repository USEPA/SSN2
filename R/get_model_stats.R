#' Assemble fitted-model statistics for a general (non-IID) Gaussian fit
#'
#' @param cov_est_object A fitted covariance estimation object
#' @param data_object The data object
#' @param estmethod The estimation method (unused; kept for a consistent
#'   dispatch signature with \code{get_model_stats_iid()})
#'
#' @noRd
get_model_stats <- function(cov_est_object, data_object, estmethod) {
  get_model_stats_core(cov_est_object, data_object, data_object$order)
}

#' Assemble fitted-model statistics for an IID (no spatial dependence, no
#' random effects) fit
#'
#' @param cov_est_object A fitted covariance estimation object, as returned by
#'   \code{get_gloglik_iid()} (or \code{use_gloglik_known()} for a fully
#'   fixed nugget)
#' @param data_object The data object
#' @param estmethod The estimation method (unused; kept for a consistent
#'   dispatch signature with \code{get_model_stats()})
#'
#' @return The same statistics as \code{get_model_stats()}, computed directly
#'   from a QR decomposition of \code{X} rather than a full covariance
#'   matrix, since with no spatial dependence or random effects the
#'   covariance matrix is simply \code{nugget * I}
#'
#' @noRd
get_model_stats_iid <- function(cov_est_object, data_object, estmethod) {
  get_model_stats_iid_core(cov_est_object, data_object, data_object$order)
}
