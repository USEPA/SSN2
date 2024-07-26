#' Specify spatial indexing covariance matrix adjustment. Not currently
#'   relevant as spatial indexing has not yet been implemented.
#'
#' @param invcov_betahat_list Placeholder.
#' @param betahat_list Placeholder.
#' @param betahat Placeholder.
#' @param eigenprods_list Placeholder.
#' @param data_object Placeholder.
#' @param params_object Placeholder.
#' @param cov_betahat_noadjust Placeholder.
#' @param var_adjust Placeholder.
#'
#' @noRd
cov_betahat_adjust <- function(invcov_betahat_list, betahat_list,
                               betahat, eigenprods_list, data_object, params_object,
                               cov_betahat_noadjust, var_adjust) {
  # for now before local implemented
  var_adjust <- "none"
  cov_betahat_adjust_val <- cov_betahat_noadjust
  cov_betahat_adjust_val
}
