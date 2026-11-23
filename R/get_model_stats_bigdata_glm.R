#' Get relevant model fit statistics and diagnostics for glms
#'
#' @param cov_est_object Covariance parameter estimation object.
#' @param data_object Data object.
#' @param estmethod Estimation method.
#'
#' @noRd
get_model_stats_bigdata_glm <- function(cov_est_object, data_object, estmethod) {
  cov_matrix_list <- get_cov_matrix_list(cov_est_object$params_object, data_object)

  if (data_object$parallel) {
    cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
      cluster_list_element <- list(
        c = cov_matrix_list[[l]],
        x = data_object$X_list[[l]],
        y = data_object$y_list[[l]],
        o = data_object$ones_list[[l]]
      )
    })
    eigenprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_eigenprods_glm_parallel)
    names(eigenprods_list) <- names(cov_matrix_list)
  } else {
    eigenprods_list <- mapply(
      c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, o = data_object$ones_list,
      function(c, x, y, o) get_eigenprods_glm(c, x, y, o),
      SIMPLIFY = FALSE
    )
  }

  get_model_stats_glm_core(cov_est_object, data_object, estmethod, data_object$order_bigdata, eigenprods_list)
}
