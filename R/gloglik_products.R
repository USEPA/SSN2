#' Get products needed to evaluate log likelihood
#'
#' @param params_object Object with covariance parameters.
#' @param data_object Data object
#' @param estmethod Estimation Method
#'
#' @noRd
gloglik_products <- function(params_object, data_object, estmethod) {
  cov_matrix_list <- get_cov_matrix_list(params_object, data_object)
  # cholesky products (no local)
  # if (data_object$parallel) {
  #   cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
  #     cluster_list_element <- list(
  #       c = cov_matrix_list[[l]],
  #       x = data_object$X_list[[l]],
  #       y = data_object$y_list[[l]]
  #     )
  #   })
  #   cholprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_cholprods_parallel)
  #   names(cholprods_list) <- names(cov_matrix_list)
  # } else {
  #   cholprods_list <- mapply(
  #     c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
  #     function(c, x, y) get_cholprods(c, x, y),
  #     SIMPLIFY = FALSE
  #   )
  # }

  cholprods_list <- mapply(
    c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
    function(c, x, y) get_cholprods(c, x, y),
    SIMPLIFY = FALSE
  )

  # storing relevant products
  ## lower chol %*% X
  SqrtSigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_X))
  ## lower chol %*% y
  SqrtSigInv_y <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_y))
  # covariance of beta hat
  ## t(X) %*% sigma_inverse %*% X
  Xt_SigInv_X <- crossprod(SqrtSigInv_X, SqrtSigInv_X)
  ## t(X) %*% sigma_inverse %*% X)^(-1)
  Xt_SigInv_X_upchol <- chol(Xt_SigInv_X)
  cov_betahat <- chol2inv(Xt_SigInv_X_upchol)
  ## t(X) %*% sigma_inverse %*% y
  Xt_SigInv_y <- crossprod(SqrtSigInv_X, SqrtSigInv_y)
  ## t(X) %*% sigma_inverse %*% X)^(-1) %*% t(X) %*% sigma_inverse %*% y
  betahat <- cov_betahat %*% Xt_SigInv_y
  ## lower chol %*% (y - X %*% beta) r stands for residual
  SqrtSigInv_r <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
  ## residual %*% sigma_inverse %*% residual
  rt_SigInv_r <- crossprod(SqrtSigInv_r, SqrtSigInv_r)

  # using wolfinger notation
  l1 <- sum(unlist(lapply(cholprods_list, function(x) 2 * sum(log(diag(x$Sig_lowchol))))))
  l2 <- as.numeric(rt_SigInv_r)

  # returning relevant quantities
  if (estmethod == "reml") {
    l3 <- 2 * sum(log(diag(Xt_SigInv_X_upchol)))
    return(list(l1 = l1, l2 = l2, l3 = l3))
  }

  if (estmethod == "ml") {
    return(list(l1 = l1, l2 = l2))
  }
}

#' Get minus twice the log likelihood
#'
#' @param gloglik_products The relevant log likelihood products
#' @param data_object Data object
#' @param estmethod Estimation method
#'
#' @noRd
get_minustwologlik <- function(gloglik_products, data_object, estmethod) {
  if (estmethod == "reml") {
    minustwologlik <- as.numeric(gloglik_products$l1 + gloglik_products$l2 + gloglik_products$l3 + (data_object$n - data_object$p) * log(2 * pi))
  } else if (estmethod == "ml") {
    minustwologlik <- as.numeric(gloglik_products$l1 + gloglik_products$l2 + data_object$n * log(2 * pi))
  }
  minustwologlik
}
