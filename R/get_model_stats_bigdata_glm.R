#' Get relevant model fit statistics and diagnostics for glms
#'
#' @param cov_est_object Covariance parameter estimation object.
#' @param data_object Data object.
#' @param estmethod Estimation method.
#'
#' @noRd
get_model_stats_bigdata_glm <- function(cov_est_object, data_object, estmethod) {
  cov_matrix_list <- get_cov_matrix_list(cov_est_object$params_object, data_object)

  # find model components
  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)

  # eigen products (put back when local impelmented)
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

  # eigenprods_list <- mapply(
  #   c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, o = data_object$ones_list,
  #   function(c, x, y, o) get_eigenprods_glm(c, x, y, o),
  #   SIMPLIFY = FALSE
  # )

  SigInv_list <- lapply(eigenprods_list, function(x) x$SigInv)
  SigInv <- Matrix::bdiag(SigInv_list)
  SigInv_X <- do.call("rbind", lapply(eigenprods_list, function(x) x$SigInv_X))

  # get inverse cov beta hat list and add together
  invcov_betahat_list <- lapply(eigenprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_X))
  invcov_betahat_sum <- Reduce("+", invcov_betahat_list)
  # find unadjusted cov beta hat matrix
  cov_betahat_noadjust <- chol2inv(chol(forceSymmetric(invcov_betahat_sum)))
  # put it in a list the number of times there are unique local indices
  cov_betahat_noadjust_list <- rep(list(cov_betahat_noadjust), times = length(invcov_betahat_list))

  # dispersion
  dispersion <- as.vector(cov_est_object$params_object$dispersion)

  # newton rhapson
  w_and_H <- get_w_and_H(data_object, dispersion,
                         SigInv_list, SigInv_X, cov_betahat_noadjust,
                         invcov_betahat_sum, estmethod,
                         ret_mHInv = TRUE
  )

  w <- as.vector(w_and_H$w)
  # H <- w_and_H$H

  # put w in eigenprods
  w_list <- split(w, sort(data_object$local_index))

  Xt_SigInv_w_list <- mapply(
    x = eigenprods_list, w = w_list,
    function(x, w) crossprod(x$SigInv_X, w),
    SIMPLIFY = FALSE
  )

  betahat_list <- mapply(
    l = cov_betahat_noadjust_list, r = Xt_SigInv_w_list,
    function(l, r) l %*% r,
    SIMPLIFY = FALSE
  )

  betahat <- as.numeric(cov_betahat_noadjust %*%
                          Reduce("+", Xt_SigInv_w_list))
  names(betahat) <- colnames(data_object$X_list[[1]])

  # adjust the global covariance of beta hat (revisit when local implemented)
  cov_betahat <- cov_betahat_adjust(
    invcov_betahat_list,
    betahat_list, betahat,
    eigenprods_list, data_object,
    cov_est_object$params_object,
    cov_betahat_noadjust, data_object$var_adjust
  )

  cov_betahat <- as.matrix(cov_betahat)
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  betawtsvarw <- wts_beta %*% w_and_H$mHInv %*% t(wts_beta)

  cov_betahat_uncorrected <- cov_betahat # save uncorrected cov beta hat
  rownames(cov_betahat_uncorrected) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat_uncorrected) <- colnames(data_object$X_list[[1]])

  cov_betahat <- as.matrix(cov_betahat + betawtsvarw)
  rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat) <- colnames(data_object$X_list[[1]])

  # return fixed and random coefficients
  coefficients <- get_coefficients_glm(betahat, cov_est_object$params_object)

  # return fitted
  fitted <- get_fitted_glm(w_list, betahat, cov_est_object$params_object, data_object, eigenprods_list)

  # return hat values
  hatvalues <- get_hatvalues_glm(w, X, data_object, dispersion)

  # return deviance i
  deviance_i <- get_deviance_glm(data_object$family, y, fitted$response, data_object$size, dispersion)
  deviance_i <- pmax(deviance_i, 0) # sometimes numerical instability can cause these to be slightly negative

  # storing relevant products
  SigInv_X_null <- do.call("rbind", lapply(eigenprods_list, function(x) x$SigInv_ones))
  ## lower chol %*% X
  SqrtSigInv_X_null <- do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_ones))
  # covariance of beta hat
  ## t(X) %*% sigma_inverse %*% X
  Xt_SigInv_X_null <- crossprod(SqrtSigInv_X_null, SqrtSigInv_X_null)
  ## t(X) %*% sigma_inverse %*% X)^(-1)
  Xt_SigInv_X_upchol_null <- chol(Xt_SigInv_X_null)
  cov_betahat_null <- chol2inv(Xt_SigInv_X_upchol_null)

  # newton rhapson
  w_and_H_null <- get_w_and_H(
    data_object, dispersion,
    SigInv_list, SigInv_X_null, cov_betahat_null, Xt_SigInv_X_null, estmethod
  )

  w_null <- as.vector(w_and_H_null$w)

  fitted_null <- get_fitted_null(w_null, data_object)

  # return deviance i
  deviance_i_null <- get_deviance_glm(data_object$family, y, fitted_null, data_object$size, dispersion)
  deviance_i_null <- pmax(deviance_i_null, 0) # sometimes numerical instability can cause these to be slightly non-negative

  deviance <- sum(deviance_i)
  deviance_null <- sum(deviance_i_null)
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)

  # should always be non-negative
  pseudoR2 <- pmax(0, pseudoR2)
  # set null model R2 equal to zero (no covariates)
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }

  # return residuals
  residuals <- get_residuals_glm(w, y, data_object, deviance_i, hatvalues, dispersion)

  # return cooks distance
  cooks_distance <- get_cooks_distance_glm(residuals, hatvalues, data_object$p)

  # return variance covariance matrices
  vcov <- get_vcov_glm(cov_betahat, cov_betahat_uncorrected) # note first argument is the adjusted one

  # reorder relevant quantities to match data order
  model_stats_glm_names <- data_object$pid[data_object$observed_index]
  ## fitted values
  fitted$response <- fitted$response[order(data_object$order_bigdata)]
  names(fitted$response) <- model_stats_glm_names
  fitted$link <- fitted$link[order(data_object$order_bigdata)]
  names(fitted$link) <- model_stats_glm_names
  fitted$tailup <- fitted$tailup[order(data_object$order_bigdata)]
  names(fitted$tailup) <- model_stats_glm_names
  fitted$taildown <- fitted$taildown[order(data_object$order_bigdata)]
  names(fitted$taildown) <- model_stats_glm_names
  fitted$euclid <- fitted$euclid[order(data_object$order_bigdata)]
  names(fitted$euclid) <- model_stats_glm_names
  fitted$nugget <- fitted$nugget[order(data_object$order_bigdata)]
  names(fitted$nugget) <- model_stats_glm_names
  ## hat values
  hatvalues <- hatvalues[order(data_object$order_bigdata)]
  names(hatvalues) <- model_stats_glm_names
  ## residuals
  residuals$response <- residuals$response[order(data_object$order_bigdata)]
  names(residuals$response) <- model_stats_glm_names
  residuals$deviance <- residuals$deviance[order(data_object$order_bigdata)]
  names(residuals$deviance) <- model_stats_glm_names
  residuals$pearson <- residuals$pearson[order(data_object$order_bigdata)]
  names(residuals$pearson) <- model_stats_glm_names
  residuals$standardized <- residuals$standardized[order(data_object$order_bigdata)]
  names(residuals$standardized) <- model_stats_glm_names
  ## cook's distance
  cooks_distance <- cooks_distance[order(data_object$order_bigdata)]
  names(cooks_distance) <- model_stats_glm_names
  y <- y[order(data_object$order_bigdata)]
  if (is.null(data_object$size)) {
    size <- NULL
  } else {
    size <- data_object$size[order(data_object$order_bigdata)]
    names(size) <- model_stats_glm_names
  }

  # return npar (number of estimated covariance parameters)
  npar <- sum(unlist(lapply(cov_est_object$is_known, function(x) length(x) - sum(x))))

  # return list
  list(
    coefficients = coefficients,
    fitted = fitted,
    hatvalues = hatvalues,
    residuals = residuals,
    cooks_distance = cooks_distance,
    vcov = vcov,
    deviance = deviance,
    pseudoR2 = pseudoR2,
    npar = npar,
    w = w,
    y = y, # problems with model.response later if not here
    size = size
  )
}
