#' @rdname loocv.SSN2
#' @method loocv ssn_glm
#' @export
loocv.ssn_glm <- function(object, cv_predict = FALSE, se.fit = FALSE, ...) {
  loocv_val <- get_loocv.ssn_glm(object, cv_predict = TRUE, se.fit = TRUE)
  response_val <- object$y
  error_val <- loocv_val$cv_predict - response_val
  se_val <- loocv_val$se.fit

  bias <- mean(error_val)
  MSPE <- loocv_val$mspe
  RMSPE <- sqrt(loocv_val$mspe)
  RAV <- sqrt(mean(se_val^2))

  # std.bias does not really make sense as se is on the link scale
  # std.RMSPE does not really make sense as se is on the link scale
  # but response is on response scale

  # coverage does not really make sense because we don't observe the
  # latent means (cover is not for the response but rather for the latent means)

  loocv_stats <- tibble(
    bias = bias,
    MSPE = MSPE,
    RMSPE = RMSPE,
    RAV = RAV
  )

  if (!cv_predict && !se.fit) {
    return(loocv_stats)
  } else {
    loocv_out <- list()
    loocv_out$stats <- loocv_stats

    if (cv_predict) {
      loocv_out$cv_predict <- loocv_val$cv_predict
    }

    if (se.fit) {
      loocv_out$se.fit <- loocv_val$se.fit
    }
    return(loocv_out)
  }
}

#' A helper to get the leave one out predictions for glms
#'
#' @param object Model object
#' @param cv_predict Whether cross validation predictions should be returned
#' @param se.fit Whether the standard error should be returned
#' @param ... Additional arguments
#'
#' @noRd
get_loocv.ssn_glm <- function(object, cv_predict = FALSE, se.fit = FALSE, ...) {
  # turn local off for now
  local <- FALSE
  local_list <- get_local_list_prediction(local)

  # store response
  y <- object$y

  if (local_list$method == "all") {
    cov_matrix_val <- covmatrix(object)
    X <- model.matrix(object)
    cholprods <- get_cholprods_glm(cov_matrix_val, X, y)
    # actually need inverse because of HW blocking
    SigInv <- chol2inv(cholprods$Sig_lowchol)
    SigInv_X <- backsolve(t(cholprods$Sig_lowchol), cholprods$SqrtSigInv_X)

    # find products
    Xt_SigInv_X <- crossprod(X, SigInv_X)
    Xt_SigInv_X_upchol <- base::chol(Xt_SigInv_X) # or Matrix::forceSymmetric()
    cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

    # glm stuff
    dispersion <- as.vector(coef(object, type = "dispersion")) # take class away
    w <- fitted(object, type = "link")
    size <- object$size

    # some products
    SigInv_w <- SigInv %*% w
    wX <- cbind(w, X)
    SigInv_wX <- cbind(SigInv_w, SigInv_X)

    # find H stuff
    wts_beta <- tcrossprod(cov_betahat, SigInv_X)
    Ptheta <- SigInv - SigInv_X %*% wts_beta
    d <- get_d(object$family, w, y, size, dispersion)
    # and then the gradient vector
    # g <-  d - Ptheta %*% w
    # Next, compute H
    D <- get_D(object$family, w, y, size, dispersion)
    H <- D - Ptheta
    mHinv <- solve(-H) # chol2inv(chol(Matrix::forceSymmetric(-H))) # solve(-H)

    # parallel stuff
    # if (local_list$parallel) {
    #   cl <- parallel::makeCluster(local_list$ncores)
    #   cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv_glm,
    #                                              Sig = cov_matrix_val,
    #                                              SigInv = SigInv, Xmat = X, w = w, wX = wX,
    #                                              SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
    #   )
    #   cl <- parallel::stopCluster(cl)
    # } else {
    #   cv_predict_val_list <- lapply(seq_len(object$n), get_loocv_glm,
    #                                 Sig = cov_matrix_val,
    #                                 SigInv = SigInv, Xmat = X, w = as.matrix(w, ncol = 1), wX = wX,
    #                                 SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
    #   )
    # }

    cv_predict_val_list <- lapply(seq_len(object$n), get_loocv_glm,
      Sig = cov_matrix_val,
      SigInv = SigInv, Xmat = X, w = as.matrix(w, ncol = 1), wX = wX,
      SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
    )

    cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    if (se.fit) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  } else {
    # get w for later
    w <- fitted(object, type = "link")

    if (local_list$parallel) {
      # turn of parallel as it is used different in predict
      local_list$parallel <- FALSE
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), loocv_local_glm, object, se.fit, local_list)
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), loocv_local_glm, object, se.fit, local_list)
    }
    if (se.fit) {
      cv_predict_val <- vapply(cv_predict_val_list, function(x) x$fit, numeric(1))
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    } else {
      cv_predict_val <- unlist(cv_predict_val_list)
    }
  }

  cv_predict_val_invlink <- invlink(cv_predict_val, object$family, object$size)

  if (cv_predict) {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), cv_predict = as.vector(cv_predict_val_invlink), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), cv_predict = as.vector(cv_predict_val_invlink))
    }
  } else {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- mean((cv_predict_val_invlink - y)^2)
    }
  }
  cv_output
}

#' Helper to get each row leave one out prediction
#'
#' @param row The row of interest
#' @param object Model object
#' @param se.fit Whether to return the standard error
#' @param local_list Local list (not yet functional)
#'
#' @noRd
loocv_local_glm <- function(row, object, se.fit, local_list) {
  object$fitted$link <- object$fitted$link[-row] # w needs to be subset
  newdata <- object$obdata[row, , drop = FALSE]
  object$obdata <- object$obdata[-row, , drop = FALSE]
  predict(object, newdata = newdata, se.fit = se.fit, local = local_list)
}

#' Get leave one out cross validation value for glm
#'
#' @param obs Index of observed data
#' @param Sig Covariance matrix
#' @param SigInv Inverse covariance matrix
#' @param Xmat Model matrix
#' @param w The predicted latent effects
#' @param wX Matrix containing w and response
#' @param SigInv_yX Product of inverse covariance matrix with w and response
#' @param mHinv The minus of the inverse Hessian of w
#' @param se.fit Whether to return standard errors of leave one out predictions
#'
#' @noRd
get_loocv_glm <- function(obs, Sig, SigInv, Xmat, w, wX, SigInv_wX, mHinv, se.fit) {
  SigInv_mm <- SigInv[obs, obs] # a constant
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  neww <- w[-obs, , drop = FALSE]
  newX <- Xmat[-obs, , drop = FALSE]
  newwX <- wX[-obs, , drop = FALSE]

  new_SigInv <- SigInv[-obs, -obs] - tcrossprod(SigInv_om, SigInv_om) / SigInv_mm
  new_SigInv_newX <- new_SigInv %*% newX
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, new_SigInv_newX))))

  new_wts_beta <- tcrossprod(new_covbetahat, new_SigInv_newX)
  obs_c <- Sig[obs, -obs, drop = FALSE]
  obs_c_new_SigInv <- obs_c %*% new_SigInv
  obs_c_new_SigInv_newX <- obs_c %*% new_SigInv_newX
  new_wts_pred <- Xmat[obs, , drop = FALSE] %*% new_wts_beta + obs_c %*% new_SigInv - obs_c_new_SigInv_newX %*% new_wts_beta


  new_pred <- new_wts_pred %*% neww

  # var
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c_new_SigInv_newX
    var_fit <- Sig[obs, obs] - tcrossprod(obs_c_new_SigInv, obs_c) + Q %*% tcrossprod(new_covbetahat, Q)
    mHinv_mm <- mHinv[obs, obs]
    mHinv_om <- mHinv[-obs, obs, drop = FALSE]
    newmHinv <- mHinv[-obs, -obs] - tcrossprod(mHinv_om, mHinv_om) / mHinv_mm
    var_adj <- as.numeric(var_fit + new_wts_pred %*% tcrossprod(newmHinv, new_wts_pred))
    se_fit <- sqrt(var_adj)
  } else {
    se_fit <- NULL
  }

  # return
  list(pred = as.numeric(new_pred), se.fit = as.numeric(se_fit))
}
