#' @rdname loocv.SSN2
#' @param type The scale (\code{response} or \code{link}) of predictions obtained
#'   when \code{cv_predict = TRUE} and using \code{ssn_glm()} objects.
#' @param delta A logical indicating whether to return delta method standard errors
#' on the response scale when \code{se.fit = TRUE} and \code{type = "response"}. The default is \code{FALSE}.
#' @method loocv ssn_glm
#' @export
loocv.ssn_glm <- function(object, cv_predict = FALSE, type = c("link", "response"),
                          se.fit = FALSE, delta = FALSE, local, ...) {
  if (missing(local)) local <- NULL

  # match type argument so the two display
  type <- match.arg(type)
  if (!is.logical(delta) || length(delta) != 1 || is.na(delta)) {
    stop("delta must be TRUE or FALSE.", call. = FALSE)
  }
  local <- resolve_cv_auto_local(local, object$n, "loocv")
  local <- resolve_loocv_local(local)

  loocv_val <- get_loocv.ssn_glm(object, cv_predict = TRUE, se.fit = TRUE, local = local)
  response_val <- object$y
  error_val <- response_val - loocv_val$cv_predict
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
      if (type == "link") {
        loocv_out$cv_predict <- loocv_val$cv_predict_link
      } else if (type == "response") {
        loocv_out$cv_predict <- loocv_val$cv_predict
      } else {
        stop("Invalid type argument.", call. = FALSE)
      }
    }

    if (se.fit) {
      loocv_out$se.fit <- if (type == "response" && delta) {
        get_delta_se(loocv_val$cv_predict_link, loocv_val$se.fit, object$family, object$size)
      } else {
        loocv_val$se.fit
      }
    }
    return(loocv_out)
  }
}

get_loocv.ssn_glm <- function(object, cv_predict = FALSE, se.fit = FALSE, local = FALSE, ...) {
  local_list <- resolve_loocv_local(local)

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
    model_frame <- model.frame(object)
    w_linpred <- fitted(object, type = "link")
    model_offset <- model.offset(model_frame)
    # Kriging uses the latent predictor without the known offset. The family
    # derivatives remain evaluated at the offset-inclusive linear predictor.
    w <- if (is.null(model_offset)) w_linpred else w_linpred - as.vector(model_offset)
    size <- object$size

    # some products
    SigInv_w <- SigInv %*% w
    wX <- cbind(w, X)
    SigInv_wX <- cbind(SigInv_w, SigInv_X)

    # find H stuff
    wts_beta <- tcrossprod(cov_betahat, SigInv_X)
    Ptheta <- SigInv - SigInv_X %*% wts_beta
    d <- get_d(object$family, w_linpred, y, size, dispersion)
    # and then the gradient vector
    # g <-  d - Ptheta %*% w
    # Next, compute H
    D <- get_D(object$family, w_linpred, y, size, dispersion)
    H <- D - Ptheta
    mHinv <- solve(-H) # chol2inv(chol(Matrix::forceSymmetric(-H))) # solve(-H)

    cv_predict_val_list <- run_pred_dispatch(
      get_loocv_glm, seq_len(object$n), local_list,
      Sig = cov_matrix_val,
      SigInv = SigInv, Xmat = X, w = as.matrix(w, ncol = 1), wX = wX,
      SigInv_wX = SigInv_wX, mHinv = mHinv, se.fit = se.fit
    )

    cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    if (!is.null(model_offset)) {
      cv_predict_val <- cv_predict_val + as.vector(model_offset)
    }
    if (se.fit) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  } else {
    cov_matrix_val <- covmatrix(object)
    X <- model.matrix(object)
    model_frame <- model.frame(object)
    model_offset <- model.offset(model_frame)
    dispersion <- as.vector(coef(object, type = "dispersion"))
    cv_predict_val_list <- run_pred_dispatch(
      get_loocv_local_glm, seq_len(object$n), local_list,
      Sig = cov_matrix_val, Xmat = X, y = y, offset = model_offset,
      size = object$size, family = object$family, dispersion = dispersion,
      local_size = local_list$size, se.fit = se.fit,
      w = fitted(object, type = "link"), betahat = coef(object),
      cov_betahat = vcov(object, var_correct = FALSE)
    )
    if (se.fit) {
      cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    } else {
      cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    }
    if (!is.null(model_offset)) {
      cv_predict_val <- cv_predict_val + as.vector(model_offset)
    }
  }

  cv_predict_val_invlink <- invlink(cv_predict_val, object$family, object$size)

  if (cv_predict) {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), cv_predict = as.vector(cv_predict_val_invlink), cv_predict_link = as.vector(cv_predict_val), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- list(mspe = mean((cv_predict_val_invlink - y)^2), cv_predict = as.vector(cv_predict_val_invlink), cv_predict_link = as.vector(cv_predict_val))
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

get_loocv_local_glm <- function(obs, Sig, Xmat, y, offset, size, family,
                                dispersion, local_size, se.fit, w, betahat, cov_betahat) {
  retain <- get_loocv_local_index(obs, Sig, local_size)
  new_Sig <- Sig[retain, retain, drop = FALSE]
  new_Sig_upchol <- chol(forceSymmetric(new_Sig))
  new_Sig_lowchol <- t(new_Sig_upchol)
  new_SigInv <- chol2inv(new_Sig_upchol)
  new_X <- Xmat[retain, , drop = FALSE]
  new_y <- y[retain]
  new_size <- if (is.null(size)) NULL else size[retain]
  new_w <- w[retain]
  new_w_free <- if (is.null(offset)) new_w else new_w - offset[retain]
  obs_c <- Sig[obs, retain, drop = FALSE]
  obs_c_new_SigInv <- obs_c %*% new_SigInv
  obs_c_new_SigInv_X <- obs_c_new_SigInv %*% new_X
  new_pred <- Xmat[obs, , drop = FALSE] %*% betahat +
    obs_c_new_SigInv %*% (new_w_free - new_X %*% betahat)
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c_new_SigInv_X
    var_fit <- Sig[obs, obs] - tcrossprod(obs_c_new_SigInv, obs_c) +
      Q %*% tcrossprod(cov_betahat, Q)
    var_adj <- as.numeric(var_fit) + get_wts_varw(
      family, new_X, new_y, new_w, new_size, dispersion,
      new_Sig_lowchol, Xmat[obs, , drop = FALSE], obs_c
    )
    se_fit <- sqrt(var_adj)
  } else {
    se_fit <- NULL
  }
  list(pred = as.numeric(new_pred), se.fit = se_fit)
}

#' Get the exact (non-local) loocv prediction and standard error for GLM-type models
#'
#' @param obs An observation to leave out
#' @param Sig The full covariance matrix
#' @param SigInv The full inverse covariance matrix
#' @param Xmat Model matrix
#' @param w The latent (link-scale) predictor vector
#' @param wX \code{cbind(w, Xmat)}
#' @param SigInv_wX \code{SigInv \%*\% wX}
#' @param mHinv The inverse of the negative Hessian of the Laplace log-likelihood
#' @param se.fit Whether to compute the standard error
#'
#' @return A list with elements \code{pred} (the link-scale loocv prediction)
#'   and \code{se.fit} (its standard error, or \code{NULL} if \code{se.fit} is \code{FALSE}),
#'   computed via a partitioned-inverse (Sherman-Morrison-type) update rather
#'   than refitting the model with the observation removed
#'
#' @noRd
get_loocv_glm <- function(obs, Sig, SigInv, Xmat, w, wX, SigInv_wX, mHinv, se.fit) {
  # like get_loocv(), this avoids literally refitting with "obs" dropped by
  # updating SigInv via a partitioned-matrix (Sherman-Morrison-type) formula;
  # the GLM case works on the link-scale latent predictor w (from the Laplace
  # approximation) rather than the response y directly
  SigInv_mm <- SigInv[obs, obs] # a constant
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  neww <- w[-obs, , drop = FALSE]
  newX <- Xmat[-obs, , drop = FALSE]
  newwX <- wX[-obs, , drop = FALSE]

  # SigInv for the data with "obs" removed, via the partitioned inverse update
  new_SigInv <- SigInv[-obs, -obs] - tcrossprod(SigInv_om, SigInv_om) / SigInv_mm
  new_SigInv_newX <- new_SigInv %*% newX
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, new_SigInv_newX))))

  new_wts_beta <- tcrossprod(new_covbetahat, new_SigInv_newX)
  obs_c <- Sig[obs, -obs, drop = FALSE]
  obs_c_new_SigInv <- obs_c %*% new_SigInv
  obs_c_new_SigInv_newX <- obs_c %*% new_SigInv_newX
  # weights that map the remaining observations' latent w to the kriging
  # prediction at "obs" (universal kriging equation on the link scale)
  new_wts_pred <- Xmat[obs, , drop = FALSE] %*% new_wts_beta + obs_c %*% new_SigInv - obs_c_new_SigInv_newX %*% new_wts_beta

  new_pred <- new_wts_pred %*% neww

  # var
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c_new_SigInv_newX
    var_fit <- Sig[obs, obs] - tcrossprod(obs_c_new_SigInv, obs_c) + Q %*% tcrossprod(new_covbetahat, Q)
    # mHinv (inverse negative Hessian of the Laplace loglik) captures the
    # extra uncertainty in w from the Laplace approximation itself; update it
    # the same way as SigInv above, then fold that uncertainty into var_fit
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
