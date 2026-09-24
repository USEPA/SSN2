#' @param type The scale (\code{response} or \code{link}) of predictions obtained
#' when \code{cv_predict = TRUE} and using \code{ssn_lm()} or \code{ssn_glm} objects.
#' @param delta A logical indicating whether to return delta method standard errors
#' on the response scale when \code{se.fit = TRUE} and \code{type = "response"}. The default is \code{FALSE}.
#' @rdname kcv.SSN2
#' @method kcv ssn_glm
#' @export
kcv.ssn_glm <- function(object, k = 5, cv_predict = FALSE,
                        type = c("link", "response"), se.fit = FALSE,
                        delta = FALSE, local, folds_index, ...) {
  if (missing(folds_index)) folds_index <- NULL
  if (missing(local)) local <- NULL
  type <- match.arg(type)
  if (!is.logical(delta) || length(delta) != 1 || is.na(delta)) {
    stop("delta must be TRUE or FALSE.", call. = FALSE)
  }
  local <- resolve_cv_auto_local(local, object$n, "kcv")
  local <- resolve_kcv_local(local)

  X <- model.matrix(object)
  folds <- resolve_kcv_folds(object, k, folds_index, X)
  if (is_loocv_folds(folds$fold_list)) {
    return(loocv(
      object, cv_predict = cv_predict, type = type, se.fit = se.fit,
      delta = delta, local = local, ...
    ))
  }

  # RAV is part of SSN2's established GLM CV statistics, so link-scale
  # standard errors are needed even when the caller does not return them.
  if (local$method == "all") {
    cv <- get_kcv.ssn_glm(object, folds$fold_list, se.fit = TRUE, local = local)
  } else {
    cv <- get_kcv_local.ssn_glm(object, folds$fold_list, se.fit = TRUE, local = local)
  }
  stats <- get_kcv_glm_stats(cv$cv_predict, object$y, cv$se.fit)

  if (!cv_predict && !se.fit) {
    return(stats)
  }
  output <- list(stats = stats)
  if (cv_predict) {
    output$cv_predict <- if (type == "link") cv$cv_predict_link else cv$cv_predict
  }
  if (se.fit) {
    output$se.fit <- if (type == "response" && delta) {
      get_delta_se(cv$cv_predict_link, cv$se.fit, object$family, object$size)
    } else {
      cv$se.fit
    }
  }
  output
}

#' Compute GLM k-fold cross-validation error statistics
#'
#' @param cv_predict A vector of k-fold cross-validation predictions.
#' @param response The observed response vector.
#' @param se.fit A vector of link-scale k-fold cross-validation prediction
#'   standard errors.
#'
#' @return A tibble with \code{bias}, \code{MSPE}, \code{RMSPE}, and
#'   \code{RAV}.
#'
#' @noRd
get_kcv_glm_stats <- function(cv_predict, response, se.fit) {
  error <- response - cv_predict
  tibble(
    bias = mean(error),
    MSPE = mean(error^2),
    RMSPE = sqrt(mean(error^2)),
    RAV = sqrt(mean(se.fit^2))
  )
}

#' Compute exact (block-update) k-fold cross-validation predictions for a GLM
#'
#' Computes the full covariance matrix, its inverse, and the negative
#' Laplace-likelihood Hessian's inverse once, then applies
#' \code{\link{get_kcv_glm}()}'s partitioned-inverse block update per fold.
#'
#' @param object A fitted \code{ssn_glm} model object.
#' @param fold_list A named list of per-fold row-index vectors.
#' @param se.fit Whether to also return prediction standard errors.
#' @param local A resolved \code{local} list.
#'
#' @return A list with \code{cv_predict_link}, \code{cv_predict} (response
#'   scale), and \code{se.fit} (link scale, \code{NULL} if not requested).
#'
#' @noRd
get_kcv.ssn_glm <- function(object, fold_list, se.fit, local) {
  y <- object$y
  Sigma <- covmatrix(object)
  X <- model.matrix(object)
  cholprods <- get_cholprods_glm(Sigma, X, y)
  Sigma_inv <- chol2inv(cholprods$Sig_lowchol)
  Sigma_inv_X <- backsolve(t(cholprods$Sig_lowchol), cholprods$SqrtSigInv_X)
  cov_betahat <- chol2inv(chol(forceSymmetric(crossprod(X, Sigma_inv_X))))

  dispersion <- as.vector(coef(object, type = "dispersion"))
  frame <- model.frame(object)
  offset <- model.offset(frame)
  w_link <- fitted(object, type = "link")
  w <- if (is.null(offset)) w_link else w_link - as.vector(offset)
  wX <- cbind(w, X)
  Sigma_inv_wX <- Sigma_inv %*% wX
  weights_beta <- tcrossprod(cov_betahat, Sigma_inv_X)
  Ptheta <- Sigma_inv - Sigma_inv_X %*% weights_beta
  H <- get_D(object$family, w_link, y, object$size, dispersion) - Ptheta
  mHinv <- solve(-H)

  values <- run_pred_dispatch(
    get_kcv_glm, fold_list, local,
    Sig = Sigma, SigInv = Sigma_inv, Xmat = X, w = matrix(w, ncol = 1),
    wX = wX, SigInv_wX = Sigma_inv_wX, mHinv = mHinv, se.fit = se.fit
  )
  pred_link <- numeric(object$n)
  if (se.fit) se <- numeric(object$n)
  for (i in seq_along(fold_list)) {
    held <- fold_list[[i]]
    pred_link[held] <- values[[i]]$pred
    if (se.fit) se[held] <- values[[i]]$se.fit
  }
  if (!is.null(offset)) pred_link <- pred_link + as.vector(offset)
  list(
    cv_predict_link = pred_link,
    cv_predict = invlink(pred_link, object$family, object$size),
    se.fit = if (se.fit) se else NULL
  )
}

#' Get the exact (non-local) kcv prediction and standard error for GLM-type models
#'
#' Block generalization of \code{\link{get_loocv_glm}()}: holds out a whole
#' fold (a vector of row indices) at once instead of a single observation, via
#' the same partitioned-matrix update with \code{solve()} against the
#' \code{m x m} blocks \code{SigInv[fold, fold]}/\code{mHinv[fold, fold]}
#' (\code{m} = fold size) in place of division by a scalar. At \code{m = 1}
#' this reduces to exactly \code{get_loocv_glm()}'s formula.
#'
#' @param fold A vector of row indices to leave out together
#' @param Sig The full covariance matrix
#' @param SigInv The full inverse covariance matrix
#' @param Xmat Model matrix
#' @param w The latent (link-scale) predictor vector
#' @param wX \code{cbind(w, Xmat)}
#' @param SigInv_wX \code{SigInv \%*\% wX}
#' @param mHinv The inverse of the negative Hessian of the Laplace log-likelihood
#' @param se.fit Whether to compute the standard error
#'
#' @return A list with elements \code{pred} (the link-scale kcv predictions
#'   for the fold) and \code{se.fit} (their standard errors, or \code{NULL} if
#'   \code{se.fit} is \code{FALSE}), computed via a partitioned-inverse
#'   (Sherman-Morrison-type) update rather than refitting the model with the
#'   fold removed
#'
#' @noRd
get_kcv_glm <- function(fold, Sig, SigInv, Xmat, w, wX, SigInv_wX, mHinv, se.fit) {
  train <- -fold
  SigInv_mm <- SigInv[fold, fold, drop = FALSE] # an m x m block (a scalar when m = 1)
  SigInv_om <- SigInv[train, fold, drop = FALSE]
  new_w <- w[train, , drop = FALSE]
  new_X <- Xmat[train, , drop = FALSE]
  # SigInv for the data with "fold" removed, via the partitioned inverse update
  new_SigInv <- SigInv[train, train, drop = FALSE] -
    SigInv_om %*% solve(SigInv_mm, t(SigInv_om))
  new_SigInv_X <- new_SigInv %*% new_X
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(new_X, new_SigInv_X))))
  new_weights_beta <- tcrossprod(new_covbetahat, new_SigInv_X)
  held_c <- Sig[fold, train, drop = FALSE]
  held_c_SigInv <- held_c %*% new_SigInv
  held_c_SigInv_X <- held_c %*% new_SigInv_X
  # weights that map the remaining observations' latent w to the kriging
  # prediction at "fold" (universal kriging equation on the link scale)
  weights_pred <- Xmat[fold, , drop = FALSE] %*% new_weights_beta +
    held_c_SigInv - held_c_SigInv_X %*% new_weights_beta
  pred <- weights_pred %*% new_w

  if (se.fit) {
    Q <- Xmat[fold, , drop = FALSE] - held_c_SigInv_X
    variance <- Sig[fold, fold, drop = FALSE] - tcrossprod(held_c_SigInv, held_c) +
      Q %*% tcrossprod(new_covbetahat, Q)
    # mHinv (inverse negative Hessian of the Laplace loglik) captures the
    # extra uncertainty in w from the Laplace approximation itself; update it
    # the same way as SigInv above, then fold that uncertainty into variance
    mHinv_mm <- mHinv[fold, fold, drop = FALSE]
    mHinv_om <- mHinv[train, fold, drop = FALSE]
    new_mHinv <- mHinv[train, train, drop = FALSE] -
      mHinv_om %*% solve(mHinv_mm, t(mHinv_om))
    variance <- variance + weights_pred %*% tcrossprod(new_mHinv, weights_pred)
    # only the diagonal (marginal SE per point) is returned, matching
    # get_loocv_glm()'s per-observation output shape
    se <- sqrt(diag(variance))
  } else {
    se <- NULL
  }
  list(pred = as.numeric(pred), se.fit = se)
}

#' Refit an \code{ssn_glm} model with one fold held out, for local k-fold cross-validation
#'
#' @param object A fitted \code{ssn_glm} model object.
#' @param fold A vector of row indices (into the fitted rows) to hold out.
#' @param local A resolved \code{local} list.
#'
#' @return The refit \code{ssn_glm} model object, with the fold's rows as its
#'   \code{".missing"} prediction set and covariance/dispersion parameters
#'   fixed at the original fit's estimates.
#'
#' @noRd
get_kcv_local_glm_refit <- function(object, fold, local) {
  fold_data <- get_kcv_fold_ssn(object, fold)
  initial <- get_kcv_known_initials(object, glm = TRUE)
  do.call(
    ssn_glm,
    c(
      list(
        formula = object$formula, ssn.object = fold_data$ssn, family = object$family,
        tailup_type = initial$tailup_type, taildown_type = initial$taildown_type,
        euclid_type = initial$euclid_type, nugget_type = initial$nugget_type,
        tailup_initial = initial$tailup_initial, taildown_initial = initial$taildown_initial,
        euclid_initial = initial$euclid_initial, nugget_initial = initial$nugget_initial,
        dispersion_initial = initial$dispersion_initial, additive = object$additive,
        estmethod = object$estmethod, anisotropy = object$anisotropy,
        random = object$random, randcov_initial = initial$randcov_initial,
        partition_factor = object$partition_factor,
        local = get_kcv_estimation_local(object, fold, local), contrasts = object$contrasts
      )
    )
  )
}

#' Compute local (big-data) k-fold cross-validation predictions for a GLM
#'
#' For each fold, refits via \code{\link{get_kcv_local_glm_refit}()} and
#' predicts the held-out rows on the link scale, extracting/aligning the
#' result via \code{\link{extract_kcv_fold_prediction}()}.
#'
#' @param object A fitted \code{ssn_glm} model object.
#' @param fold_list A named list of per-fold row-index vectors.
#' @param se.fit Whether to also return prediction standard errors.
#' @param local A resolved \code{local} list.
#'
#' @return A list with \code{cv_predict_link}, \code{cv_predict} (response
#'   scale), and \code{se.fit} (link scale, \code{NULL} if not requested).
#'
#' @noRd
get_kcv_local.ssn_glm <- function(object, fold_list, se.fit, local) {
  pred_link <- numeric(object$n)
  if (se.fit) se <- numeric(object$n)
  for (i in seq_along(fold_list)) {
    held <- fold_list[[i]]
    refit <- get_kcv_local_glm_refit(object, held, local)
    prediction <- predict(
      refit, newdata = ".missing", type = "link", se.fit = se.fit,
      interval = "none", local = local
    )
    fold_prediction <- extract_kcv_fold_prediction(
      prediction, object$observed_index[held], se.fit
    )
    pred_link[held] <- fold_prediction$pred
    if (se.fit) se[held] <- fold_prediction$se.fit
  }
  list(
    cv_predict_link = pred_link,
    cv_predict = invlink(pred_link, object$family, object$size),
    se.fit = if (se.fit) se else NULL
  )
}
