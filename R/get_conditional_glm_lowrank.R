#' Prepare the base-sized GLM conditional approximation
#' @param object A fitted ssn_glm object.
#' @param index Observed rows selected for the simulation base.
#' @param covariance Optional prepared base covariance context.
#' @return Whitened conditional factors and a base-sized covariance context.
#' @noRd
get_conditional_glm_base <- function(object, index, covariance = NULL) {
  X <- model.matrix(object)[index, , drop = FALSE]
  eta <- fitted(object, type = "link")[index]
  offset <- model.offset(model.frame(object))
  w <- if (is.null(offset)) eta else eta - offset[index]
  if (is.null(covariance)) covariance <- get_conditional_lowrank_covariance(object, index)
  cov_lowchol_base <- t(chol(conditional_vecchia_covariance(covariance, seq_along(index))))
  SqrtSigInv_X <- forwardsolve(cov_lowchol_base, X)
  D <- get_D(object$family, eta, object$y[index], object$size[index],
    as.vector(coef(object, type = "dispersion")))
  # Congruence by the base factor avoids forming the spatial precision.
  mH <- diag(length(index)) - crossprod(cov_lowchol_base, D %*% cov_lowchol_base)
  cov_lowchol_mH <- tryCatch(t(chol(as.matrix(mH))), error = function(e) {
    stop("The simulation base latent precision is not positive definite; check the fitted model or select a different base.", call. = FALSE)
  })
  cov_betahat <- vcov(object)
  wts_latent <- backsolve(t(cov_lowchol_mH), forwardsolve(cov_lowchol_mH, SqrtSigInv_X))
  list(covariance = covariance, cov_lowchol_base = cov_lowchol_base,
    cov_lowchol_mH = cov_lowchol_mH, cov_betahat_lowchol = t(chol(cov_betahat)),
    SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_w = forwardsolve(cov_lowchol_base, w),
    wts_latent = wts_latent, betahat = coef(object))
}

#' Draw low-rank GLM simulations with cached factors
#' @param object A fitted ssn_glm object.
#' @param newdata Processed prediction data.
#' @param Xnew Aligned prediction design matrix.
#' @param local_list Resolved low-rank settings.
#' @param samples Number of draws.
#' @param latent Whether to retain full-base observed latent draws.
#' @param random_design Optional shared random-effect design.
#' @param covariance Optional prepared base covariance context.
#' @return Matrices of coefficients and offset-free prediction latent values.
#' @noRd
get_conditional_glm_lowrank <- function(object, newdata, Xnew, local_list, samples, latent = FALSE, random_design = NULL, covariance = NULL) {
  if (!is.null(object$random) && is.null(random_design)) random_design <- get_conditional_random_design(object, newdata, local_list$index$base)
  chunk_size <- get_conditional_chunk_size(local_list)
  base <- get_conditional_glm_base(object, local_list$index$base, covariance)
  effect <- if (!is.null(object$random)) {
    get_conditional_lowrank_effect(base$covariance, random_design)
  } else NULL
  n_base <- NROW(base$cov_lowchol_base)
  p <- length(base$betahat)
  beta <- matrix(NA_real_, p, samples, dimnames = list(names(base$betahat), NULL))
  residual <- matrix(NA_real_, n_base, samples)
  residual_mean <- base$SqrtSigInv_w - as.numeric(base$SqrtSigInv_X %*% base$betahat)
  residual_weights <- base$wts_latent - base$SqrtSigInv_X
  observed <- if (latent) matrix(NA_real_, NROW(model.matrix(object)), samples) else NULL
  effects <- if (!is.null(effect)) matrix(NA_real_, length(random_design$variance), samples) else NULL
  for (start in seq.int(1L, samples, by = chunk_size)) {
    columns <- seq.int(start, min(samples, start + chunk_size - 1L))
    n_draw <- length(columns)
    delta <- base$cov_betahat_lowchol %*% matrix(rnorm(p * n_draw), p, n_draw)
    beta[, columns] <- sweep(delta, 1L, base$betahat, "+")
    z <- backsolve(t(base$cov_lowchol_mH), matrix(rnorm(n_base * n_draw), n_base, n_draw))
    residual[, columns] <- sweep(residual_weights %*% delta + z, 1L, residual_mean, "+")
    if (latent) observed[local_list$index$base, columns] <- base$cov_lowchol_base %*% residual[, columns, drop = FALSE] +
      model.matrix(object)[local_list$index$base, , drop = FALSE] %*% beta[, columns, drop = FALSE]
    if (!is.null(effect)) {
      natural <- base$cov_lowchol_base %*% residual[, columns, drop = FALSE]
      a <- draw_conditional_random(effect, natural)
      effects[, columns] <- a
      residual[, columns] <- forwardsolve(t(effect$blocks[[1L]]$factor), natural - as.matrix(random_design$observed %*% a))
    }
  }
  if (!is.null(effect)) {
    base$covariance$randcov <- numeric()
    base$cov_lowchol_base <- t(effect$blocks[[1L]]$factor)
  }
  new_val <- draw_conditional_lowrank_blocks(base$covariance, newdata, Xnew, beta, residual,
    base$cov_lowchol_base, local_list, effects, if (!is.null(effect)) random_design$newdata else NULL)
  list(beta = beta, newdata = new_val, latent = observed)
}
