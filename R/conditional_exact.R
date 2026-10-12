prepare_conditional_exact_prediction <- function(object, newdata_name, lower = NULL) {
  if (is.null(lower)) lower <- t(chol(as.matrix(covmatrix(object))))
  C0 <- as.matrix(covmatrix(object, newdata_name, cov_type = "obs.pred"))
  whitened <- forwardsolve(lower, C0)
  context <- list(Sigma22 = as.matrix(covmatrix(object, newdata_name, cov_type = "pred.pred")),
    SqrtSigInv_C0 = whitened, newdata_name = newdata_name)
  conditional <- get_conditional_cov(context)
  list(weights = t(backsolve(t(lower), whitened)),
       factor = conditional$chol_cond_cov, covariance = conditional$Sigma_cond)
}

draw_conditional_exact_prediction <- function(prediction, residual) {
  m <- NROW(prediction$weights)
  samples <- NCOL(residual)
  prediction$weights %*% residual +
    prediction$factor %*% matrix(rnorm(m * samples), m, samples)
}

draw_conditional_exact <- function(object, newdata_name, samples, output,
                                   type = "link", newdata_size = NULL) {
  glm <- inherits(object, "ssn_glm")
  context <- if (glm) get_conditional_context_glm(object, newdata_name, local = TRUE) else
    get_conditional_context(object, newdata_name, local = TRUE)
  val <- list()
  if ("object" %in% output) {
    val$object <- matrix(if (glm) context$w else context$y_object,
                         nrow = NROW(context$Xmat), ncol = samples)
  }
  if (!any(c("newdata", "beta", "latent", "random") %in% output)) return(val)
  need_field <- any(c("newdata", "latent", "random") %in% output)
  joint <- if (glm && need_field) get_conditional_glm_joint(object) else NULL
  effect <- prediction <- NULL
  if (need_field && !is.null(object$random) && any(c("newdata", "random") %in% output)) {
    design <- get_conditional_random_design(object, context$newdata)
    effect <- prepare_conditional_random(object, design)
  }
  if ("newdata" %in% output) {
    spatial <- if (is.null(effect)) object else effect$spatial
    lower <- if (!is.null(effect)) t(effect$blocks[[1L]]$factor) else
      if (!is.null(joint)) joint$cov_lowchol else NULL
    prediction <- prepare_conditional_exact_prediction(spatial, newdata_name, lower)
  }
  if (!is.null(joint)) {
    draws <- draw_conditional_glm_joint(joint, samples, residual = TRUE)
    beta <- draws$beta
    residual <- draws$residual
  } else {
    p <- length(context$betahat)
    beta <- sweep(t(chol(context$cov_betahat)) %*% matrix(rnorm(p * samples), p, samples),
                  1L, context$betahat, "+")
    rownames(beta) <- names(context$betahat)
    residual <- if (need_field) context$y - context$Xmat %*% beta else NULL
  }
  if ("beta" %in% output) val$beta <- beta
  if ("latent" %in% output) {
    val$latent <- residual + context$Xmat %*% beta
    if (!is.null(context$offset)) val$latent <- sweep(val$latent, 1L, context$offset, "+")
    rownames(val$latent) <- rownames(context$Xmat)
  }
  if (!is.null(effect)) {
    a <- draw_conditional_random(effect, residual)
    residual <- residual - as.matrix(effect$design$observed %*% a)
    if ("random" %in% output) val$random <- conditional_random_output(a, effect$design)
  }
  if ("newdata" %in% output) {
    value <- context$x0 %*% beta + draw_conditional_exact_prediction(prediction, residual)
    if (!is.null(effect)) value <- value + as.matrix(effect$design$newdata %*% a)
    if (!is.null(context$newdata_offset)) value <- sweep(value, 1L, context$newdata_offset, "+")
    if (glm && type != "link") {
      if (identical(context$family, "binomial") && is.null(newdata_size)) newdata_size <- rep(1, NROW(value))
      mu <- invlink(value, context$family, size = 1)
      value <- if (type == "response") {
        if (identical(context$family, "binomial")) mu * newdata_size else mu
      } else draw_glm_response(context$family, mu, context$dispersion, newdata_size)
    }
    if (context$add_newdata_rows) rownames(value) <- object$missing_index
    val$newdata <- value
  }
  val
}
