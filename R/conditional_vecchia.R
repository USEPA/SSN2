conditional_vecchia_random_design <- function(covariance) {
  list(observed = covariance$random[seq_len(covariance$n), , drop = FALSE],
    newdata = covariance$random[covariance$n + seq_len(NROW(covariance$newdata)), , drop = FALSE],
    variance = covariance$randcov, info = covariance$effect_info)
}

get_conditional_vecchia_lm <- function(object, context, covariance, local_list, samples, random = FALSE) {
  design <- if (!is.null(object$random)) conditional_vecchia_random_design(covariance) else NULL
  effect <- if (!is.null(design)) prepare_conditional_random(object, design, local_list,
    covariance_context = covariance) else NULL
  covariance$randcov <- numeric()
  prediction <- if (NROW(context$newdata)) prepare_conditional_vecchia(covariance, local_list) else NULL
  p <- length(context$betahat)
  factor <- t(chol(context$cov_betahat))
  beta <- matrix(NA_real_, p, samples, dimnames = list(names(context$betahat), NULL))
  value <- matrix(NA_real_, NROW(context$x0), samples, dimnames = list(rownames(context$x0), NULL))
  effects <- if (random && !is.null(effect)) matrix(NA_real_, length(design$variance), samples) else NULL
  chunk_size <- get_conditional_chunk_size(local_list)
  for (start in seq.int(1L, samples, by = chunk_size)) {
    columns <- seq.int(start, min(samples, start + chunk_size - 1L))
    count <- length(columns)
    beta[, columns] <- sweep(factor %*% matrix(rnorm(p * count), p, count), 1L, context$betahat, "+")
    residual <- context$y - context$Xmat %*% beta[, columns, drop = FALSE]
    if (!is.null(effect)) {
      a <- draw_conditional_random(effect, residual)
      if (random) effects[, columns] <- a
      residual <- residual - as.matrix(design$observed %*% a)
    }
    if (!is.null(prediction)) {
      residual <- draw_conditional_vecchia(prediction, residual, count)
      if (!is.null(effect)) residual <- residual + as.matrix(design$newdata %*% a)
      value[, columns] <- context$x0 %*% beta[, columns, drop = FALSE] + residual
    }
  }
  list(beta = beta, newdata = value,
    random = if (!is.null(effects)) conditional_random_output(effects, design) else NULL)
}

draw_conditional_vecchia_model <- function(object, newdata_name, local_list, samples, output,
                                           type = "link", newdata_size = NULL) {
  if (!any(c("newdata", "latent", "random") %in% output)) {
    return(draw_conditional_exact(object, NULL, samples, output, type, newdata_size))
  }
  glm <- inherits(object, "ssn_glm")
  context <- if (glm) get_conditional_context_glm(object, newdata_name, local = TRUE) else
    get_conditional_context(object, newdata_name, local = TRUE)
  covariance <- get_conditional_vecchia_covariance(object, newdata_name)
  value <- if (glm) {
    get_conditional_vecchia_glm(object, newdata_name, context$x0, local_list, samples,
      latent = "latent" %in% output, random = "random" %in% output, covariance = covariance)
  } else get_conditional_vecchia_lm(object, context, covariance, local_list, samples, random = "random" %in% output)
  if ("object" %in% output) value$object <- matrix(if (glm) context$w else context$y_object,
    nrow = NROW(context$Xmat), ncol = samples)
  if ("latent" %in% output && !is.null(context$offset)) {
    value$latent <- sweep(value$latent, 1L, context$offset, "+")
  }
  if ("newdata" %in% output) {
    if (!is.null(context$newdata_offset)) value$newdata <- sweep(value$newdata, 1L, context$newdata_offset, "+")
    if (glm && type != "link") {
      if (identical(context$family, "binomial") && is.null(newdata_size)) newdata_size <- rep(1, NROW(context$newdata))
      chunk_size <- get_conditional_chunk_size(local_list)
      for (start in seq.int(1L, samples, by = chunk_size)) {
        columns <- seq.int(start, min(samples, start + chunk_size - 1L))
        mu <- invlink(value$newdata[, columns, drop = FALSE], context$family, size = 1)
        value$newdata[, columns] <- if (type == "response") {
          if (identical(context$family, "binomial")) mu * newdata_size else mu
        } else draw_glm_response(context$family, mu, context$dispersion, newdata_size)
      }
    }
    if (context$add_newdata_rows) rownames(value$newdata) <- object$missing_index
  }
  value
}
