get_conditional_lowrank_covariance <- function(object, index, newdata_name = NULL) {
  object$ssn.object$obs <- object$ssn.object$obs[index, , drop = FALSE]
  object$n <- length(index)
  get_conditional_vecchia_covariance(object, newdata_name)
}

get_conditional_lowrank_effect <- function(covariance, design) {
  spatial <- covariance
  spatial$randcov <- numeric()
  S <- conditional_vecchia_covariance(spatial, seq_len(spatial$n))
  prepare_conditional_random(covariance$object, design, covariance = S)
}

get_conditional_glm_block_cov <- function(rows, covariance) {
  prediction <- covariance$n + rows
  list(cross = conditional_vecchia_covariance(covariance, seq_len(covariance$n), prediction),
       new = conditional_vecchia_covariance(covariance, prediction))
}

factor_conditional_glm_block <- function(covariance, cov_lowchol_base) {
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol_base, covariance$cross)
  cond_cov <- covariance$new - crossprod(SqrtSigInv_c0)
  list(SqrtSigInv_c0 = SqrtSigInv_c0, chol_cond_cov = t(chol(cond_cov)))
}

get_conditional_glm_block <- function(rows, covariance, cov_lowchol_base) {
  factor_conditional_glm_block(get_conditional_glm_block_cov(rows, covariance), cov_lowchol_base)
}

get_conditional_lm_lowrank <- function(object, newdata, Xnew, local_list, samples,
                                      covariance, random_design = NULL) {
  index <- local_list$index$base
  X <- model.matrix(object)[index, , drop = FALSE]
  y <- model.response(model.frame(object))
  offset <- model.offset(model.frame(object))
  if (!is.null(offset)) y <- y - offset
  y <- y[index]
  effect <- if (!is.null(random_design)) get_conditional_lowrank_effect(covariance, random_design) else NULL
  lower <- if (is.null(effect)) {
    t(chol(conditional_vecchia_covariance(covariance, seq_along(index))))
  } else t(effect$blocks[[1L]]$factor)
  if (!is.null(effect)) covariance$randcov <- numeric()
  betahat <- coef(object)
  factor <- t(chol(vcov(object)))
  beta <- matrix(NA_real_, length(betahat), samples, dimnames = list(names(betahat), NULL))
  residual <- matrix(NA_real_, length(index), samples)
  effects <- if (!is.null(effect)) matrix(NA_real_, length(random_design$variance), samples) else NULL
  chunk_size <- get_conditional_chunk_size(local_list)
  for (start in seq.int(1L, samples, by = chunk_size)) {
    columns <- seq.int(start, min(samples, start + chunk_size - 1L))
    beta[, columns] <- betahat + factor %*% matrix(rnorm(length(betahat) * length(columns)), length(betahat))
    natural <- y - X %*% beta[, columns, drop = FALSE]
    if (!is.null(effect)) {
      a <- draw_conditional_random(effect, natural)
      effects[, columns] <- a
      natural <- natural - as.matrix(random_design$observed %*% a)
    }
    residual[, columns] <- forwardsolve(lower, natural)
  }
  value <- draw_conditional_lowrank_blocks(covariance, newdata, Xnew, beta, residual,
    lower, local_list, effects, if (!is.null(effect)) random_design$newdata else NULL)
  list(beta = beta, newdata = value)
}

draw_conditional_lowrank <- function(object, newdata_name, local_list, samples, output,
                                     type = "link", newdata_size = NULL,
                                     center = if (inherits(object, "ssn_glm")) get_latent_center(object) else NULL) {
  if ("random" %in% output) {
    stop("Random-effect output is unavailable on low-rank paths; effects are conditioned only on the sampled base and used internally for prediction.", call. = FALSE)
  }
  if ("latent" %in% output && length(local_list$index$base) != NROW(model.matrix(object))) {
    stop("Full latent output is unavailable for a reduced low-rank base; only base latent values are drawn.", call. = FALSE)
  }
  get_conditional_chunk_size(local_list)
  if (!any(c("newdata", "latent") %in% output)) {
    return(draw_conditional_exact(object, NULL, samples, output, type, newdata_size, center = center))
  }
  glm <- inherits(object, "ssn_glm")
  context <- if (glm) get_conditional_context_glm(object, newdata_name, local = TRUE, center = center) else
    get_conditional_context(object, newdata_name, local = TRUE)
  covariance <- get_conditional_lowrank_covariance(object, local_list$index$base, newdata_name)
  design <- if (!is.null(object$random)) get_conditional_random_design(object, context$newdata, local_list$index$base) else NULL
  value <- if (glm) {
    get_conditional_glm_lowrank(object, context$newdata, context$x0, local_list, samples,
      latent = "latent" %in% output, random_design = design, covariance = covariance, center = center)
  } else {
    get_conditional_lm_lowrank(object, context$newdata, context$x0, local_list, samples,
      covariance, design)
  }
  if ("object" %in% output) value$object <- matrix(if (glm) context$w else context$y_object,
    nrow = NROW(context$Xmat), ncol = samples)
  if ("latent" %in% output) {
    if (!is.null(context$offset)) value$latent <- sweep(value$latent, 1L, context$offset, "+")
    rownames(value$latent) <- rownames(context$Xmat)
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
