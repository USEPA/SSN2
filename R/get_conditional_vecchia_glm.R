get_conditional_vecchia_glm_context <- function(object, newdata_name, local_list, covariance = NULL) {
  X <- model.matrix(object)
  n <- NROW(X)
  eta <- fitted(object, type = "link")
  offset <- model.offset(model.frame(object))
  w <- as.numeric(if (is.null(offset)) eta else eta - offset)
  D <- diag(get_D(object$family, eta, object$y, object$size, as.vector(coef(object, type = "dispersion"))))
  if (any(!is.finite(w)) || any(!is.finite(D))) stop("Observed latent values and response curvature must be finite.", call. = FALSE)
  if (is.null(covariance)) covariance <- get_conditional_vecchia_covariance(object, newdata_name)
  observed <- covariance$observed
  if (isTRUE(object$anisotropy)) {
    params <- coef(object, type = "euclid")
    coords <- sf::st_coordinates(observed)
    xy <- transform_anis(coords[, 1L], coords[, 2L], params[["rotate"]], params[["scale"]])
    sf::st_geometry(observed) <- sf::st_sfc(lapply(seq_len(n), function(i) {
      sf::st_point(c(xy$xcoord_val[i], xy$ycoord_val[i]))
    }))
  }
  ord <- get_decorrelate_order(covariance$rows[seq_len(n), , drop = FALSE], local_list$ordering_observed, observed)
  operators <- vector("list", n)
  all_sites <- seq_len(n)
  for (k in seq_len(n)) {
    i <- ord[k]
    # N and K reuse the same observed-only scores without refilling overlaps.
    scores <- if (local_list$method != "all") conditional_vecchia_scores(covariance, i, local_list$method, all_sites) else NULL
    N <- conditional_vecchia_neighbors(covariance, i, ord[seq_len(k - 1L)], local_list$size, local_list$method, scores)
    K <- c(i, conditional_vecchia_neighbors(covariance, i, all_sites, local_list$size - 1L, local_list$method, scores))
    A <- unique(c(i, N, K))
    label <- paste("Observed site", i, "neighborhood")
    gaussian <- prepare_conditional_local_gaussian(conditional_vecchia_covariance(covariance, A),
      X[A, , drop = FALSE], w[A], D[A], label)
    op <- prepare_conditional_latent_operator(gaussian, 1L, match(N, A), label)
    op$i <- i
    op$neighbors <- N
    operators[[k]] <- op
  }
  spatial <- covariance
  spatial$randcov <- numeric()
  prediction <- if (NROW(covariance$newdata)) prepare_conditional_vecchia(spatial, local_list) else NULL
  list(operators = operators, order = ord, prediction = prediction, X = X,
    betahat = coef(object), coefficient_factor = t(chol(vcov(object))))
}

get_conditional_vecchia_glm <- function(object, newdata_name, Xnew, local_list, samples,
                                       latent = FALSE, random_design = NULL, random = FALSE, covariance = NULL) {
  chunk_size <- get_conditional_chunk_size(local_list)
  if (is.null(covariance)) covariance <- get_conditional_vecchia_covariance(object, newdata_name)
  if (!is.null(object$random) && is.null(random_design)) random_design <- conditional_vecchia_random_design(covariance)
  effect <- if (!is.null(object$random)) prepare_conditional_random(object, random_design, local_list,
    covariance_context = covariance) else NULL
  context <- get_conditional_vecchia_glm_context(object, newdata_name, local_list, covariance)
  n <- NROW(context$X)
  m <- NROW(Xnew)
  p <- length(context$betahat)
  beta <- matrix(NA_real_, p, samples, dimnames = list(names(context$betahat), NULL))
  new_val <- matrix(NA_real_, m, samples, dimnames = list(rownames(Xnew), NULL))
  observed <- if (latent) matrix(NA_real_, n, samples, dimnames = list(rownames(context$X), NULL)) else NULL
  effects <- if (random && !is.null(effect)) matrix(NA_real_, length(random_design$variance), samples) else NULL
  for (start in seq.int(1L, samples, by = chunk_size)) {
    columns <- seq.int(start, min(samples, start + chunk_size - 1L))
    count <- length(columns)
    delta <- context$coefficient_factor %*% matrix(rnorm(p * count), p, count)
    beta[, columns] <- sweep(delta, 1L, context$betahat, "+")
    w <- matrix(NA_real_, n, count)
    for (op in context$operators) {
      mu <- op$intercept + as.numeric(crossprod(op$coefficient, delta))
      if (length(op$neighbors)) mu <- mu + as.numeric(crossprod(op$weights, w[op$neighbors, , drop = FALSE]))
      w[op$i, ] <- mu + sqrt(op$variance) * rnorm(count)
    }
    if (latent) observed[, columns] <- w
    residual <- w - context$X %*% beta[, columns, drop = FALSE]
    if (!is.null(effect)) {
      a <- draw_conditional_random(effect, residual)
      if (random) effects[, columns] <- a
      residual <- residual - as.matrix(random_design$observed %*% a)
    }
    if (m) {
      residual <- draw_conditional_vecchia(context$prediction, residual, count)
      if (!is.null(effect)) residual <- residual + as.matrix(random_design$newdata %*% a)
      new_val[, columns] <- Xnew %*% beta[, columns, drop = FALSE] + residual
    }
  }
  list(beta = beta, newdata = new_val, latent = observed,
    random = if (!is.null(effects)) conditional_random_output(effects, random_design) else NULL)
}
