conditional_spatial_object <- function(object) {
  object["random"] <- list(NULL)
  object$coefficients$params_object$randcov <- NULL
  object
}

get_conditional_random_design <- function(object, newdata = NULL, index = NULL) {
  get_randcov_design(object$ssn.object$obs, object$coefficients$params_object$randcov,
                    object$partition_factor, newdata, index)
}

conditional_random_blocks <- function(object, size = 1000L, groups = NULL, method = "kmeans") {
  observed <- object$ssn.object$obs
  n <- NROW(observed)
  if (is.null(size)) size <- 1000L
  if (is.null(method)) method <- "kmeans"
  if (is.null(groups)) groups <- ceiling(n / size)
  if (groups == 1L) return(list(seq_len(n)))
  if (method == "kmeans") {
    xy <- sf::st_coordinates(observed)
    if (NROW(unique(xy)) < groups) stop("Too few distinct observed locations for the requested random-effect blocks; increase local$size_random or reduce local$groups_random, or use local$method_random = \"random\".", call. = FALSE)
    if (groups == n) return(as.list(seq_len(n)))
  }
  labels <- get_local_estimation_index(list(method = method, groups = groups, size = ceiling(n / groups)), observed, n)
  unname(split(seq_len(n), labels))
}

prepare_conditional_random <- function(object, design, local_list = NULL, covariance = NULL, spatial_factor = NULL, covariance_context = NULL) {
  Z <- design$observed
  n <- NROW(Z)
  active <- which(design$variance > 0 & design$info$conditioned)
  if (any(!is.finite(design$variance) | design$variance < 0)) stop("Invalid random-effect variances.", call. = FALSE)
  if (length(active)^2 > 1e8) stop("The random-effect effects covariance is too large; specify 'local' instead.", call. = FALSE)
  spatial <- conditional_spatial_object(object)
  if (!is.null(covariance_context)) covariance_context$randcov <- numeric()
  blocks <- list(seq_len(n))
  if (!is.null(local_list) && local_list$approximation == "vecchia") {
    blocks <- local_list$random_blocks
    if (is.null(blocks)) blocks <- conditional_random_blocks(object, local_list$size_random,
      local_list$groups_random, local_list$method_random)
  }
  precision <- diag(1 / design$variance[active], length(active))
  prepared <- vector("list", length(blocks))
  for (k in seq_along(blocks)) {
    rows <- blocks[[k]]
    if (!is.null(covariance_context)) {
      S <- conditional_vecchia_covariance(covariance_context, rows)
    } else if (is.null(covariance)) {
      sub <- spatial
      S <- get_decorrelate_observed_covariance(sub, spatial$ssn.object$obs[rows, , drop = FALSE])
    } else S <- covariance[rows, rows, drop = FALSE]
    R <- if (is.null(spatial_factor)) conditional_local_factor(S, "Random-effect spatial block")$U else spatial_factor
    Sinv <- chol2inv(R)
    cols <- which(as.numeric(Matrix::colSums(abs(Z[rows, active, drop = FALSE]))) > 0)
    Zp <- as.matrix(Z[rows, active[cols], drop = FALSE])
    weights <- Sinv %*% Zp
    precision[cols, cols] <- precision[cols, cols, drop = FALSE] + crossprod(Zp, weights)
    prepared[[k]] <- list(rows = rows, cols = cols, weights = weights,
      factor = if (length(blocks) == 1L) R else NULL)
  }
  # One shared effect covariance: H = (Gamma^-1 + sum Zp' Sp^-1 Zp)^-1.
  H <- if (length(active)) chol2inv(conditional_local_factor(precision, "Random-effect precision")$U) else matrix(numeric(), 0, 0)
  L <- if (length(active)) t(conditional_local_factor(H, "Random-effect conditional covariance")$U) else matrix(numeric(), 0, 0)
  list(design = design, active = active, blocks = prepared, covariance = H, factor = L, spatial = spatial)
}

draw_conditional_random <- function(context, residual) {
  q <- length(context$design$variance)
  samples <- NCOL(residual)
  z <- matrix(rnorm(q * samples), q, samples)
  a <- sqrt(context$design$variance) * z
  if (length(context$active)) {
    rhs <- matrix(0, length(context$active), samples)
    for (block in context$blocks) {
      rhs[block$cols, ] <- rhs[block$cols, , drop = FALSE] +
        crossprod(block$weights, residual[block$rows, , drop = FALSE])
    }
    a[context$active, ] <- context$covariance %*% rhs + context$factor %*% z[context$active, , drop = FALSE]
  }
  a
}

conditional_random_output <- function(a, design) {
  rownames(a) <- rownames(design$info)
  attr(a, "effect_info") <- design$info
  a
}

check_conditional_random_output <- function(value) {
  if (is.null(value)) stop("Requested random-effect draws were not generated.", call. = FALSE)
  value
}
