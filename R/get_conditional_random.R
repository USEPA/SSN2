conditional_spatial_object <- function(object) {
  object["random"] <- list(NULL)
  object$coefficients$params_object$randcov <- NULL
  object
}

get_conditional_random_design <- function(object, newdata = NULL, index = NULL) {
  get_randcov_design(object$ssn.object$obs, object$coefficients$params_object$randcov,
                    object$partition_factor, newdata, index)
}

prepare_conditional_random <- function(object, design, covariance = NULL, spatial_factor = NULL) {
  Z <- design$observed
  n <- NROW(Z)
  active <- which(design$variance > 0 & design$info$conditioned)
  if (any(!is.finite(design$variance) | design$variance < 0)) stop("Invalid random-effect variances.", call. = FALSE)
  if (length(active)^2 > 1e8) warning(sprintf(
    "Conditional simulation requires factoring a covariance matrix for %s random-effect coefficients. This may result in exceedingly long computational times.",
    format(length(active), big.mark = ",", scientific = FALSE)
  ), call. = FALSE, immediate. = TRUE)
  spatial <- conditional_spatial_object(object)
  blocks <- list(seq_len(n))
  precision <- diag(1 / design$variance[active], length(active))
  prepared <- vector("list", length(blocks))
  for (k in seq_along(blocks)) {
    rows <- blocks[[k]]
    if (is.null(covariance)) {
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
