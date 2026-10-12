prepare_conditional_vecchia <- function(context, local_list) {
  n <- context$n
  ord <- local_list$order
  pool <- c(seq_len(n), n + ord)
  # Prediction ties follow the reordered sequence, not original prediction IDs.
  context$tie_index <- match(seq_along(pool), pool)
  operators <- vector("list", length(ord))
  for (k in seq_along(ord)) {
    i <- n + ord[k]
    candidates <- pool[seq_len(n + k - 1L)]
    neighbors <- conditional_vecchia_neighbors(context, i, candidates, local_list$size, local_list$method)
    covariance <- conditional_vecchia_covariance(context, c(i, neighbors))
    op <- prepare_conditional_local_weights(covariance, 1L, seq_along(neighbors) + 1L,
      paste("Prediction site", ord[k]))
    op$neighbors <- match(neighbors, pool)
    operators[[k]] <- op
  }
  list(operators = operators, order = ord, n = n)
}

draw_conditional_vecchia <- function(context, base_val, samples) {
  m <- length(context$order)
  pool <- matrix(NA_real_, context$n + m, samples)
  pool[seq_len(context$n), ] <- base_val
  Z <- matrix(rnorm(m * samples), m, samples)
  for (k in seq_len(m)) {
    op <- context$operators[[k]]
    pool[context$n + k, ] <- as.numeric(crossprod(op$weights, pool[op$neighbors, , drop = FALSE])) +
      sqrt(op$variance) * Z[k, ]
  }
  value <- matrix(NA_real_, m, samples)
  value[context$order, ] <- pool[context$n + seq_len(m), , drop = FALSE]
  value
}
