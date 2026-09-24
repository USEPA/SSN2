
#' Validate and normalize the \code{local} argument for conditional simulation
#'
#' Dispatches on \code{local$approximation} to one of two big-data
#' approximations, matching spmodel's \code{conditional()}: \code{"low-rank"}
#' (the new default; see \code{\link{get_conditional_local_lowrank}()}) or
#' \code{"vecchia"} (SSN2's original sole engine; see
#' \code{\link{get_local_vecchia_settings}()}, shared with
#' \code{\link{get_ssn_simulate_local}()}).
#'
#' @param local A logical or list; see the \code{local} argument to
#'   \code{\link{conditional}()}.
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set being simulated.
#'
#' If \code{local} is omitted, it is set to \code{TRUE} (triggering the
#' \code{"low-rank"} default below) whenever the observed data size or the
#' number of prediction locations exceeds 5,000, matching spmodel's
#' \code{conditional()} (\code{get_local_list_conditional()},
#' \code{R/get_local_list.R}); otherwise it is set to \code{FALSE}.
#'
#' @return A resolved list with \code{method = "exact"} (the dense path), or
#'   \code{approximation = "vecchia"}/\code{"low-rank"} plus that
#'   approximation's own settings.
#'
#' @noRd
get_conditional_local <- function(local, object, newdata_name) {
  newdata <- object$ssn.object$preds[[newdata_name]]
  if (is.null(local)) {
    if (object$n > 5000 || NROW(newdata) > 5000) {
      local <- TRUE
      message(
        "Because the observed data size or the number of prediction locations exceeds 5,000, ",
        "we are using a low-rank big-data approximation. Set local = FALSE to use the exact ",
        "simulation, or local = list(approximation = \"vecchia\") for neighbor-truncated ",
        "sequential conditioning instead."
      )
    } else {
      local <- FALSE
    }
  }
  if (identical(local, FALSE)) {
    return(list(method = "exact"))
  }
  if (identical(local, TRUE)) local <- list()
  if (!is.list(local)) {
    stop("local must be FALSE, TRUE, or a list.", call. = FALSE)
  }
  check_local_approximation_ambiguity(local)
  if (is.null(local$approximation)) local$approximation <- "low-rank"
  if (length(local$approximation) != 1L || !local$approximation %in% c("low-rank", "vecchia")) {
    stop("local$approximation must be \"low-rank\" or \"vecchia\".", call. = FALSE)
  }
  if (identical(local$approximation, "vecchia")) {
    return(get_local_vecchia_settings(local))
  }
  get_conditional_local_lowrank(local, object, newdata, object$n, NROW(newdata))
}

#' Build the \code{"low-rank"} big data approximation settings for
#' \code{\link{conditional}()}
#'
#' Matches spmodel's \code{get_local_list_conditional_lowrank()}
#' (\code{R/get_local_list.R}) exactly: the observed-data base sample
#' (\code{method_base}/\code{size_base}/\code{reorder_base}) and the
#' \code{newdata} blocking (\code{method_new}/\code{size_new}/
#' \code{reorder_new}/\code{kmeans_new}) are independent settings, unlike
#' \code{\link{get_ssn_simulate_local_lowrank}()} (\code{R/ssn_simulate_local.R}),
#' which only blocks one side.
#'
#' @param local The partially-resolved \code{local} list (already has
#'   \code{approximation == "low-rank"}).
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata The prediction data frame for the requested prediction set.
#' @param n The observed sample size.
#' @param n_pred The number of \code{newdata} rows.
#'
#' @return \code{local}, with every \code{"low-rank"} default filled in and
#'   \code{index = list(base = ..., new = ...)} set (\code{index$new} is a
#'   plain row-index vector when \code{method_new == "all"}, or a list of
#'   blocks otherwise, matching spmodel's own convention).
#'
#' @noRd
get_conditional_local_lowrank <- function(local, object, newdata, n, n_pred) {
  names_local <- names(local)

  if (!"method_base" %in% names_local) local$method_base <- "base"
  if (length(local$method_base) != 1L || !local$method_base %in% c("base", "all")) {
    stop("local$method_base must be \"base\" or \"all\".", call. = FALSE)
  }
  if (!"method_new" %in% names_local) local$method_new <- "base"
  if (length(local$method_new) != 1L || !local$method_new %in% c("base", "all")) {
    stop("local$method_new must be \"base\" or \"all\".", call. = FALSE)
  }
  if (!"size_base" %in% names_local) local$size_base <- 5000L
  if (!"size_new" %in% names_local) local$size_new <- 1000L
  if (!"reorder_base" %in% names_local) local$reorder_base <- "grts"
  if (length(local$reorder_base) != 1L || !local$reorder_base %in% c("none", "random", "grts")) {
    stop("local$reorder_base must be \"none\", \"random\", or \"grts\".", call. = FALSE)
  }
  if (!"reorder_new" %in% names_local) local$reorder_new <- "random"
  if (length(local$reorder_new) != 1L || !local$reorder_new %in% c("none", "random")) {
    stop("local$reorder_new must be \"none\" or \"random\".", call. = FALSE)
  }
  if (!"kmeans_new" %in% names_local) {
    local$kmeans_new <- !identical(local$reorder_new, "none")
  }
  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$size_base >= n) local$method_base <- "all"
  if (local$size_new >= n_pred) local$method_new <- "all"

  index_base <- seq_len(n)
  index_new <- seq_len(n_pred)

  if (!identical(local$method_base, "all")) {
    observed <- object$ssn.object$obs
    rows_base <- get_decorrelate_rows(observed)
    index_base <- get_decorrelate_order(rows_base, local$reorder_base, observed)
    index_base <- index_base[seq_len(local$size_base)]
  }

  if (!identical(local$method_new, "all")) {
    rows_new <- get_decorrelate_rows(newdata)
    index_new <- get_decorrelate_order(rows_new, local$reorder_new, newdata)
    groups <- ceiling(n_pred / local$size_new)
    index_new <- get_local_lowrank_blocks(index_new, newdata, groups, local$kmeans_new)
  }

  local$index <- list(base = index_base, new = index_new)

  if (local$parallel) {
    n_index <- if (identical(local$method_new, "all")) 1L else length(local$index$new)
    local$ncores <- if ("ncores" %in% names_local) {
      min(n_index, local$ncores, parallel::detectCores())
    } else {
      min(n_index, parallel::detectCores())
    }
  }

  local
}

#' Simulate one block of a low-rank base-and-block conditional simulation
#'
#' Given observed residuals against simulated \code{beta} draws
#' (\code{base_residual}, restricted to the base sample's rows), draws
#' values at \code{block_index}'s \code{newdata} rows from their conditional
#' distribution given the base sample alone -- the workhorse behind
#' \code{\link{get_conditional_lowrank_ssn}()}, mirroring spmodel's
#' \code{get_conditional_new_from_base_adjust()}/
#' \code{get_conditional_new_from_base_adjust_glm()}
#' (\code{R/get_conditional_new_from_base.R}), simplified to materialize the
#' full residual matrix up front (matching
#' \code{\link{get_conditional_vecchia_ssn}()}'s own convention) rather than
#' spmodel's factored \code{SqrtSigInv_y}/\code{SqrtSigInv_X} optimization.
#' Takes every input as an explicit argument (rather than closing over them)
#' so it can be dispatched via \code{parallel::parLapply()} without requiring
#' \code{clusterExport()}.
#'
#' @param block_index Row indices (into \code{newdata}) of this block's
#'   locations.
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set being simulated.
#' @param base_index Row indices (into the observed data) of the base
#'   sample.
#' @param base_residual An \code{n_base x samples} matrix of observed
#'   residuals against simulated \code{beta} draws, already restricted to
#'   \code{base_index}'s rows.
#' @param base_lowchol The lower triangular Cholesky factor of the base
#'   sample's covariance matrix.
#' @param samples The number of simulations (columns of \code{base_residual}).
#' @param var_adj_pieces GLM latent-process uncertainty pieces from
#'   \code{\link{get_conditional_local_var_adj_pieces}()}, or \code{NULL} for
#'   Gaussian models.
#'
#' @return An \code{length(block_index) x samples} matrix of simulated
#'   residuals for this block.
#'
#' @noRd
get_conditional_new_from_base_ssn <- function(block_index, object, newdata_name, base_index,
                                               base_residual, base_lowchol, samples,
                                               var_adj_pieces = NULL) {
  n_block <- length(block_index)
  n_base <- length(base_index)

  cross_full <- get_block_obs_covariance(object, newdata_name, block_index)
  cross_covariance <- t(cross_full[, base_index, drop = FALSE])

  block_covariance <- get_block_pred_covariance(object, newdata_name, block_index, block_index)
  if (!is.null(var_adj_pieces)) {
    sqrt_pred <- var_adj_pieces$sqrt_mhinv_wts[, block_index, drop = FALSE]
    block_covariance <- block_covariance + crossprod(sqrt_pred, sqrt_pred)
  }

  sqrt_siginv_cross <- forwardsolve(base_lowchol, cross_covariance)
  sqrt_siginv_base <- forwardsolve(base_lowchol, base_residual)

  cond_cov <- block_covariance - crossprod(sqrt_siginv_cross, sqrt_siginv_cross)
  cond_cov <- as.matrix(Matrix::forceSymmetric(cond_cov))
  cond_lowchol <- chol_lower_with_pivot_fallback(
    cond_cov,
    "The block conditional covariance matrix for local low-rank conditional simulation is not numerically positive semidefinite; check covariance parameters and duplicate locations."
  )

  z <- matrix(rnorm(n_block * samples), n_block, samples)
  cond_mu <- crossprod(sqrt_siginv_cross, sqrt_siginv_base)
  (cond_lowchol %*% z) + cond_mu
}

#' Simulate conditional Gaussian draws via a low-rank base+block approximation
#'
#' Mirrors spmodel's low-rank branch of \code{conditional.splm()}/
#' \code{conditional.spglm()} (\code{R/conditional.R}): restricts the
#' observed data conditioned on to a base sample (unless
#' \code{local$method_base == "all"}), factors its covariance once, and
#' simulates every \code{newdata} block conditionally on the base sample
#' alone via \code{\link{get_conditional_new_from_base_ssn}()} (the whole of
#' \code{newdata} is treated as one block when
#' \code{local$method_new == "all"}, matching spmodel's own convention).
#' Unlike \code{\link{get_conditional_vecchia_ssn}()}, blocks do not depend
#' on each other or on the base's own simulated values (only on the fixed
#' \code{base_residual}), so \code{local$parallel} genuinely parallelizes
#' this engine.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set being simulated.
#' @param newdata The prediction data frame for \code{newdata_name}.
#' @param base_val An \code{n_obs x samples} matrix of composition-sampled
#'   observed-data residuals (one column per draw of \code{beta}).
#' @param local A resolved \code{local} list (from
#'   \code{\link{get_conditional_local_lowrank}()}).
#' @param samples The number of simulated columns to draw.
#' @param var_adj_pieces GLM latent-process uncertainty pieces from
#'   \code{\link{get_conditional_local_var_adj_pieces}()}, or \code{NULL} for
#'   Gaussian models.
#'
#' @return An \code{n_new x samples} matrix of simulated values, in
#'   \code{newdata}'s original row order.
#'
#' @noRd
get_conditional_lowrank_ssn <- function(object, newdata_name, newdata, base_val, local, samples,
                                         var_adj_pieces = NULL) {
  n_obs <- object$n
  n_new <- NROW(newdata)

  base_index <- if (identical(local$method_base, "all")) seq_len(n_obs) else local$index$base
  base_data <- object$ssn.object$obs[base_index, , drop = FALSE]
  base_residual <- base_val[base_index, , drop = FALSE]

  base_covariance <- get_decorrelate_observed_covariance(object, base_data)
  base_lowchol <- chol_lower_with_pivot_fallback(
    base_covariance,
    "The base sample covariance matrix for local low-rank conditional simulation is not numerically positive semidefinite; check covariance parameters and duplicate locations."
  )

  blocks <- if (identical(local$method_new, "all")) list(seq_len(n_new)) else local$index$new

  if (local$parallel) {
    cl <- parallel::makeCluster(local$ncores)
    new_val <- parallel::parLapply(
      cl, blocks, get_conditional_new_from_base_ssn,
      object, newdata_name, base_index, base_residual, base_lowchol, samples, var_adj_pieces
    )
    parallel::stopCluster(cl)
  } else {
    new_val <- lapply(
      blocks, get_conditional_new_from_base_ssn,
      object, newdata_name, base_index, base_residual, base_lowchol, samples, var_adj_pieces
    )
  }

  Y <- matrix(NA_real_, n_new, samples)
  for (i in seq_along(blocks)) {
    Y[blocks[[i]], ] <- as.matrix(new_val[[i]])
  }
  Y
}

#' Simulate conditional Gaussian draws via sequential covariance-neighbor conditioning
#'
#' A Vecchia-style sequential engine for conditional (kriging) simulation:
#' each \code{newdata} row is drawn, in \code{conditioning$ordering} order,
#' from its conditional distribution given (up to) \code{conditioning$size}
#' of the most-correlated points from a growing pool that starts with every
#' observed residual (\code{base_val}) and gains one simulated \code{newdata}
#' row per step. \code{conditioning$method = "all"} conditions on the entire
#' pool, recovering the exact conditional distribution. When
#' \code{var_adj_pieces} is supplied (GLM latent-process uncertainty), its
#' contribution is folded into the target/pool covariance blocks for the
#' \code{newdata}-vs-\code{newdata} portion only.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set being simulated.
#' @param newdata The prediction data frame for \code{newdata_name}.
#' @param base_val An \code{n_obs x samples} matrix of composition-sampled
#'   observed-data residuals (one column per draw of \code{beta}).
#' @param conditioning A resolved \code{"vecchia"} conditioning list (from
#'   \code{\link{get_conditional_local}()}/\code{\link{get_local_vecchia_settings}()})
#'   with \code{method}, \code{size}, and \code{ordering}.
#' @param samples The number of simulated columns to draw.
#' @param var_adj_pieces GLM latent-process uncertainty pieces from
#'   \code{\link{get_conditional_local_var_adj_pieces}()}, or \code{NULL} for
#'   Gaussian models.
#'
#' @return An \code{n_new x samples} matrix of simulated values, in
#'   \code{newdata}'s original row order.
#'
#' @noRd
get_conditional_vecchia_ssn <- function(object, newdata_name, newdata, base_val, conditioning, samples,
                                         var_adj_pieces = NULL) {
  n_obs <- object$n
  n_new <- NROW(newdata)
  observed <- object$ssn.object$obs

  rows_new <- get_decorrelate_rows(newdata)
  ord <- get_decorrelate_order(rows_new, conditioning$ordering, newdata)
  size <- if (identical(conditioning$method, "all")) (n_obs + n_new) else conditioning$size

  Y_ordered <- matrix(NA_real_, n_new, samples)
  Z <- matrix(rnorm(n_new * samples), n_new, samples)

  pool_val_full <- matrix(NA_real_, n_obs + n_new, samples)
  pool_val_full[seq_len(n_obs), ] <- base_val

  for (k in seq_len(n_new)) {
    current <- ord[[k]]
    n_new_pool <- k - 1L
    pool_n <- n_obs + n_new_pool
    pool_is_obs <- c(rep(TRUE, n_obs), rep(FALSE, n_new_pool))
    pool_idx <- c(seq_len(n_obs), if (n_new_pool > 0) ord[seq_len(n_new_pool)] else integer(0))

    obs_cross <- as.numeric(get_block_obs_covariance(object, newdata_name, current))
    pred_cross <- if (n_new_pool > 0) {
      as.numeric(get_block_pred_covariance(object, newdata_name, current, ord[seq_len(n_new_pool)]))
    } else {
      numeric(0)
    }
    cov_target_pool_full <- c(obs_cross, pred_cross)

    variance <- get_decorrelate_marginal_variance(object, newdata[current, , drop = FALSE])
    if (!is.null(var_adj_pieces)) {
      variance <- variance + var_adj_pieces$var_adj_diag[[current]]
    }

    keep <- get_decorrelate_covariance_neighbors(cov_target_pool_full, size)

    is_obs_keep <- pool_is_obs[keep]
    obs_group <- keep[is_obs_keep]
    pred_group <- keep[!is_obs_keep]
    keep_ordered <- c(obs_group, pred_group)
    cov_target_pool <- cov_target_pool_full[keep_ordered]

    obs_abs <- pool_idx[obs_group]
    pred_abs <- pool_idx[pred_group]
    a <- length(obs_abs)
    b <- length(pred_abs)

    cov_pool_pool <- matrix(0, a + b, a + b)
    if (a > 0) {
      cov_pool_pool[seq_len(a), seq_len(a)] <- get_decorrelate_observed_covariance(object, observed[obs_abs, , drop = FALSE])
    }
    if (b > 0) {
      cov_pool_pool[(a + 1):(a + b), (a + 1):(a + b)] <- get_block_pred_covariance(object, newdata_name, pred_abs, pred_abs)
    }
    if (a > 0 && b > 0) {
      # get_block_obs_covariance() returns length(pred_abs) x n_obs (pred x
      # obs); select obs_abs's columns then transpose to get the a x b
      # (obs x pred) block cov_pool_pool needs.
      op <- t(get_block_obs_covariance(object, newdata_name, pred_abs)[, obs_abs, drop = FALSE])
      cov_pool_pool[seq_len(a), (a + 1):(a + b)] <- op
      cov_pool_pool[(a + 1):(a + b), seq_len(a)] <- t(op)
    }

    if (!is.null(var_adj_pieces) && b > 0) {
      sqrt_target <- var_adj_pieces$sqrt_mhinv_wts[, current, drop = FALSE]
      sqrt_pred <- var_adj_pieces$sqrt_mhinv_wts[, pred_abs, drop = FALSE]
      cov_target_pool[(a + 1):(a + b)] <- cov_target_pool[(a + 1):(a + b)] + as.numeric(crossprod(sqrt_target, sqrt_pred))
      cov_pool_pool[(a + 1):(a + b), (a + 1):(a + b)] <- cov_pool_pool[(a + 1):(a + b), (a + 1):(a + b)] + crossprod(sqrt_pred, sqrt_pred)
    }

    cov_pool_pool <- as.matrix(Matrix::forceSymmetric(cov_pool_pool))
    chol_pool <- chol_lower_with_pivot_fallback(
      cov_pool_pool,
      paste0(
        "The neighbor covariance matrix for simulated prediction record ", k,
        " (of newdata \"", newdata_name, "\") is not numerically positive semidefinite; ",
        "conditional simulation cannot proceed for this prediction set."
      )
    )

    w <- backsolve(t(chol_pool), forwardsolve(chol_pool, cov_target_pool))
    cond_var <- max(variance - sum(w * cov_target_pool), 0)
    pool_val <- pool_val_full[keep_ordered, , drop = FALSE]
    cond_mean <- as.numeric(crossprod(w, pool_val))

    Y_ordered[k, ] <- cond_mean + sqrt(cond_var) * Z[k, ]
    pool_val_full[n_obs + k, ] <- Y_ordered[k, ]
  }

  Y <- matrix(NA_real_, n_new, samples)
  Y[ord, ] <- Y_ordered
  Y
}

#' Compute the diagonal GLM latent-process uncertainty adjustment for local conditional simulation
#'
#' Ports spmodel's \code{get_conditional_vecchia_glm()} design: the GLM
#' latent-process ("var_adj") uncertainty term is a single joint,
#' non-truncatable distribution over all observed data, so it is computed
#' once, densely, over the observed data here (never subject to
#' \code{conditioning$size} truncation), and its per-\code{newdata}-row
#' contribution is folded into \code{\link{get_conditional_vecchia_ssn}()}'s
#' target/pool covariance blocks only for the \code{newdata}-vs-\code{newdata}
#' portion.
#'
#' @param context A GLM conditional-simulation context from
#'   \code{\link{get_conditional_context_glm}()}, with \code{object},
#'   \code{Xmat}, \code{cov_betahat_uncorrected}, \code{family}, \code{w},
#'   \code{y}, \code{size}, \code{dispersion}, \code{newdata_name},
#'   \code{x0}.
#'
#' @return A list with \code{sqrt_mhinv_wts} (an \code{n_obs x n_new} matrix)
#'   and \code{var_adj_diag} (a length-\code{n_new} vector), the per-row
#'   diagonal latent-process variance contribution.
#'
#' @noRd
get_conditional_local_var_adj_pieces <- function(context) {
  object <- context$object

  cov_matrix_val <- covmatrix(object)
  cov_lowchol_base <- t(chol(cov_matrix_val))
  SigInv <- chol2inv(t(cov_lowchol_base))
  SigInv_X <- SigInv %*% context$Xmat
  wts_beta <- tcrossprod(context$cov_betahat_uncorrected, SigInv_X)
  Ptheta <- SigInv - SigInv_X %*% wts_beta # n_obs x n_obs

  D <- get_D(context$family, context$w, context$y, context$size, context$dispersion)
  H <- as.matrix(D - Ptheta)
  cov_lowchol_mH <- chol_lower_with_pivot_fallback(
    -H,
    paste0(
      "The negative Hessian of the latent process for \"", context$newdata_name,
      "\" is not numerically positive semidefinite; local conditional simulation cannot proceed for this prediction set."
    )
  )

  C0_all <- covmatrix(object, context$newdata_name, cov_type = "obs.pred") # n_obs x n_new
  c0_mat <- t(C0_all) # n_new x n_obs
  wts_pred_all <- context$x0 %*% wts_beta + c0_mat %*% SigInv - (c0_mat %*% SigInv_X) %*% wts_beta # n_new x n_obs
  wts_pred_all <- t(wts_pred_all) # n_obs x n_new
  sqrt_mhinv_wts <- forwardsolve(cov_lowchol_mH, wts_pred_all) # n_obs x n_new
  var_adj_diag <- colSums(sqrt_mhinv_wts^2)

  list(sqrt_mhinv_wts = sqrt_mhinv_wts, var_adj_diag = var_adj_diag)
}

#' Draw local (big-data) conditional Gaussian samples
#'
#' Local-conditioning analogue of \code{\link{draw_conditional_gaussian}()}:
#' composition-samples \code{beta} the same way, but replaces the exact dense
#' kriging correction with either \code{\link{get_conditional_lowrank_ssn}()}'s
#' base+block approximation or \code{\link{get_conditional_vecchia_ssn}()}'s
#' sequential, covariance-neighbor-truncated simulation, per
#' \code{conditioning$approximation}.
#'
#' @param context A Gaussian conditional-simulation context from
#'   \code{\link{get_conditional_context}()}.
#' @param conditioning A resolved conditioning list (from
#'   \code{\link{get_conditional_local}()}).
#' @param samples The number of simulated columns to draw.
#' @param output A character vector of requested outputs; any of
#'   \code{"object"}, \code{"beta"}, \code{"newdata"}.
#'
#' @return A list with elements named by \code{output}, each an
#'   \code{n x samples} matrix.
#'
#' @noRd
draw_conditional_gaussian_local <- function(context, conditioning, samples, output) {
  p <- NCOL(context$Xmat)

  val <- list()
  if ("object" %in% output) {
    val$object <- matrix(rep(context$y_object, times = samples), ncol = samples)
  }

  need_beta <- any(c("beta", "newdata") %in% output)
  if (!need_beta) {
    return(val)
  }

  cov_betahat_lowchol <- t(chol(context$cov_betahat))
  beta_draws <- as.vector(context$betahat) + cov_betahat_lowchol %*% matrix(rnorm(p * samples), p, samples)
  rownames(beta_draws) <- names(context$betahat)

  if ("beta" %in% output) {
    val$beta <- beta_draws
  }

  if ("newdata" %in% output) {
    base_val <- as.vector(context$y) - context$Xmat %*% beta_draws
    kriging_correction <- if (identical(conditioning$approximation, "low-rank")) {
      get_conditional_lowrank_ssn(
        context$object, context$newdata_name, context$newdata, base_val, conditioning, samples
      )
    } else {
      get_conditional_vecchia_ssn(
        context$object, context$newdata_name, context$newdata, base_val, conditioning, samples
      )
    }
    draws <- context$x0 %*% beta_draws + kriging_correction

    if (!is.null(context$newdata_offset)) {
      draws <- draws + context$newdata_offset
    }
    if (context$add_newdata_rows) {
      rownames(draws) <- context$object$missing_index
    }
    val$newdata <- draws
  }

  val
}

#' Draw local (big-data) conditional GLM samples
#'
#' Local-conditioning analogue of \code{\link{draw_conditional_glm}()}:
#' composition-samples \code{beta} the same way, folds in the GLM
#' latent-process uncertainty adjustment via
#' \code{\link{get_conditional_local_var_adj_pieces}()}, and replaces the
#' exact dense kriging correction with either
#' \code{\link{get_conditional_lowrank_ssn}()}'s base+block approximation or
#' \code{\link{get_conditional_vecchia_ssn}()}'s sequential,
#' covariance-neighbor-truncated simulation (per
#' \code{conditioning$approximation}) on the link scale before applying the
#' requested \code{type} transform.
#'
#' @param context A GLM conditional-simulation context from
#'   \code{\link{get_conditional_context_glm}()}.
#' @param conditioning A resolved conditioning list (from
#'   \code{\link{get_conditional_local}()}).
#' @param samples The number of simulated columns to draw.
#' @param type One of \code{"link"}, \code{"response"}, or \code{"new"};
#'   see \code{\link{conditional}()}.
#' @param output A character vector of requested outputs; any of
#'   \code{"object"}, \code{"beta"}, \code{"newdata"}.
#' @param newdata_size The binomial trial size used when \code{type} implies
#'   response-family draws.
#'
#' @return A list with elements named by \code{output}, each an
#'   \code{n x samples} matrix.
#'
#' @noRd
draw_conditional_glm_local <- function(context, conditioning, samples, type, output, newdata_size) {
  p <- NCOL(context$Xmat)

  val <- list()
  if ("object" %in% output) {
    val$object <- matrix(rep(context$w, times = samples), ncol = samples)
  }

  need_beta <- any(c("beta", "newdata") %in% output)
  if (!need_beta) {
    return(val)
  }

  cov_betahat_lowchol <- t(chol(context$cov_betahat))
  beta_draws <- as.vector(context$betahat) + cov_betahat_lowchol %*% matrix(rnorm(p * samples), p, samples)
  rownames(beta_draws) <- names(context$betahat)

  if ("beta" %in% output) {
    val$beta <- beta_draws
  }

  if ("newdata" %in% output) {
    base_val <- as.vector(context$w_free) - context$Xmat %*% beta_draws
    var_adj_pieces <- get_conditional_local_var_adj_pieces(context)

    kriging_correction <- if (identical(conditioning$approximation, "low-rank")) {
      get_conditional_lowrank_ssn(
        context$object, context$newdata_name, context$newdata, base_val, conditioning, samples,
        var_adj_pieces = var_adj_pieces
      )
    } else {
      get_conditional_vecchia_ssn(
        context$object, context$newdata_name, context$newdata, base_val, conditioning, samples,
        var_adj_pieces = var_adj_pieces
      )
    }
    link_draws <- context$x0 %*% beta_draws + kriging_correction

    if (!is.null(context$newdata_offset)) {
      link_draws <- link_draws + context$newdata_offset
    }

    draws <- if (identical(type, "link")) {
      link_draws
    } else {
      mu <- invlink(link_draws, context$family, size = 1)
      if (identical(type, "response")) {
        if (identical(context$family, "binomial")) mu * newdata_size else mu
      } else {
        draw_glm_response(context$family, mu, context$dispersion, newdata_size)
      }
    }

    if (context$add_newdata_rows) {
      rownames(draws) <- context$object$missing_index
    }
    val$newdata <- draws
  }

  val
}
