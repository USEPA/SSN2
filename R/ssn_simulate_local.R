#' Validate and normalize the shared \code{"vecchia"} sub-settings of \code{local}
#'
#' Shared by \code{\link{get_ssn_simulate_local}()} and \code{conditional()}'s
#' own local resolver: fills/validates \code{method}, \code{size}, and
#' \code{ordering} -- the sequential, neighbor-truncated big-data
#' approximation, reached via \code{local$approximation = "vecchia"} (or, for
#' \code{local = TRUE}, no longer the default -- see
#' \code{\link{get_ssn_simulate_local}()}).
#'
#' @param local A list already known to have \code{approximation == "vecchia"}
#'   (or about to be given it).
#' @param default_ordering The default \code{ordering} value.
#'
#' @return \code{local}, with \code{approximation = "vecchia"} and every
#'   \code{"vecchia"} setting filled in.
#'
#' @noRd
get_local_vecchia_settings <- function(local, default_ordering = "pid") {
  if (is.null(local$method)) {
    local$method <- "covariance"
  } else if (length(local$method) != 1L || !local$method %in% c("all", "covariance")) {
    stop("local$method must be \"all\" or \"covariance\".", call. = FALSE)
  }
  if (identical(local$method, "covariance")) {
    if (is.null(local$size)) {
      local$size <- 30L
    } else if (length(local$size) != 1L || !is.numeric(local$size) ||
        !is.finite(local$size) || local$size < 1) {
      stop("local$size must be one positive number for covariance conditioning.", call. = FALSE)
    }
    local$size <- as.integer(local$size)
  }
  if (is.null(local$ordering)) {
    local$ordering <- default_ordering
  } else if (length(local$ordering) != 1L || is.na(local$ordering) ||
      !local$ordering %in% c(
        "pid", "none", "random", "maxmin", "middleout", "outsidein",
        "coordinate", "grts"
      )) {
    stop(
      "local$ordering must be \"pid\", \"none\", \"random\", \"maxmin\", \"middleout\", ",
      "\"outsidein\", \"coordinate\", or \"grts\".",
      call. = FALSE
    )
  }
  local$approximation <- "vecchia"
  local
}

#' Reject a pre-\code{approximation} \code{local} list
#'
#' The \code{method}/\code{size} controls belong to \code{"vecchia"}.
#' Require an explicit approximation when either is supplied so the default
#' \code{"low-rank"} path cannot silently ignore them.
#'
#' @param local A list \code{local} argument.
#'
#' @return \code{NULL}, invisibly, if \code{local} is unambiguous; otherwise
#'   an error.
#'
#' @noRd
check_local_approximation_ambiguity <- function(local) {
  names_local <- names(local)
  if (!"approximation" %in% names_local && (("method" %in% names_local) || ("size" %in% names_local))) {
    stop(
      "local now defaults to a low-rank approximation and no longer accepts method/size directly. ",
      "Set local$approximation = \"vecchia\" to keep using neighbor-truncated sequential conditioning, ",
      "or remove method/size to use the new low-rank defaults.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Split a big-data approximation's remaining locations into blocks
#'
#' Shared by \code{\link{get_ssn_simulate_local_lowrank}()} and
#' \code{conditional()}'s own low-rank resolver: either \code{kmeans()}
#' clusters \code{index_new}'s own coordinates (reusing the same
#' \code{kmeans()}-on-\code{sf::st_coordinates()} pattern
#' \code{\link{get_local_estimation_index}()} (\code{R/local.R}) already
#' established for \code{ssn_lm()}/\code{ssn_glm()}'s own local-likelihood
#' fitting), or splits them into contiguous blocks of (nearly) equal size,
#' matching spmodel's \code{get_local_list_simulation_lowrank()}/
#' \code{get_local_list_conditional_lowrank()} exactly.
#'
#' @param index_new The (already ordered) row indices to split into blocks.
#' @param data The \code{sf} rows \code{index_new} indexes into. Only used
#'   (and only required) when \code{kmeans_new} is \code{TRUE}.
#' @param groups The number of blocks.
#' @param kmeans_new Whether to cluster on coordinates (\code{TRUE}) or split
#'   contiguously (\code{FALSE}).
#'
#' @return A list of \code{groups} row-index vectors, together covering every
#'   element of \code{index_new} exactly once.
#'
#' @noRd
get_local_lowrank_blocks <- function(index_new, data, groups, kmeans_new) {
  n_new <- length(index_new)
  if (kmeans_new) {
    coords <- sf::st_coordinates(data[index_new, , drop = FALSE])
    cluster <- stats::kmeans(coords, centers = groups, iter.max = 30)$cluster
    return(split(index_new, cluster))
  }
  sizes <- rep(n_new %/% groups, groups)
  extra <- n_new %% groups
  if (extra > 0L) sizes[seq_len(extra)] <- sizes[seq_len(extra)] + 1L
  split(index_new, rep(seq_len(groups), times = sizes))
}

#' Validate and normalize the \code{local} argument for unconditional simulation
#'
#' Dispatches on \code{local$approximation} to one of two big-data
#' approximations, matching spmodel's \code{sprnorm()}: \code{"low-rank"}
#' (the new default; see \code{\link{get_ssn_simulate_local_lowrank}()}) or
#' \code{"vecchia"} (SSN2's original sole engine; see
#' \code{\link{get_local_vecchia_settings}()}).
#'
#' @param local A logical or list; see the \code{local} argument to
#'   \code{\link{ssn_rnorm}()}/\code{\link{ssn_simulate}()}.
#' @param n The number of locations to simulate.
#' @param data The \code{sf} rows to simulate at (\code{ssn.object$obs}).
#'
#' @return A resolved list with \code{method = "exact"} (the dense path), or
#'   \code{approximation = "vecchia"}/\code{"low-rank"} plus that
#'   approximation's own settings.
#'
#' @noRd
get_ssn_simulate_local <- function(local, n, data) {
  if (is.null(local) || identical(local, FALSE)) {
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
  get_ssn_simulate_local_lowrank(local, n, data)
}

#' Build the \code{"low-rank"} big data approximation settings for
#' \code{\link{ssn_rnorm}()}/\code{\link{ssn_simulate}()}
#'
#' Matches spmodel's \code{get_local_list_simulation_lowrank()}
#' (\code{R/get_local_list.R}) exactly: a base sample is drawn (ordered via
#' \code{reorder_base}, sized \code{size_base}), and the remaining locations
#' are split into blocks of (approximately) \code{size_new}, either by
#' \code{kmeans_new} clustering on coordinates or contiguously (see
#' \code{\link{get_local_lowrank_blocks}()}).
#'
#' @param local The partially-resolved \code{local} list (already has
#'   \code{approximation == "low-rank"}).
#' @param n The number of locations to simulate.
#' @param data The \code{sf} rows to simulate at.
#'
#' @return \code{local}, with every \code{"low-rank"} default filled in and,
#'   when \code{method_base != "all"}, \code{index = list(base = ..., new =
#'   ...)} set.
#'
#' @noRd
get_ssn_simulate_local_lowrank <- function(local, n, data) {
  names_local <- names(local)

  if (!"method_base" %in% names_local) local$method_base <- "base"
  if (length(local$method_base) != 1L || !local$method_base %in% c("base", "all")) {
    stop("local$method_base must be \"base\" or \"all\".", call. = FALSE)
  }
  if (!"size_base" %in% names_local) local$size_base <- 5000L
  if (!"size_new" %in% names_local) local$size_new <- 1000L
  if (!"reorder_base" %in% names_local) local$reorder_base <- "grts"
  if (length(local$reorder_base) != 1L || !local$reorder_base %in% c("none", "random", "grts")) {
    stop("local$reorder_base must be \"none\", \"random\", or \"grts\".", call. = FALSE)
  }
  if (!"kmeans_new" %in% names_local) {
    local$kmeans_new <- !identical(local$reorder_base, "none")
  }
  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$size_base >= n) local$method_base <- "all"

  if (!identical(local$method_base, "all")) {
    rows <- get_decorrelate_rows(data)
    index <- get_decorrelate_order(rows, local$reorder_base, data)
    index_base <- index[seq_len(local$size_base)]
    index_new <- index[-seq_len(local$size_base)]
    groups <- ceiling(length(index_new) / local$size_new)

    local$index <- list(
      base = index_base,
      new = get_local_lowrank_blocks(index_new, data, groups, local$kmeans_new)
    )

    if (local$parallel) {
      n_index <- length(local$index$new)
      local$ncores <- if ("ncores" %in% names_local) {
        min(n_index, local$ncores, parallel::detectCores())
      } else {
        min(n_index, parallel::detectCores())
      }
    }
  }

  local
}

#' Build a lightweight covariance-fit-shaped list for unconditional simulation
#'
#' Assembles the minimal \code{covariance_fit}-shaped object (matching the
#' shape produced by \code{\link{get_decorrelate_context}()}) needed by the
#' shared \code{get_decorrelate_*} covariance helpers, so those helpers can be
#' reused unchanged by \code{\link{get_ssn_simulate_vecchia}()}.
#'
#' @param ssn.object A fitted-model-ready SSN object.
#' @param params_object A joint covariance parameter object (tailup, taildown,
#'   euclid, nugget, and optionally randcov).
#' @param additive The additive function value column name, or \code{NULL}.
#' @param anisotropy Whether Euclidean anisotropy is active.
#' @param partition_factor A one-sided partition factor formula, or
#'   \code{NULL}.
#'
#' @return A list with \code{ssn.object}, \code{coefficients$params_object},
#'   \code{additive}, \code{anisotropy}, \code{partition_factor},
#'   \code{diagtol}, \code{random_xlev}, and \code{partition_xlev}.
#'
#' @noRd
get_ssn_simulate_covariance_fit <- function(ssn.object, params_object, additive, anisotropy, partition_factor) {
  randcov_names <- names(params_object$randcov)
  random_xlev <- if (is.null(randcov_names)) NULL else get_randcov_xlev(randcov_names, ssn.object$obs)
  partition_xlev <- get_partition_xlev(partition_factor, ssn.object$obs)
  list(
    ssn.object = ssn.object,
    coefficients = list(params_object = params_object),
    additive = additive,
    anisotropy = anisotropy,
    partition_factor = partition_factor,
    diagtol = 0,
    random_xlev = random_xlev,
    partition_xlev = partition_xlev
  )
}

#' Simulate unconditional Gaussian draws via sequential covariance-neighbor conditioning
#'
#' A Vecchia-style sequential simulation engine: observations are ordered via
#' \code{\link{get_decorrelate_order}()}, and each point after the first is
#' drawn from its conditional distribution given (up to) \code{size} of the
#' most-correlated earlier-ordered points, reusing the same neighbor/
#' covariance helpers as \code{\link{ssn_decorrelate_data}()}'s big-data path.
#' \code{conditioning$method = "all"} conditions on every earlier point,
#' recovering the exact joint Gaussian distribution.
#'
#' @param covariance_fit A \code{covariance_fit}-shaped list (from
#'   \code{\link{get_ssn_simulate_covariance_fit}()} or
#'   \code{\link{get_decorrelate_context}()}) describing the joint covariance.
#' @param samples The number of simulated columns to draw.
#' @param conditioning A resolved \code{"vecchia"} conditioning list (from
#'   \code{\link{get_ssn_simulate_local}()}/\code{\link{get_local_vecchia_settings}()})
#'   with \code{method}, \code{size}, and \code{ordering}.
#'
#' @return An \code{n x samples} matrix of simulated values, in the original
#'   (unpermuted) row order.
#'
#' @noRd
get_ssn_simulate_vecchia <- function(covariance_fit, samples, conditioning) {
  observed <- covariance_fit$ssn.object$obs
  n <- NROW(observed)
  rows <- get_decorrelate_rows(observed)
  permutation <- get_decorrelate_order(rows, conditioning$ordering, observed)
  size <- if (identical(conditioning$method, "all")) n else conditioning$size

  Z <- matrix(rnorm(n * samples), n, samples)
  Y_ordered <- matrix(NA_real_, n, samples)

  for (position in seq_len(n)) {
    current <- permutation[[position]]
    variance <- get_decorrelate_marginal_variance(covariance_fit, observed[current, , drop = FALSE])

    if (position == 1L) {
      Y_ordered[position, ] <- sqrt(variance) * Z[position, ]
      next
    }

    pool_positions <- seq_len(position - 1L)
    candidates <- permutation[pool_positions]
    cross_covariance <- get_decorrelate_observed_cross_covariance(
      covariance_fit, observed[current, , drop = FALSE], observed[candidates, , drop = FALSE]
    )
    keep <- get_decorrelate_covariance_neighbors(cross_covariance, size)
    neighbor_positions <- pool_positions[keep]
    neighbor_index <- candidates[keep]

    covariance_neighbors <- get_decorrelate_observed_covariance(
      covariance_fit, observed[neighbor_index, , drop = FALSE]
    )
    chol_pool <- tryCatch(chol(covariance_neighbors), error = function(error) NULL)
    if (is.null(chol_pool)) {
      stop(
        "The neighbor covariance matrix for simulated record ", position,
        " (of the requested ordering) is not numerically positive definite; check covariance parameters and duplicate locations.",
        call. = FALSE
      )
    }

    c_vec <- cross_covariance[keep]
    w <- backsolve(chol_pool, forwardsolve(t(chol_pool), c_vec))
    cond_var <- max(variance - sum(w * c_vec), 0)
    pool_val <- Y_ordered[neighbor_positions, , drop = FALSE]
    cond_mean <- as.numeric(crossprod(w, pool_val))

    Y_ordered[position, ] <- cond_mean + sqrt(cond_var) * Z[position, ]
  }

  Y <- matrix(NA_real_, n, samples)
  Y[permutation, ] <- Y_ordered
  Y
}

#' Simulate one block of a low-rank base-and-block unconditional simulation
#'
#' Given a Gaussian process already simulated at a "base" set of locations
#' (\code{base_val}), draws values at \code{block_index}'s locations from
#' their conditional distribution given the base draws alone -- the
#' workhorse behind \code{\link{get_ssn_simulate_lowrank}()}'s big-data
#' approximation, mirroring spmodel's \code{get_conditional_new_from_base()}
#' (\code{R/get_conditional_new_from_base.R}). Takes every input as an
#' explicit argument (rather than closing over them) so it can be dispatched
#' via \code{parallel::parLapply()} without requiring \code{clusterExport()}.
#'
#' @param block_index Row indices (into \code{observed}) of this block's
#'   locations.
#' @param covariance_fit A \code{covariance_fit}-shaped list.
#' @param observed The full observed data (\code{covariance_fit$ssn.object$obs}).
#' @param base_data The base sample's own rows (\code{observed[index$base, ]}).
#' @param base_val An \code{n_base x samples} matrix of simulated base-sample
#'   values.
#' @param base_lowchol The lower triangular Cholesky factor of the base
#'   sample's covariance matrix.
#' @param samples The number of simulations (columns of \code{base_val}).
#'
#' @return An \code{length(block_index) x samples} matrix of simulated
#'   values for this block.
#'
#' @noRd
get_ssn_simulate_new_from_base <- function(block_index, covariance_fit, observed, base_data,
                                            base_val, base_lowchol, samples) {
  block_data <- observed[block_index, , drop = FALSE]
  n_block <- length(block_index)
  n_base <- NROW(base_data)

  cross_covariance <- get_decorrelate_observed_cross_covariance(covariance_fit, base_data, block_data)
  cross_covariance <- matrix(cross_covariance, nrow = n_base, ncol = n_block)
  block_covariance <- get_decorrelate_observed_covariance(covariance_fit, block_data)

  sqrt_siginv_cross <- forwardsolve(base_lowchol, cross_covariance)
  sqrt_siginv_base <- forwardsolve(base_lowchol, base_val)

  cond_cov <- block_covariance - crossprod(sqrt_siginv_cross, sqrt_siginv_cross)
  cond_cov <- as.matrix(Matrix::forceSymmetric(cond_cov))
  cond_lowchol <- chol_lower_with_pivot_fallback(
    cond_cov,
    "The block conditional covariance matrix for local low-rank simulation is not numerically positive semidefinite; check covariance parameters and duplicate locations."
  )

  z <- matrix(rnorm(n_block * samples), n_block, samples)
  cond_mu <- crossprod(sqrt_siginv_cross, sqrt_siginv_base)
  (cond_lowchol %*% z) + cond_mu
}

#' Simulate unconditional Gaussian draws via a low-rank base+block approximation
#'
#' Mirrors spmodel's low-rank branch of \code{sprnorm.exponential()}
#' (\code{R/sprnorm.R}): draws the base sample's values jointly (its own
#' Cholesky factor), then simulates every other location's block
#' conditionally on the base sample alone via
#' \code{\link{get_ssn_simulate_new_from_base}()}, treating distinct blocks
#' as conditionally independent given the base. Unlike
#' \code{\link{get_ssn_simulate_vecchia}()}, blocks do not depend on each
#' other, so \code{local$parallel} genuinely parallelizes this engine.
#'
#' @param covariance_fit A \code{covariance_fit}-shaped list.
#' @param samples The number of simulated columns to draw.
#' @param local A resolved \code{local} list (from
#'   \code{\link{get_ssn_simulate_local_lowrank}()}) with \code{method_base}
#'   and, when \code{method_base != "all"}, \code{index}.
#'
#' @return An \code{n x samples} matrix of simulated values, in the original
#'   row order.
#'
#' @noRd
get_ssn_simulate_lowrank <- function(covariance_fit, samples, local) {
  observed <- covariance_fit$ssn.object$obs
  n <- NROW(observed)

  if (identical(local$method_base, "all")) {
    base_index <- seq_len(n)
    new_index <- list()
  } else {
    base_index <- local$index$base
    new_index <- local$index$new
  }
  base_data <- observed[base_index, , drop = FALSE]
  n_base <- NROW(base_data)

  base_covariance <- get_decorrelate_observed_covariance(covariance_fit, base_data)
  base_lowchol <- chol_lower_with_pivot_fallback(
    base_covariance,
    "The base sample covariance matrix for local low-rank simulation is not numerically positive semidefinite; check covariance parameters and duplicate locations."
  )
  base_val <- base_lowchol %*% matrix(rnorm(n_base * samples), n_base, samples)

  Y <- matrix(NA_real_, n, samples)
  Y[base_index, ] <- as.matrix(base_val)

  if (length(new_index) == 0L) {
    return(Y)
  }

  if (local$parallel) {
    cl <- parallel::makeCluster(local$ncores)
    new_val <- parallel::parLapply(
      cl, new_index, get_ssn_simulate_new_from_base,
      covariance_fit, observed, base_data, base_val, base_lowchol, samples
    )
    parallel::stopCluster(cl)
  } else {
    new_val <- lapply(
      new_index, get_ssn_simulate_new_from_base,
      covariance_fit, observed, base_data, base_val, base_lowchol, samples
    )
  }

  for (i in seq_along(new_index)) {
    Y[new_index[[i]], ] <- as.matrix(new_val[[i]])
  }
  Y
}
