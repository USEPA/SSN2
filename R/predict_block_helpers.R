#' Select the distance backend for block/chunked prediction covariance
#'
#' Returns \code{"none"} when there is no stream covariance to compute.
#' Otherwise defers entirely to
#' \code{\link{select_square_dist_backend}()}/\code{\link{select_pred_dist_backend}()}.
#' Block prediction, conditional simulation, and decorrelation request many
#' small submatrices. Prefer \code{.bmat} for indexed reads; the dense
#' \code{.RData} fallback loads each whole network matrix before subsetting
#' and does not cache it between chunks.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set, or \code{".missing"}.
#' @param square Whether the block being computed is prediction-to-prediction
#'   (\code{TRUE}) rather than prediction-to-observed (\code{FALSE}).
#' @param prefer Which backend to prefer when both exist; see
#'   \code{\link{select_square_dist_backend}()}. Defaults to \code{"bigdata"}.
#'
#' @return \code{"none"}, or the backend selected by the relevant
#'   \code{select_*_dist_backend()} helper (\code{"dense"} or
#'   \code{"bigdata"}).
#'
#' @noRd
get_block_pred_backend <- function(object, newdata_name, square = FALSE, prefer = "bigdata") {
  initial_object <- get_initial_object_from_coef(object)
  tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
  taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")
  if (tailup_none && taildown_none) return("none")

  if (square) {
    ext <- if (identical(newdata_name, ".missing")) "obs" else newdata_name
    select_square_dist_backend(object$ssn.object, ext, tailup_none, taildown_none, prefer = prefer)
  } else {
    select_pred_dist_backend(object$ssn.object, newdata_name, tailup_none, taildown_none, prefer = prefer)
  }
}

#' Split a row index into fixed-size chunks
#'
#' @param index An integer index to split.
#' @param chunk_size The maximum chunk size.
#'
#' @return A list of index chunks, each of length at most \code{chunk_size}.
#'
#' @noRd
get_prediction_chunks <- function(index, chunk_size) {
  split(index, ceiling(seq_along(index) / chunk_size))
}

#' Apply a function over prediction chunks, optionally in parallel
#'
#' Falls back to \code{lapply()} when \code{parallel} is \code{FALSE} or
#' there is only one chunk. Otherwise starts a cluster and loads SSN2 on each
#' worker (via \code{pkgload::load_all()} if running from a development
#' install, otherwise \code{loadNamespace()}) before applying \code{fun}.
#'
#' @param chunks A list of index chunks, from
#'   \code{\link{get_prediction_chunks}()}.
#' @param fun A function to apply to each chunk.
#' @param parallel Whether to use parallel processing.
#' @param ncores The number of cluster workers to use when \code{parallel} is
#'   \code{TRUE}.
#'
#' @return A list of \code{fun}'s results, one per chunk.
#'
#' @noRd
get_block_chunk_apply <- function(chunks, fun, parallel = FALSE, ncores = NULL) {
  if (!isTRUE(parallel) || length(chunks) < 2L) {
    return(lapply(chunks, fun))
  }
  nworkers <- min(as.integer(ncores), length(chunks))
  cl <- parallel::makeCluster(nworkers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  dev_path <- NULL
  if (requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("SSN2")) {
    source_file <- utils::getSrcFilename(utils::getSrcref(get_block_chunk_apply), full.names = TRUE)
    candidate <- normalizePath(file.path(dirname(source_file), ".."), mustWork = TRUE)
    if (file.exists(file.path(candidate, "DESCRIPTION"))) dev_path <- candidate
  }
  parallel::clusterCall(cl, function(dev_path) {
    if (is.null(dev_path)) {
      loadNamespace("SSN2")
    } else {
      pkgload::load_all(dev_path, quiet = TRUE)
    }
    NULL
  }, dev_path)
  parallel::parLapply(cl, chunks, fun)
}

#' Build the covariance block between a chunk of prediction rows and every observed row
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set, or \code{".missing"}.
#' @param rows The prediction rows (within \code{newdata_name}) to include.
#' @param backend The distance backend to read with (\code{"dense"} or
#'   \code{"bigdata"}). When \code{NULL} (the default), resolved via
#'   \code{\link{get_block_pred_backend}()}, which prefers \code{"bigdata"}
#'   with dense fallback. Pass an already-resolved backend to keep chunk reads
#'   consistent with the caller's selection.
#' @param randcov_context A precomputed random-effect training-side context
#'   from \code{\link{get_randcov_context}()}, or \code{NULL} (re-derived
#'   internally). Pass one in when calling this repeatedly for the same fit
#'   (e.g. once per row-chunk) to avoid rebuilding it every call.
#' @param partition_context A precomputed partition-factor training-side
#'   context from \code{\link{get_partition_context}()}, or \code{NULL}
#'   (re-derived internally). Pass one in when calling this repeatedly for the
#'   same fit (e.g. once per row-chunk) to avoid rebuilding it every call.
#'
#' @return A \code{length(rows) x object$n} covariance matrix.
#'
#' @noRd
get_block_obs_covariance <- function(object, newdata_name, rows, backend = NULL, randcov_context = NULL, partition_context = NULL) {
  chunk_object <- object
  chunk_object$ssn.object$preds[[newdata_name]] <- object$ssn.object$preds[[newdata_name]][rows, , drop = FALSE]
  params_object <- object$coefficients$params_object
  initial_object <- get_initial_object_from_coef(object)
  if (is.null(backend)) backend <- get_block_pred_backend(object, newdata_name)
  dist_object <- get_dist_pred_object(
    chunk_object, newdata_name, initial_object,
    backend = backend
  )
  covariance <- as.matrix(get_cov_vector(
    params_object, dist_object, object$ssn.object$obs,
    chunk_object$ssn.object$preds[[newdata_name]], object$partition_factor,
    object$anisotropy, randcov_context = randcov_context, partition_context = partition_context
  ))
  expected_dim <- c(NROW(chunk_object$ssn.object$preds[[newdata_name]]), object$n)
  if (!identical(dim(covariance), expected_dim) && length(covariance) == 1L) {
    covariance <- matrix(covariance, nrow = expected_dim[[1]], ncol = expected_dim[[2]])
  }
  covariance
}

#' Build the covariance block between two chunks of prediction rows
#'
#' Assembles stream/Euclidean covariance plus, if active, random-effect and
#' partition-factor contributions, and adds the nugget variance on any
#' diagonal entries (matching rows/columns representing the same
#' \code{NetworkID}/\code{pid}).
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set, or \code{".missing"}.
#' @param rows,columns The prediction rows (within \code{newdata_name}) for
#'   each side of the block.
#'
#' @return A \code{length(rows) x length(columns)} covariance matrix.
#'
#' @noRd
get_block_pred_covariance <- function(object, newdata_name, rows, columns) {
  newdata <- object$ssn.object$preds[[newdata_name]]
  d1 <- newdata[rows, , drop = FALSE]
  d2 <- newdata[columns, , drop = FALSE]
  params_object <- object$coefficients$params_object
  data_object <- list(
    ssn.object = object$ssn.object,
    additive = object$additive,
    anisotropy = object$anisotropy
  )
  ext <- if (identical(newdata_name, ".missing")) "obs" else newdata_name
  dist_object <- get_dist_object_bigdata_cross(
    d1, d2, params_object, data_object, ext = ext,
    backend = get_block_pred_backend(object, newdata_name, square = TRUE)
  )
  covariance <- get_cov_matrix_cross(
    params_object, dist_object, data_object = data_object
  )
  if (!identical(dim(covariance), c(NROW(d1), NROW(d2))) && length(covariance) == 1L) {
    covariance <- matrix(covariance, nrow = NROW(d1), ncol = NROW(d2))
  }
  if (!is.null(params_object$randcov)) {
    randcov_names <- get_randcov_names(object$random)
    xlev_list <- extend_randcov_xlev(object$random_xlev, newdata, randcov_names)
    covariance <- covariance + randcov_vector(params_object$randcov, d2, d1, xlev_list = xlev_list)
  }
  partition_xlev <- extend_partition_xlev(
    object$partition_factor, object$partition_xlev, newdata
  )
  partition <- partition_vector(object$partition_factor, d2, d1, xlev = partition_xlev)
  if (!is.null(partition)) {
    covariance <- covariance * partition
  }

  d1_netgeom <- ssn_get_netgeom(d1, reformat = TRUE)
  d2_netgeom <- ssn_get_netgeom(d2, reformat = TRUE)
  d1_key <- paste(d1_netgeom$NetworkID, d1_netgeom$pid, sep = "\r")
  d2_key <- paste(d2_netgeom$NetworkID, d2_netgeom$pid, sep = "\r")
  diagonal <- match(d1_key, d2_key)
  has_diagonal <- which(!is.na(diagonal))
  if (length(has_diagonal)) {
    dependent_variance <- sum(
      params_object$tailup[["de"]], params_object$taildown[["de"]],
      params_object$euclid[["de"]]
    )
    nugget_variance <- get_spatial_nugget_var(params_object, object$diagtol) - dependent_variance
    multiplier <- if (is.null(partition)) 1 else partition[cbind(has_diagonal, diagonal[has_diagonal])]
    covariance[cbind(has_diagonal, diagonal[has_diagonal])] <-
      covariance[cbind(has_diagonal, diagonal[has_diagonal])] + nugget_variance * multiplier
  }
  as.matrix(covariance)
}

#' Compute block-prediction averages of observed- and within-block covariance
#'
#' Computes \code{c0} (the average covariance between every observed
#' location and the block's rows) and, if requested, \code{s0} (the average
#' within-block covariance, approximated using \code{nodes}, a representative
#' subset of block rows). Both are computed in chunks (optionally in
#' parallel) to bound peak memory.
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param newdata_name The name of the prediction set, or \code{".missing"}.
#' @param c0_rows The block's rows used to average \code{c0} over.
#' @param s0_rows The block's rows used to average \code{s0} over.
#' @param nodes The representative node rows used to approximate \code{s0}'s
#'   within-block covariance.
#' @param chunk_size The maximum number of rows processed per covariance
#'   chunk.
#' @param compute_s0 Whether to also compute \code{s0}.
#' @param parallel Whether to use parallel processing.
#' @param ncores The number of cluster workers to use when \code{parallel} is
#'   \code{TRUE}.
#'
#' @return A list with \code{c0} (a length-\code{object$n} vector) and
#'   \code{s0} (a scalar, or \code{NULL} if \code{compute_s0} is
#'   \code{FALSE}).
#'
#' @noRd
get_block_quantities <- function(object, newdata_name, c0_rows, s0_rows, nodes,
                                 chunk_size, compute_s0 = TRUE,
                                 parallel = FALSE, ncores = NULL) {
  n_obs <- object$n
  c0_chunks <- get_prediction_chunks(c0_rows, chunk_size)
  c0_parts <- get_block_chunk_apply(c0_chunks, function(rows) {
    colSums(get_block_obs_covariance(object, newdata_name, rows))
  }, parallel = parallel, ncores = ncores)
  c0_sum <- Reduce(`+`, c0_parts, init = numeric(n_obs))
  c0 <- c0_sum / length(c0_rows)
  if (!compute_s0) return(list(c0 = c0, s0 = NULL))

  row_chunks <- get_prediction_chunks(s0_rows, chunk_size)
  node_chunks <- get_prediction_chunks(nodes, chunk_size)
  node_sums <- get_block_chunk_apply(node_chunks, function(node_chunk) {
    node_sum <- 0
    for (row_chunk in row_chunks) {
      node_sum <- node_sum + sum(
        get_block_pred_covariance(object, newdata_name, row_chunk, node_chunk)
      )
    }
    node_sum
  }, parallel = parallel, ncores = ncores)
  s0_sum <- sum(unlist(node_sums, use.names = FALSE))
  list(c0 = c0, s0 = s0_sum / (length(s0_rows) * length(nodes)))
}

#' Sum the marginal (spatial + nugget + random-effect) variance over block rows
#'
#' Used for the within-block variance's diagonal contribution in block
#' prediction (the average pairwise within-block covariance excludes each
#' row's covariance with itself, whose marginal variance is added back
#' separately).
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param newdata_name The name of the prediction set, or \code{".missing"}.
#' @param rows The block's rows (within \code{newdata_name}) to sum over.
#' @param randcov_context A precomputed random-effect training-side context
#'   from \code{\link{get_randcov_context}()}, or \code{NULL} (re-derived
#'   per row).
#'
#' @return The sum of marginal variances over \code{rows}.
#'
#' @noRd
get_block_diagonal_sum <- function(object, newdata_name, rows, randcov_context = NULL) {
  newdata <- object$ssn.object$preds[[newdata_name]][rows, , drop = FALSE]
  spatial_nugget <- get_spatial_nugget_var(
    object$coefficients$params_object, object$diagtol
  )
  random <- vapply(seq_len(NROW(newdata)), function(i) {
    randcov_newvar(object$coefficients$params_object$randcov, newdata[i, , drop = FALSE], context = randcov_context)
  }, numeric(1))
  sum(spatial_nugget + random)
}

#' Select representative node rows for block-prediction approximation
#'
#' Orders rows deterministically by \code{pid}, then \code{NetworkID},
#' then input row position and takes an evenly spread subset across that order
#' when \code{ordering = "pid"}. Other orderings select the first nodes in
#' the requested order; \code{"random"} samples nodes uniformly.
#'
#' @param newdata The block's prediction data frame.
#' @param size_new The number of node rows to select.
#' @param ordering The node-selection order; see \code{\link{ssn_decorrelate}()}.
#'
#' @return An integer index (length \code{min(size_new, NROW(newdata))})
#'   selecting node rows.
#'
#' @noRd
get_block_nodes <- function(newdata, size_new, ordering) {
  ordering <- get_decorrelate_ordering(ordering)
  n <- NROW(newdata)
  size_new <- min(as.integer(size_new), n)
  if (size_new == n) return(seq_len(n))
  if (identical(ordering, "none")) return(seq_len(size_new))
  if (identical(ordering, "random")) return(sample.int(n, size_new))

  ordered <- get_decorrelate_order(get_decorrelate_rows(newdata), ordering, newdata)
  if (!identical(ordering, "pid")) return(ordered[seq_len(size_new)])
  ordered[floor((seq_len(size_new) - 0.5) * n / size_new) + 1L]
}
