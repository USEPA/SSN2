get_conditional_vecchia_covariance <- function(object, newdata_name = NULL, spatial = FALSE) {
  if (spatial) object <- conditional_spatial_object(object)
  observed <- object$ssn.object$obs
  newdata <- if (is.null(newdata_name)) observed[FALSE, , drop = FALSE] else
    object$ssn.object$preds[[newdata_name]]
  if (!is.null(newdata_name) && is.null(newdata)) stop("Invalid prediction set name.", call. = FALSE)
  design <- get_conditional_random_design(object, newdata)
  partition <- if (is.null(object$partition_factor)) NULL else {
    f <- reformulate(labels(terms(object$partition_factor)), intercept = FALSE, env = asNamespace("SSN2"))
    randcov_design_key(rbind(model.frame(f, observed, na.action = na.pass),
                            model.frame(f, newdata, na.action = na.pass)))
  }
  backends <- list(observed = get_decorrelate_observed_backend(object))
  if (NROW(newdata)) {
    backends$cross <- get_block_pred_backend(object, newdata_name)
    backends$newdata <- get_block_pred_backend(object, newdata_name, square = TRUE)
  }
  rows <- get_decorrelate_rows(observed)
  coords <- sf::st_coordinates(observed)
  if (NROW(newdata)) {
    rows <- rbind(rows, get_decorrelate_rows(newdata))
    coords <- rbind(coords, sf::st_coordinates(newdata))
  }
  spatial_object <- conditional_spatial_object(object)
  spatial_object$partition_factor <- NULL
  list(object = object, spatial = spatial_object, observed = observed, newdata = newdata,
       newdata_name = newdata_name, n = NROW(observed), rows = rows, coords = coords,
       backends = backends, random = rbind(design$observed, design$newdata),
       randcov = design$variance, effect_info = design$info, partition = partition,
       initial = get_initial_object_from_coef(object))
}

conditional_vecchia_covariance <- function(context, rows, columns = rows) {
  n <- context$n
  total <- n + NROW(context$newdata)
  for (index in list(rows, columns)) {
    if (!is.numeric(index) || anyNA(index) || any(!is.finite(index)) ||
        any(index < 1 | index > total | index != floor(index))) {
      stop("Covariance indices must identify rows in the conditioning support.", call. = FALSE)
    }
  }
  covariance <- matrix(0, length(rows), length(columns))
  if (!length(rows) || !length(columns)) return(covariance)
  object <- context$spatial
  params <- object$coefficients$params_object
  data_object <- list(ssn.object = object$ssn.object, additive = object$additive,
                      anisotropy = object$anisotropy)
  for (pred_rows in c(FALSE, TRUE)) for (pred_cols in c(FALSE, TRUE)) {
    ir <- which((rows > n) == pred_rows)
    ic <- which((columns > n) == pred_cols)
    if (!length(ir) || !length(ic)) next
    r <- rows[ir] - if (pred_rows) n else 0L
    c <- columns[ic] - if (pred_cols) n else 0L
    if (pred_rows == pred_cols) {
      data <- if (pred_rows) context$newdata else context$observed
      ext <- if (!pred_rows || identical(context$newdata_name, ".missing")) "obs" else context$newdata_name
      backend <- if (pred_rows) context$backends$newdata else context$backends$observed
      distance <- get_dist_object_bigdata_cross(data[r, , drop = FALSE], data[c, , drop = FALSE],
                                               params, data_object, ext = ext, backend = backend)
      block <- get_cov_matrix_cross(params, distance, data_object = data_object)
    } else {
      obs <- if (pred_rows) c else r
      pred <- if (pred_rows) r else c
      sub <- object
      sub$ssn.object$obs <- context$observed[obs, , drop = FALSE]
      sub$ssn.object$preds[[context$newdata_name]] <- context$newdata[pred, , drop = FALSE]
      distance <- get_dist_pred_object(sub, context$newdata_name, context$initial,
                                       backend = context$backends$cross)
      block <- get_cov_vector(params, distance, sub$ssn.object$obs,
        sub$ssn.object$preds[[context$newdata_name]], anisotropy = object$anisotropy)
      if (!pred_rows) block <- t(block)
    }
    covariance[ir, ic] <- as.matrix(block)
  }
  if (length(context$randcov)) {
    covariance <- covariance + as.matrix(context$random[rows, , drop = FALSE] %*%
      Matrix::Diagonal(x = context$randcov) %*% t(context$random[columns, , drop = FALSE]))
  }
  dependent <- sum(params$tailup[["de"]], params$taildown[["de"]], params$euclid[["de"]])
  nugget <- get_spatial_nugget_var(params, object$diagtol) - dependent
  covariance <- covariance + nugget * outer(rows, columns, "==")
  if (!is.null(context$partition)) {
    covariance <- covariance * outer(context$partition[rows], context$partition[columns], "==")
  }
  dimnames(covariance) <- list(as.character(rows), as.character(columns))
  covariance
}

conditional_vecchia_scores <- function(context, i, method, columns = seq_len(context$n + NROW(context$newdata))) {
  if (method != "covariance") stop('local$method must be "all" or "covariance".', call. = FALSE)
  score <- stats::setNames(rep(Inf, length(columns)), as.character(columns))
  eligible <- columns != i
  if (!is.null(context$partition)) eligible <- eligible & context$partition[columns] == context$partition[i]
  score[eligible] <- -abs(as.numeric(conditional_vecchia_covariance(context, i, columns[eligible])))
  score
}

conditional_vecchia_neighbors <- function(context, i, pool, size, method, scores = NULL) {
  if (size == 0L && method != "all") return(integer())
  if (method == "all") return(pool[pool != i])
  if (is.null(scores)) scores <- conditional_vecchia_scores(context, i, method, pool)
  score <- scores[as.character(pool)]
  eligible <- is.finite(score)
  pool <- pool[eligible]
  score <- score[eligible]
  if (length(pool) <= size) return(pool)
  cutoff <- sort(score, partial = size)[size]
  candidates <- which(score <= cutoff)
  tie_index <- if (is.null(context$tie_index)) pool[candidates] else context$tie_index[pool[candidates]]
  pool[candidates[order(score[candidates], tie_index)][seq_len(size)]]
}
