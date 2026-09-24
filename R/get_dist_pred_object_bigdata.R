get_distjunc_pred_matlist_bigdata_one_network <- function(x, obs_pid, pred_pid, ssn.object, newdata_name) {
  n_obs <- length(obs_pid)
  n_pred <- length(pred_pid)
  if (n_obs == 0 || n_pred == 0) {
    return(list(
      distjunca = matrix(0, nrow = n_obs, ncol = n_pred),
      distjuncb = matrix(0, nrow = n_pred, ncol = n_obs)
    ))
  }

  obs_pid_sorted <- as.character(sort(as.numeric(obs_pid)))
  pred_pid_sorted <- as.character(sort(as.numeric(pred_pid)))

  workspace.name.a <- paste("dist.net", x, ".a.bmat", sep = "")
  workspace.name.b <- paste("dist.net", x, ".b.bmat", sep = "")
  path.a <- file.path(ssn.object$path, "distance", newdata_name, workspace.name.a)
  path.b <- file.path(ssn.object$path, "distance", newdata_name, workspace.name.b)
  if (!file.exists(path.a) || !file.exists(path.b)) {
    stop("Unable to locate required distance matrix", call. = FALSE)
  }

  # .a: bmat stores pred rows x obs cols, the transpose of the dense .a
  fm_a <- fm.open(path.a)
  on.exit(close(fm_a), add = TRUE)
  pred_match_a <- match(pred_pid_sorted, rownames(fm_a))
  obs_match_a <- match(obs_pid_sorted, colnames(fm_a))
  missing_a <- c(pred_pid_sorted[is.na(pred_match_a)], obs_pid_sorted[is.na(obs_match_a)])
  if (length(missing_a) > 0) {
    stop(
      "Unable to locate stored distance information for the following pid(s): ",
      paste(missing_a, collapse = ", "), call. = FALSE
    )
  }
  distjunca <- t(fm_a[pred_match_a, obs_match_a])
  rownames(distjunca) <- obs_pid_sorted
  colnames(distjunca) <- pred_pid_sorted

  # .b: bmat is already pred rows x obs cols, same orientation as dense .b
  fm_b <- fm.open(path.b)
  on.exit(close(fm_b), add = TRUE)
  pred_match_b <- match(pred_pid_sorted, rownames(fm_b))
  obs_match_b <- match(obs_pid_sorted, colnames(fm_b))
  missing_b <- c(pred_pid_sorted[is.na(pred_match_b)], obs_pid_sorted[is.na(obs_match_b)])
  if (length(missing_b) > 0) {
    stop(
      "Unable to locate stored distance information for the following pid(s): ",
      paste(missing_b, collapse = ", "), call. = FALSE
    )
  }
  distjuncb <- fm_b[pred_match_b, obs_match_b]
  rownames(distjuncb) <- pred_pid_sorted
  colnames(distjuncb) <- obs_pid_sorted

  list(distjunca = distjunca, distjuncb = distjuncb)
}

get_distjunc_pred_matlist_bigdata <- function(ssn.object, newdata_name, order_list_pred) {
  network_index_obs <- as.numeric(as.character(order_list_pred$network_index))
  network_index_pred <- as.numeric(as.character(order_list_pred$network_index_pred))
  network_index_vals <- sort(unique(c(network_index_obs, network_index_pred)))
  network_pid_obs <- as.character(order_list_pred$pid)
  network_pid_pred <- as.character(order_list_pred$pid_pred)
  is_missing <- identical(newdata_name, ".missing")

  distjunc_pred_matlist <- lapply(network_index_vals, function(x) {
    obs_pid_x <- network_pid_obs[network_index_obs == x]
    pred_pid_x <- network_pid_pred[network_index_pred == x]

    if (length(obs_pid_x) == 0 || length(pred_pid_x) == 0) {
      return(list(
        distjunca = matrix(0, nrow = length(obs_pid_x), ncol = length(pred_pid_x)),
        distjuncb = matrix(0, nrow = length(pred_pid_x), ncol = length(obs_pid_x))
      ))
    }

    if (is_missing) {
      # The caller restores input order from PID-sorted network blocks.
      obs_pid_x <- as.character(sort(as.numeric(obs_pid_x)))
      pred_pid_x <- as.character(sort(as.numeric(pred_pid_x)))
      distjunca <- get_distjunc_matlist_bigdata_cross(
        rep(x, length(obs_pid_x)), obs_pid_x,
        rep(x, length(pred_pid_x)), pred_pid_x,
        ssn.object,
        ext = "obs"
      )
      distjuncb <- get_distjunc_matlist_bigdata_cross(
        rep(x, length(pred_pid_x)), pred_pid_x,
        rep(x, length(obs_pid_x)), obs_pid_x,
        ssn.object,
        ext = "obs"
      )
      list(distjunca = distjunca, distjuncb = distjuncb)
    } else {
      get_distjunc_pred_matlist_bigdata_one_network(x, obs_pid_x, pred_pid_x, ssn.object, newdata_name)
    }
  })

  distjunca <- lapply(distjunc_pred_matlist, function(x) x$distjunca)
  distjuncb <- lapply(distjunc_pred_matlist, function(x) x$distjuncb)

  list(distjunca = distjunca, distjuncb = distjuncb)
}

get_local_cov_means_bigdata <- function(object, newdata_name, initial_object, params_object) {
  netgeom_pred <- ssn_get_netgeom(object$ssn.object$preds[[newdata_name]], reformat = TRUE)
  network_index_pred <- as.numeric(as.character(netgeom_pred$NetworkID))
  network_vals <- sort(unique(network_index_pred))

  sum_acc <- NULL
  m_total <- 0

  for (x in network_vals) {
    pred_index_x <- which(network_index_pred == x)
    if (length(pred_index_x) == 0) {
      next
    }

    # one network's prediction rows at a time is the bounded block: the full
    # obs.pred cross block is never materialized across all networks at once
    chunk_object <- object
    chunk_object$ssn.object$preds[[newdata_name]] <- object$ssn.object$preds[[newdata_name]][pred_index_x, , drop = FALSE]

    dist_pred_object_x <- get_dist_pred_object(chunk_object, newdata_name, initial_object, backend = "bigdata")
    cov_vector_x <- get_cov_vector(
      params_object, dist_pred_object_x, object$ssn.object$obs,
      chunk_object$ssn.object$preds[[newdata_name]],
      object$partition_factor, object$anisotropy
    )
    col_sums_x <- colSums(as.matrix(cov_vector_x))
    sum_acc <- if (is.null(sum_acc)) col_sums_x else sum_acc + col_sums_x
    m_total <- m_total + length(pred_index_x)
  }

  sum_acc / m_total
}

get_local_cov_vector_means <- function(object, newdata_name, params_object) {
  tailup_none <- inherits(params_object$tailup, "tailup_none")
  taildown_none <- inherits(params_object$taildown, "taildown_none")
  backend <- select_pred_dist_backend(object$ssn.object, newdata_name, tailup_none, taildown_none)

  if (identical(backend, "bigdata")) {
    initial_object <- get_initial_object_from_coef(object)
    get_local_cov_means_bigdata(object, newdata_name, initial_object, params_object)
  } else {
    colMeans(covmatrix(object, newdata_name))
  }
}
