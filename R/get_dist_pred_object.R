#' @param prefer Which backend to prefer when both \code{.RData} and
#'   \code{.bmat} distance matrices exist: \code{"dense"} (the default) reads
#'   \code{.RData}, \code{"bigdata"} reads \code{.bmat}. Either backend
#'   returns numerically identical distances; the two differ only in whether
#'   an entire network's matrix is deserialized into memory and then
#'   subsetted (\code{"dense"}) or the requested submatrix is read directly
#'   off disk via a memory-mapped \code{filematrix} (\code{"bigdata"}). Local
#'   big-data machinery (model-fitting \code{local}, block prediction, and
#'   the local/Vecchia/low-rank decorrelation, simulation, and conditional
#'   simulation neighbor-pool code) typically only ever needs small
#'   submatrices at a time and passes \code{"bigdata"}, so it is not left
#'   unreachable whenever \code{.RData} also happens to exist. Exact/default
#'   paths, which always need an entire network's matrix anyway, keep the
#'   \code{"dense"} default.
#' @noRd
select_pred_dist_backend <- function(ssn.object, newdata_name, tailup_none, taildown_none, prefer = "dense") {
  ext <- if (identical(newdata_name, ".missing")) "obs" else newdata_name
  rdata_pattern <- if (identical(newdata_name, ".missing")) "^dist\\.net[0-9]+\\.RData$" else "^dist\\.net[0-9]+\\.a\\.RData$"
  bmat_pattern <- if (identical(newdata_name, ".missing")) "^dist\\.net[0-9]+\\.bmat$" else "^dist\\.net[0-9]+\\.a\\.bmat$"
  select_dist_backend_from_patterns(
    ssn.object, ext, rdata_pattern, bmat_pattern, tailup_none, taildown_none,
    label = newdata_name, prefer = prefer
  )
}

select_square_dist_backend <- function(ssn.object, ext, tailup_none, taildown_none, prefer = "dense") {
  select_dist_backend_from_patterns(
    ssn.object, ext, "^dist\\.net[0-9]+\\.RData$", "^dist\\.net[0-9]+\\.bmat$", tailup_none, taildown_none,
    label = ext, prefer = prefer
  )
}

select_dist_backend_from_patterns <- function(ssn.object, ext, rdata_pattern, bmat_pattern,
                                              tailup_none, taildown_none, label, prefer = "dense") {
  if (tailup_none && taildown_none) {
    return("none")
  }

  dist_dir <- file.path(ssn.object$path, "distance", ext)
  rdata_exists <- length(list.files(dist_dir, pattern = rdata_pattern)) > 0
  bmat_exists <- length(list.files(dist_dir, pattern = bmat_pattern)) > 0

  if (identical(prefer, "bigdata")) {
    if (bmat_exists) return("bigdata")
    if (rdata_exists) return("dense")
  } else {
    if (rdata_exists) return("dense")
    if (bmat_exists) return("bigdata")
  }
  stop(
    "Unable to locate distance matrices for \"", label, "\". Checked for dense ",
    "(.RData) and filematrix (.bmat) distance files in \"", dist_dir, "\" and ",
    "found neither. Run ssn_create_distmat() or ssn_create_bigdist() for this ",
    "dataset first.",
    call. = FALSE
  )
}

get_initial_object_from_coef <- function(object) {
  get_initial_object(
    tailup_type = remove_covtype(class(coef(object, type = "tailup"))),
    taildown_type = remove_covtype(class(coef(object, type = "taildown"))),
    euclid_type = remove_covtype(class(coef(object, type = "euclid"))),
    nugget_type = remove_covtype(class(coef(object, type = "nugget"))),
    tailup_initial = NULL,
    taildown_initial = NULL,
    euclid_initial = NULL,
    nugget_initial = NULL
  )
}

get_dist_pred_object <- function(object, newdata_name, initial_object, backend = NULL) {
  # get netgeom
  netgeom <- ssn_get_netgeom(object$ssn.object$obs, reformat = TRUE)

  # get network index
  network_index <- netgeom$NetworkID

  # get pid
  pid <- netgeom$pid

  # distance order
  dist_order <- order(network_index, pid)

  # inverse of distance order
  inv_dist_order <- order(dist_order)

  # get netgeom
  netgeom_pred <- ssn_get_netgeom(object$ssn.object$preds[[newdata_name]], reformat = TRUE)

  # get network pred index
  network_index_pred <- netgeom_pred$NetworkID

  # get pid
  pid_pred <- netgeom_pred$pid

  # distance order
  dist_order_pred <- order(network_index_pred, pid_pred)

  # inverse of distance order
  inv_dist_order_pred <- order(dist_order_pred)

  # create "order" list for predictions
  order_list_pred <- list(
    network_index = network_index,
    pid = pid, dist_order = dist_order,
    inv_dist_order = inv_dist_order,
    network_index_pred = network_index_pred,
    pid_pred = pid_pred, dist_order_pred = dist_order_pred,
    inv_dist_order_pred = inv_dist_order_pred
  )

  if (is.null(backend)) {
    tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
    taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")
    backend <- select_pred_dist_backend(object$ssn.object, newdata_name, tailup_none, taildown_none)
  }

  # get list of prediction distance matrices in order of the original data
  dist_pred_matlist <- get_dist_pred_matlist(
    object$ssn.object, newdata_name, initial_object, object$additive,
    order_list_pred,
    backend = backend
  )

  # see whether euclid is none to avoid unnecessary computations
  euclid_none <- inherits(initial_object$euclid_initial, "euclid_none")
  if (euclid_none) {
    # not needed if no euclid covariance
    dist_pred_matlist <- c(dist_pred_matlist, list(euclid_mat = NULL))
    # dist_pred_matlist$euclid_matix <- NULL does not return anything
  } else {
    # find coordinates in pid (data) order
    obs_coords <- sf::st_coordinates(object$ssn.object$obs)
    pred_coords <- sf::st_coordinates(object$ssn.object$preds[[newdata_name]])

    # check for anisotropy and store accordingly
    if (object$anisotropy) {
      # store as vectors with drop
      dist_pred_matlist$.obs_xcoord <- obs_coords[, 1, drop = TRUE]
      dist_pred_matlist$.obs_ycoord <- obs_coords[, 2, drop = TRUE]
      dist_pred_matlist$.pred_xcoord <- pred_coords[, 1, drop = TRUE]
      dist_pred_matlist$.pred_ycoord <- pred_coords[, 2, drop = TRUE]
    } else {
      dist_vector_x <- outer(X = obs_coords[, 1], Y = pred_coords[, 1], FUN = function(X, Y) (X - Y)^2)
      dist_vector_y <- outer(X = obs_coords[, 2], Y = pred_coords[, 2], FUN = function(X, Y) (X - Y)^2)
      dist_vector <- sqrt(dist_vector_x + dist_vector_y)
      dist_pred_matlist$euclid_pred_mat <- Matrix::Matrix(dist_vector, sparse = TRUE)
    }
  }

  # transpose the matrices so dimensions are usable with predict()
  dist_pred_matlist <- lapply(dist_pred_matlist, function(x) if (is.null(x)) NULL else t(x))

  # return relevant prediction data object and order list
  dist_pred_object <- c(dist_pred_matlist, order_list_pred)

  # return prediction distance object
  dist_pred_object
}

# vectorized version of get_dist_pred_object
get_dist_pred_matlist <- function(ssn.object, newdata_name, initial_object, additive,
                                  order_list_pred, backend = "dense") {
  # store network indices and orders
  network_index <- order_list_pred$network_index
  dist_order <- order_list_pred$dist_order
  inv_dist_order <- order_list_pred$inv_dist_order
  inv_dist_order_pred <- order_list_pred$inv_dist_order_pred

  # see whether tailup and taildown are none to avoid unnecessary computations
  tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
  taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")

  # return all NULL if they are both none (no stream distance needed)
  if (tailup_none && taildown_none) {
    dist_pred_matlist <- list(
      distjunc_pred_matlist = NULL,
      mask_pred_matlist = NULL,
      a_pred_matlist = NULL,
      b_pred_matlist = NULL,
      hydro_pred_matlist = NULL,
      w_pred_matlist = NULL
    )
  } else {
    # otherwise

    distjunc_pred_matlist <- if (identical(backend, "bigdata")) {
      get_distjunc_pred_matlist_bigdata(ssn.object, newdata_name, order_list_pred)
    } else {
      get_distjunc_pred_matlist(ssn.object, newdata_name, order_list_pred)
    }

    # get other matrices as a list
    dist_pred_matlist <- list(
      distjunc_pred_matlist = distjunc_pred_matlist,
      mask_pred_matlist = get_mask_pred_matlist(distjunc_pred_matlist),
      a_pred_matlist = get_a_pred_matlist(distjunc_pred_matlist),
      b_pred_matlist = get_b_pred_matlist(distjunc_pred_matlist),
      hydro_pred_matlist = get_hydro_pred_matlist(distjunc_pred_matlist)
    )

    # if only taildown covariance, do not need additive matrix
    if (tailup_none) {
      # create distance pred matrix list (0's implied by bdiag get zeroed out
      # in covariance by mask matrix)
      dist_pred_matlist <- list(
        distjunca_pred_mat = Matrix::bdiag(dist_pred_matlist$distjunc_pred_matlist$distjunca),
        distjuncb_pred_mat = Matrix::bdiag(dist_pred_matlist$distjunc_pred_matlist$distjuncb),
        mask_pred_mat = Matrix::bdiag(dist_pred_matlist$mask_pred_matlist),
        a_pred_mat = Matrix::bdiag(dist_pred_matlist$a_pred_matlist),
        b_pred_mat = Matrix::bdiag(dist_pred_matlist$b_pred_matlist),
        hydro_pred_mat = Matrix::bdiag(dist_pred_matlist$hydro_pred_matlist)
      )

      dist_pred_matlist <- mapply(
        reorder_dist_pred_field, names(dist_pred_matlist), dist_pred_matlist,
        MoreArgs = list(inv_dist_order = inv_dist_order, inv_dist_order_pred = inv_dist_order_pred),
        SIMPLIFY = FALSE
      )

      # store additive matrix as NULL
      dist_pred_matlist <- c(dist_pred_matlist, list(w_pred_mat = NULL))
    } else {
      # compute additive matrix
      dist_pred_matlist$w_pred_matlist <- get_w_pred_matlist(
        ssn.object,
        newdata_name,
        order_list_pred,
        additive,
        dist_pred_matlist$b_pred_matlist,
        dist_pred_matlist$mask_pred_matlist
      )

      # create distance pred matrix list (0's implied by bdiag get zeroed out
      # in covariance by mask matrix)
      dist_pred_matlist <- list(
        distjunca_pred_mat = Matrix::bdiag(dist_pred_matlist$distjunc_pred_matlist$distjunca),
        distjuncb_pred_mat = Matrix::bdiag(dist_pred_matlist$distjunc_pred_matlist$distjuncb),
        mask_pred_mat = Matrix::bdiag(dist_pred_matlist$mask_pred_matlist),
        a_pred_mat = Matrix::bdiag(dist_pred_matlist$a_pred_matlist),
        b_pred_mat = Matrix::bdiag(dist_pred_matlist$b_pred_matlist),
        hydro_pred_mat = Matrix::bdiag(dist_pred_matlist$hydro_pred_matlist),
        w_pred_mat = Matrix::bdiag(dist_pred_matlist$w_pred_matlist)
      )

      dist_pred_matlist <- mapply(
        reorder_dist_pred_field, names(dist_pred_matlist), dist_pred_matlist,
        MoreArgs = list(inv_dist_order = inv_dist_order, inv_dist_order_pred = inv_dist_order_pred),
        SIMPLIFY = FALSE
      )
    }
  }
  # return distance prediction object
  dist_pred_matlist
}

reorder_dist_pred_field <- function(name, x, inv_dist_order, inv_dist_order_pred) {
  if (is.null(x)) {
    return(NULL)
  }
  if (identical(name, "distjuncb_pred_mat")) {
    x[inv_dist_order_pred, inv_dist_order, drop = FALSE]
  } else {
    x[inv_dist_order, inv_dist_order_pred, drop = FALSE]
  }
}

get_distjunc_pred_matlist <- function(ssn.object, newdata_name, order_list_pred) {
  # check and make sure there is missing data to predict
  if (newdata_name %in% names(ssn.object$preds) && NROW(ssn.object$preds[[newdata_name]]) == 0) {
    stop("No missing data to predict", call. = FALSE)
  }

  # get network index values and their unique entries
  network_index_obs <- as.numeric(as.character(order_list_pred$network_index))
  network_index_pred <- as.numeric(as.character(order_list_pred$network_index_pred))
  network_index_vals <- sort(unique(c(network_index_obs, network_index_pred)))
  # network_index_integer <- seq_along(network_index_vals)

  # get network pid
  network_pid_obs <- as.character(order_list_pred$pid)
  network_pid_pred <- as.character(order_list_pred$pid_pred)

  read_matrix <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    # Full distance matrix loaded before chunk subsetting; consider per-call caching for block kriging.
    unserialize(con)
  }

  # find distance junction prediction matrices (as a list) separately for each
  # network index
  distjunc_pred_matlist <- lapply(network_index_vals, function(x) {
    # find observations for each network index
    ind_obs <- which(network_index_obs == x)
    # find the number of observations having that index
    n_obs <- length(ind_obs)
    # find observations for each prediction network index
    ind_pred <- which(network_index_pred == x)
    # find number of predictions having that index
    n_pred <- length(ind_pred)

    # loop through as long as there are at least some observations for both
    if (n_obs != 0 && n_pred != 0) {
      # operate differently if observations are raw prediction data or induced
      # by NA values in the response
      if (newdata_name == ".missing") {
        # on the disk, distance matrices are stored by network
        workspace_name <- paste("dist.net", x, ".RData", sep = "")
        # path to the distance matrices on disk
        path <- file.path(ssn.object$path, "distance", "obs", workspace_name)
        # check to see if the file exists on the disk
        if (!file.exists(path)) {
          stop("Unable to locate required distance matrix", call. = FALSE)
        }
        distmat <- read_matrix(path)
        obs_match <- match(network_pid_obs[network_index_obs == x], rownames(distmat))
        pred_match <- match(network_pid_pred[network_index_pred == x], colnames(distmat))
        if (anyNA(obs_match) || anyNA(pred_match)) {
          stop("Unable to locate stored distance information for requested observation or prediction pid.", call. = FALSE)
        }
        distmata <- distmat[obs_match, pred_match, drop = FALSE]
        distmatb <- distmat[pred_match, obs_match, drop = FALSE]
      } else {
        # on the disk, distance matrices are stored by network
        workspace.name.a <- paste("dist.net", x,
          ".a.RData",
          sep = ""
        )
        workspace.name.b <- paste("dist.net", x,
          ".b.RData",
          sep = ""
        )
        # path to the distance matrices on disk
        path.a <- file.path(
          ssn.object$path,
          "distance", newdata_name, workspace.name.a
        )
        # check to see if the file exists on the disk
        if (!file.exists(path.a)) {
          stop("Unable to locate required distance matrix", call. = FALSE)
        }
        path.b <- file.path(
          ssn.object$path,
          "distance", newdata_name, workspace.name.b
        )
        # check to see if the file exists on the disk
        if (!file.exists(path.b)) {
          stop("Unable to locate required distance matrix", call. = FALSE)
        }
        distmata_all <- read_matrix(path.a)
        distmatb_all <- read_matrix(path.b)
        obs_match <- match(network_pid_obs[network_index_obs == x], rownames(distmata_all))
        pred_match <- match(network_pid_pred[network_index_pred == x], colnames(distmata_all))
        if (anyNA(obs_match) || anyNA(pred_match)) {
          stop("Unable to locate stored distance information for requested observation or prediction pid.", call. = FALSE)
        }
        distmata <- distmata_all[obs_match, pred_match, drop = FALSE]
        distmatb <- distmatb_all[pred_match, obs_match, drop = FALSE]
      }

      # find pid order
      pid_order_obs <- order(as.numeric(rownames(distmata)))
      pid_order_pred <- order(as.numeric(rownames(distmatb)))
      # return distance junction matrices
      distjunca <- distmata[pid_order_obs, pid_order_pred, drop = FALSE]
      distjuncb <- distmatb[pid_order_pred, pid_order_obs, drop = FALSE]
      # distjunca <- distmata # assumes they are ordered
      # distjuncb <- distmatb # assumes they are ordered
    } else {
      distjunca <- Matrix::Matrix(0, nrow = n_obs, ncol = n_pred)
      distjuncb <- t(distjunca)
    }
    list(distjunca = distjunca, distjuncb = distjuncb)
  })

  distjunca <- lapply(distjunc_pred_matlist, function(x) x$distjunca)
  distjuncb <- lapply(distjunc_pred_matlist, function(x) x$distjuncb)

  distjunc_pred_matlist <- list(distjunca = distjunca, distjuncb = distjuncb)
}

get_mask_pred_matlist <- function(distjunc_pred_matlist) {
  mask_pred_list <- lapply(distjunc_pred_matlist$distjunca, function(x) {
    Matrix::Matrix(1, nrow = dim(x)[1], ncol = dim(x)[2], sparse = TRUE)
  })
}

get_a_pred_matlist <- function(distjunc_pred_matlist) {
  a_matrix_list <- mapply(
    a = distjunc_pred_matlist$distjunca,
    b = distjunc_pred_matlist$distjuncb,
    function(a, b) {
      Matrix::Matrix(pmax(as.matrix(a), as.matrix(t(b))), sparse = TRUE)
    },
    SIMPLIFY = FALSE
  )
}

get_b_pred_matlist <- function(distjunc_pred_matlist) {
  a_matrix_list <- mapply(
    a = distjunc_pred_matlist$distjunca,
    b = distjunc_pred_matlist$distjuncb,
    function(a, b) {
      Matrix::Matrix(pmin(as.matrix(a), as.matrix(t(b))), sparse = TRUE)
    },
    SIMPLIFY = FALSE
  )
}

get_hydro_pred_matlist <- function(distjunc_pred_matlist) {
  a_matrix_list <- mapply(
    a = distjunc_pred_matlist$distjunca,
    b = distjunc_pred_matlist$distjuncb,
    function(a, b) {
      a + t(b)
    },
    SIMPLIFY = FALSE
  )
}

get_w_pred_matlist <- function(ssn.object, newdata_name, order_list_pred, additive, b_pred_matlist, mask_pred_matlist) {
  # make list
  network_index_obs <- as.numeric(as.character(order_list_pred$network_index))
  network_index_pred <- as.numeric(as.character(order_list_pred$network_index_pred))
  network_index_vals <- sort(unique(c(network_index_obs, network_index_pred)))
  # network_index_integer <- seq_along(network_index_vals)

  dist_order <- order_list_pred$dist_order

  # order weights by dist_order
  additive_val <- as.numeric(ssn.object$obs[[additive]]) # remember a character here
  additive_val_order <- additive_val[dist_order]
  additive_pred_val <- as.numeric(ssn.object$preds[[newdata_name]][[additive]]) # remember a character here
  dist_pred_order <- order(network_index_pred, order_list_pred$pid_pred)
  additive_pred_val_order <- additive_pred_val[dist_pred_order]

  # make additive unmasked
  additive_pred_matlist <- lapply(network_index_vals, function(x) {
    vals_obs <- network_index_obs == x
    ind_obs <- which(vals_obs[dist_order])
    n_obs <- length(ind_obs)
    addfval_obs <- additive_val_order[ind_obs]
    vals_pred <- network_index_pred == x
    ind_pred <- which(vals_pred[dist_pred_order])
    n_pred <- length(ind_pred)
    addfval_pred <- additive_pred_val_order[ind_pred]
    additive_obs_val <- do.call(cbind, replicate(n_pred, addfval_obs, FALSE))
    additive_pred_val <- do.call(cbind, replicate(n_obs, addfval_pred, simplify = FALSE))
    if (n_obs != 0 && n_pred != 0) {
      w_pred_val <- pmin(additive_obs_val, t(additive_pred_val)) / pmax(additive_obs_val, t(additive_pred_val))
    } else {
      w_pred_val <- Matrix::Matrix(0, nrow = n_obs, ncol = n_pred)
    }
    Matrix::Matrix(sqrt(w_pred_val), sparse = TRUE)
  })

  # make w
  w_pred_matlist_val <- mapply(
    FUN = function(additive, b, m) {
      if (NROW(additive) > 0) {
        return(additive * (b == 0) * m) # b == 0 is flow connected)
      } else {
        return(additive) # return zero matrix if it is there
      }
    },
    additive = additive_pred_matlist,
    b = b_pred_matlist,
    m = mask_pred_matlist,
    SIMPLIFY = FALSE
  )
}
