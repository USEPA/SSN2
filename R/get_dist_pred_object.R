#' Get prediction distance object
#'
#' @param object Data object.
#' @param newdata_name Name of the prediction data set.
#' @param initial_object Initial value object.
#'
#' @noRd
get_dist_pred_object <- function(object, newdata_name, initial_object) {
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



  # get list of prediction distance matrices in order of the original data
  dist_pred_matlist <- get_dist_pred_matlist(
    object$ssn.object, newdata_name, initial_object, object$additive,
    order_list_pred
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
                                  order_list_pred) {

  # store network indices and orders
  network_index <- order_list_pred$network_index
  dist_order <- order_list_pred$dist_order
  inv_dist_order <- order_list_pred$inv_dist_order
  inv_dist_order_pred <- order_list_pred$inv_dist_order_pred

  # see whether tailup and taildown are none to avoid unnecssary computations
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

    # get dist junction matrices as a list (for efficiency, do things
    # network by network and then combine so zeroes populate accordingly)
    distjunc_pred_matlist <- get_distjunc_pred_matlist(ssn.object, newdata_name, order_list_pred)

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

      # get distance matrices in pid (data) order (rows by data order and
      # columns by prediction data order)
      dist_pred_matlist <- lapply(dist_pred_matlist, function(x) {
        preds_val <- tryCatch(x[inv_dist_order, inv_dist_order_pred, drop = FALSE],
          error = function(e) x[inv_dist_order_pred, inv_dist_order, drop = FALSE]
        )
        preds_val
      })

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

      # get distance matrices in pid (data) order (rows by data order and
      # columns by prediction data order)
      dist_pred_matlist <- lapply(dist_pred_matlist, function(x) {
        preds_val <- tryCatch(x[inv_dist_order, inv_dist_order_pred, drop = FALSE],
          error = function(e) x[inv_dist_order_pred, inv_dist_order, drop = FALSE]
        )
        preds_val
      })
    }
  }
  # return distance prediction object
  dist_pred_matlist
}

#' Get list of prediction distance junction matrices
#'
#' @param ssn.object SSN object.
#' @param newdata_name Name of the prediction data set.
#' @param order_list_pred The order for observations in the prediction data set.
#'
#' @noRd
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
        # some code to read from disk (binary representation)
        file_handle <- file(path, open = "rb")
        # get the distance-to-nearest junction matrix
        distmat <- unserialize(file_handle)
        # close the file on disk
        close(file_handle)

        # find observations that are used to build model
        which_obs <- rownames(distmat) %in% network_pid_obs[network_index_obs == x]
        # find prediction observations
        which_pred <- !which_obs
        distmata <- distmat[which_obs, which_pred, drop = FALSE]
        distmatb <- distmat[which_pred, which_obs, drop = FALSE]
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
        # distance matrix a
        file_handle <- file(path.a, open = "rb")
        distmata <- unserialize(file_handle)
        close(file_handle)
        # distance matrix b
        file_handle <- file(path.b, open = "rb")
        distmatb <- unserialize(file_handle)
        close(file_handle)

        # only keep observed from ssn object
        which_obs <- rownames(distmata) %in% network_pid_obs[network_index_obs == x]
        distmata <- distmata[which_obs, , drop = FALSE]
        distmatb <- distmatb[, which_obs, drop = FALSE]
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


#' Get mask prediction distance matrices
#'
#' @param distjunc_pred_matlist
#'
#' @noRd
get_mask_pred_matlist <- function(distjunc_pred_matlist) {
  mask_pred_list <- lapply(distjunc_pred_matlist$distjunca, function(x) {
    Matrix::Matrix(1, nrow = dim(x)[1], ncol = dim(x)[2], sparse = TRUE)
  })
}

#' Get a prediction distance matrices
#'
#' @param distjunc_pred_matlist
#'
#' @noRd
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

#' Get b prediction distance matrices
#'
#' @param distjunc_pred_matlist
#'
#' @noRd
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

#' Get hydrologic prediction distance matrices
#'
#' @param distjunc_pred_matlist
#'
#' @noRd
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

#' Get w prediction distance matrices
#'
#' @param distjunc_pred_matlist
#'
#' @noRd
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
