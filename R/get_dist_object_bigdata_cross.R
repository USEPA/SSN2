get_dist_object_bigdata_cross <- function(d1, d2, params_object, data_object) {


  # get list of distance matrices in order of the original data
  dist_matlist <- get_dist_matlist_bigdata_cross(d1, d2, params_object, data_object)

  # see whether euclid is none to avoid unnecessary computations
  euclid_none <- inherits(params_object$euclid, "euclid_none")
  if (euclid_none) {
    # not needed if no euclid covariance
    dist_matlist <-c(dist_matlist, list(euclid_mat = NULL))
    # dist_matlist$euclid_matix <- NULL does not return anything
  } else {
    # store euclid coordinates (in original data order)
    d1_obs_coords <- sf::st_coordinates(d1)
    d2_obs_coords <- sf::st_coordinates(d2)
    if (data_object$anisotropy) {
      dist_matlist$.xcoord_d1 <- d1_obs_coords[, 1, drop = TRUE]
      dist_matlist$.ycoord_d1 <- d1_obs_coords[, 2, drop = TRUE]
      dist_matlist$.xcoord_d2 <- d2_obs_coords[, 1, drop = TRUE]
      dist_matlist$.ycoord_d2 <- d2_obs_coords[, 2, drop = TRUE]
    } else {
      x_diffs <- outer(d1_obs_coords[, 1, drop = TRUE], d2_obs_coords[, 1, drop = TRUE], "-")
      y_diffs <- outer(d1_obs_coords[, 2, drop = TRUE], d2_obs_coords[, 2, drop = TRUE], "-")
      euclid_mat <- sqrt(x_diffs^2 + y_diffs^2)
      dist_matlist$euclid_mat <- Matrix::Matrix(euclid_mat, sparse = TRUE)
    }
  }

  # append lists
  dist_object_bigdata_cross_oblist <- dist_matlist

  # return distance object
  dist_object_bigdata_cross_oblist
}

get_dist_matlist_bigdata_cross <- function(d1, d2, params_object, data_object) {


  d1_netgeom <- ssn_get_netgeom(d1)
  d2_netgeom <- ssn_get_netgeom(d2)
  d1_network_index <- d1_netgeom$NetworkID
  d2_network_index <- d2_netgeom$NetworkID
  d1_pid <- d1_netgeom$pid
  d2_pid <- d2_netgeom$pid

  if (is.null(data_object$additive)) {
    d1_additive_val <- NULL
    d2_additive_val <- NULL
  } else {
    d1_additive_val <- d1[[data_object$additive]]
    d2_additive_val <- d2[[data_object$additive]]
  }



  # see whether tailup and taildown are none to avoid unnecssary computations
  tailup_none <- inherits(params_object$tailup, "tailup_none")
  taildown_none <- inherits(params_object$taildown, "taildown_none")

  # return all NULL if they are both none (no stream distance needed)
  if (tailup_none && taildown_none) {
    dist_matlist <- list(
      distjunc_mat = NULL,
      mask_mat = NULL,
      a_mat = NULL,
      b_mat = NULL,
      hydro_mat = NULL,
      w_mat = NULL
    )
  } else {
    # otherwise

    # get dist junction matrices as a list (for efficiency, do things
    # network by network and then combine so zeroes populate accordingly)
    distjunc_matlist1 <- get_distjunc_matlist_bigdata_cross(
      d1_network_index, d1_pid, d2_network_index, d2_pid, data_object$ssn.object
    )

    # to get the transpose
    distjunc_matlist2 <- get_distjunc_matlist_bigdata_cross(
      d2_network_index, d2_pid, d1_network_index, d1_pid, data_object$ssn.object
    )

    # turn these into a matrix object when its 1x1 the get_a, get_b below fail
    dist_matlist <- list(
      distjunc_mat = Matrix::Matrix(distjunc_matlist1, sparse = TRUE),
      mask_mat = Matrix::Matrix(get_mask_matlist_bigdata_cross(d1_network_index, d2_network_index, distjunc_matlist1), sparse = TRUE),
      a_mat = Matrix::Matrix(get_a_matlist_bigdata_cross(distjunc_matlist1, distjunc_matlist2), sparse = TRUE),
      b_mat = Matrix::Matrix(get_b_matlist_bigdata_cross(distjunc_matlist1, distjunc_matlist2), sparse = TRUE),
      hydro_mat = Matrix::Matrix(get_hydro_matlist_bigdata_cross(distjunc_matlist1, distjunc_matlist2), sparse = TRUE)
    )

    # if only taildown covariacne, do not need additive matrix
    if (tailup_none) {
      # store additive matrix as NULL
      dist_matlist <- c(dist_matlist, list(w_mat = NULL))
    } else {
      dist_matlist$w_mat <- get_w_matlist_bigdata_cross(d1_additive_val, d2_additive_val, dist_matlist$b_mat, dist_matlist$mask_mat)
    }
  }
  dist_matlist
}

#' Get distance to the nearest junction matrix
#'
#' @param network_index Network index
#' @param ssn.object SSN object
#' @param newdata_name Name of the newdata matrix (if relevant)
#'
#' @noRd
get_distjunc_matlist_bigdata_cross <- function(d1_network_index, d1_pid, d2_network_index, d2_pid, ssn.object) {

  # REMEMBER OBDATA_LIST ALREADY ORDERED BY SPLIT ID, NETWORK ID, PID
  ext <- "obs"

  # get unique network index vals (from STARS)
  # only looking for unique networks on rows as cols will zero out from masking anyways
  network_index_vals <- sort(as.numeric(as.character(unique(d1_network_index))))

  # regular matrix object as this causes 1x1 matrices to fail in spatial indexing
  # distjunc <- Matrix::Matrix(0, nrow = length(d1_pid), ncol = length(d2_pid), sparse = TRUE)
  distjunc <- matrix(0, nrow = length(d1_pid), ncol = length(d2_pid))
  # already ordered by PID
  rownames(distjunc) <- d1_pid
  colnames(distjunc) <- d2_pid
  # find distance junction list by iterating through each network
  for (x in network_index_vals) { # order by net ID
    # on the disk, distance matrices are stored by network
    workspace_name <- paste("dist.net", x, ".bmat", sep = "")
    # path to the distance matrices on disk
    path <- file.path(ssn.object$path, "distance", ext, workspace_name)
    # check to see if the file exists on the disk
    if (!file.exists(path)) {
      stop("Unable to locate required distance matrix", call. = FALSE)
    }
    distjunc_fm <- fm.open(path)
    rows_keep <- d1_pid[rownames(distjunc) %in% rownames(distjunc_fm)]
    cols_keep <- d2_pid[colnames(distjunc) %in% colnames(distjunc_fm)]
    rows_keep_distjunc_fm <- which(rownames(distjunc_fm) %in% rows_keep)
    cols_keep_distjunc_fm <- which(colnames(distjunc_fm) %in% cols_keep)
    if (length(rows_keep_distjunc_fm) > 0 && length(cols_keep_distjunc_fm)) {
      rows_keep_distjunc <- which(rownames(distjunc) %in% rows_keep)
      cols_keep_distjunc <- which(colnames(distjunc) %in% cols_keep)
      distjunc_fm_sub <- distjunc_fm[rows_keep_distjunc_fm, cols_keep_distjunc_fm]
      # order by PID
      distjunc_fm_sub <- distjunc_fm_sub[order(as.numeric(rows_keep)), order(as.numeric(cols_keep))]
      # put in distance junction matrix which is already ordered by local ID, net ID, PID
      distjunc[rows_keep_distjunc, cols_keep_distjunc] <- distjunc_fm_sub
    }
    close(distjunc_fm)
  }
  distjunc
}

get_mask_matlist_bigdata_cross <- function(d1_network_index, d2_network_index, distjunc_mat) {
  mask_mat <- outer(d1_network_index, d2_network_index, FUN = "==")
  rownames(mask_mat) <- rownames(distjunc_mat)
  colnames(mask_mat) <- colnames(distjunc_mat)
  1 * mask_mat
}

get_a_matlist_bigdata_cross <- function(distjunc_matlist1, distjunc_matlist2) {
  pmax(distjunc_matlist1, t(distjunc_matlist2))
}

get_b_matlist_bigdata_cross <- function(distjunc_matlist1, distjunc_matlist2) {
  pmin(distjunc_matlist1, t(distjunc_matlist2))
}

get_hydro_matlist_bigdata_cross <- function(distjunc_matlist1, distjunc_matlist2) {
  distjunc_matlist1 + t(distjunc_matlist2)
}

get_w_matlist_bigdata_cross <- function(d1_additive_val, d2_additive_val, b_matlist_bigdata_cross, mask_matlist_bigdata_cross) {
  d1_additive_matrix_val <- replicate(length(d2_additive_val), d1_additive_val, simplify = "array")
  d2_additive_matrix_val <- replicate(length(d1_additive_val), d2_additive_val, simplify = "array")
  w_matrix_val <- pmin(d1_additive_matrix_val, t(d2_additive_matrix_val)) / pmax(d1_additive_matrix_val, t(d2_additive_matrix_val))
  # sqrt(w_matrix_val) # bug obviously
  sqrt(w_matrix_val) * (b_matlist_bigdata_cross == 0) * mask_matlist_bigdata_cross
}
