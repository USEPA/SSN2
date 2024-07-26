#' Get the distance matrix oibject
#'
#' @param ssn.object SSN object.
#' @param initial_object Initial value object.
#' @param additive Name of the additive function value column.
#' @param anisotropy Whether there is anisotropy.
#'
#' @return A distance matrix object that contains various distance matrices used in modeling.
#' @noRd
get_dist_object <- function(ssn.object, initial_object, additive, anisotropy) {
  # get netgeom
  netgeom <- ssn_get_netgeom(ssn.object$obs, reformat = TRUE)

  # get network index
  network_index <- netgeom$NetworkID

  # get pid
  pid <- netgeom$pid # not needed now but can reorder by it later

  # distance order
  dist_order <- order(network_index, pid)

  # inverse of distance order (i.e., original data order)
  inv_dist_order <- order(dist_order)
  # get back pid by inv_dist_order[dist_order]

  # create "order" list
  order_list <- list(
    network_index = network_index,
    pid = pid,
    dist_order = dist_order,
    inv_dist_order = inv_dist_order
  )

  # get list of distance matrices in order of the original data
  dist_matlist <- get_dist_matlist(
    ssn.object, initial_object, additive,
    order_list
  )



  # see whether euclid is none to avoid unnecessary computations
  euclid_none <- inherits(initial_object$euclid_initial, "euclid_none")
  if (euclid_none) {
    # not needed if no euclid covariance
    dist_matlist <- c(dist_matlist, list(euclid_mat = NULL))
    # dist_matlist$euclid_matix <- NULL does not return anything
  } else {
    # store euclid coordinates (in original data order)
    obs_coords <- sf::st_coordinates(ssn.object$obs)

    # check for anisotropy and store accordingly
    if (anisotropy) {
      # store as vectors with drop
      dist_matlist$.xcoord <- obs_coords[, 1, drop = TRUE]
      dist_matlist$.ycoord <- obs_coords[, 2, drop = TRUE]
    } else {
      # compute distance matrix
      dist_matlist$euclid_mat <- Matrix::Matrix(as.matrix(dist(obs_coords)), sparse = TRUE)
    }
  }

  # append lists
  dist_object <- c(dist_matlist, order_list)

  # return distance object
  dist_object
}

# a vectorized version of get_dist_object
get_dist_object_oblist <- function(dist_object, observed_index, local_index) {
  # find names of the distance object
  dist_object_names <- names(dist_object)
  # apply appropriate subsetting to each element and rename
  dist_object <- lapply(dist_object_names, subset_dist_object, dist_object, observed_index)
  names(dist_object) <- dist_object_names

  # sort by each local index
  unq_local_index_sort <- sort(unique(local_index))

  # if there is only one local index, create a template list with one element
  if (length(unq_local_index_sort) == 1) {
    dist_object <- list(dist_object)
    names(dist_object) <- unq_local_index_sort
  } else { # commented out w/ local not available
    # if there is more than one local index, create a template list with many elements
    # dist_object <- lapply(unq_local_index_sort, function(x) {
    #   index_val <- which(local_index == x)
    #   new_val <- lapply(dist_object, subset_dist_object, index_val)
    # })
    # names(dist_object) <- names(split(unq_local_index_sort, unq_local_index_sort))
  }
  dist_object_oblist <- dist_object
  dist_object_oblist
}

#' Get list of full distance matrices
#'
#' @param ssn.object SSN object.
#' @param initial_object Initial value object.
#' @param additive Name of the additive function value column.
#' @param order_list A list of order by pid and network.
#'
#' @noRd
get_dist_matlist <- function(ssn.object, initial_object, additive,
                             order_list) {
  network_index <- order_list$network_index
  dist_order <- order_list$dist_order
  inv_dist_order <- order_list$inv_dist_order

  # see whether tailup and taildown are none to avoid unnecssary computations
  tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
  taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")

  # return all NULL if they are both none (no stream distance needed)
  if (tailup_none && taildown_none) {
    dist_matlist <- list(
      distjunc_matlist = NULL,
      mask_matlist = NULL,
      a_matlist = NULL,
      b_matlist = NULL,
      hydro_matlist = NULL,
      w_matlist = NULL
    )
  } else {
    # otherwise

    # get dist junction matrices as a list (for efficiency, do things
    # network by network and then combine so zeroes populate accordingly)
    distjunc_matlist <- get_distjunc_matlist(order_list$network_index, ssn.object)

    # get other matrices as a list
    dist_matlist <- list(
      distjunc_matlist = distjunc_matlist,
      mask_matlist = get_mask_matlist(distjunc_matlist),
      a_matlist = get_a_matlist(distjunc_matlist),
      b_matlist = get_b_matlist(distjunc_matlist),
      hydro_matlist = get_hydro_matlist(distjunc_matlist)
    )

    # if only taildown covariacne, do not need additive matrix
    if (tailup_none) {
      # store as single sparse Matrix
      dist_matlist <- list(
        distjunc_mat = Matrix::bdiag(dist_matlist$distjunc_matlist),
        mask_mat = Matrix::bdiag(dist_matlist$mask_matlist),
        a_mat = Matrix::bdiag(dist_matlist$a_matlist),
        b_mat = Matrix::bdiag(dist_matlist$b_matlist),
        hydro_mat = Matrix::bdiag(dist_matlist$hydro_matlist)
      )

      # reorder in terms of the original data
      dist_matlist <- lapply(dist_matlist, function(x) {
        x[inv_dist_order, inv_dist_order, drop = FALSE]
      })

      # store additive matrix as NULL
      dist_matlist <- c(dist_matlist, list(w_mat = NULL))
    } else {
      # compute additive matrix
      dist_matlist$w_matlist <- get_w_matlist(
        ssn.object,
        additive,
        order_list$network_index,
        order_list$dist_order,
        dist_matlist$b_matlist,
        dist_matlist$mask_matlist
      )

      # store as single sparse Matrix
      dist_matlist <- list(
        distjunc_mat = Matrix::bdiag(dist_matlist$distjunc_matlist),
        mask_mat = Matrix::bdiag(dist_matlist$mask_matlist),
        a_mat = Matrix::bdiag(dist_matlist$a_matlist),
        b_mat = Matrix::bdiag(dist_matlist$b_matlist),
        hydro_mat = Matrix::bdiag(dist_matlist$hydro_matlist),
        w_mat = Matrix::bdiag(dist_matlist$w_matlist)
      )

      # reorder in terms of the original data
      dist_matlist <- lapply(dist_matlist, function(x) {
        x[inv_dist_order, inv_dist_order, drop = FALSE]
      })
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
get_distjunc_matlist <- function(network_index, ssn.object, newdata_name = NULL) {
  if (is.null(newdata_name)) {
    ext <- "obs"
  } else {
    ext <- newdata_name # for block Kriging prediction
  }

  # get unique network index vals (from STARS)
  network_index_vals <- sort(as.numeric(as.character(unique(network_index))))

  # find distance junction list by iterating through each network
  distjunc_list <- lapply(network_index_vals, function(x) {
    # on the disk, distance matrices are stored by network
    workspace_name <- paste("dist.net", x, ".RData", sep = "")
    # path to the distance matrices on disk
    path <- file.path(ssn.object$path, "distance", ext, workspace_name)
    # check to see if the file exists on the disk
    if (!file.exists(path)) {
      stop("Unable to locate required distance matrix", call. = FALSE)
    }
    # some code to read from disk (binary representation)
    file_handle <- file(path, open = "rb")
    # get the distance-to-nearest junction matrix
    dist_mat <- unserialize(file_handle)
    # get pid order
    # could also get this as dist_order[network_index == x]
    pid_order <- order(as.numeric(rownames(dist_mat)))
    # close the file on disk
    close(file_handle)
    # get distance juncture
    distjunc <- dist_mat[pid_order, pid_order, drop = FALSE]
  })
}


#' Subset relevant objects in the distance object
#'
#' @param dist_object_name Name of the relevant element in the distance object
#' @param dist_object Distance object
#' @param index Index of the distance matrix elements for ordering
#'
#' @noRd
subset_dist_object <- function(dist_object_name, dist_object, index) {
  # store the distance object element
  dist_object_element <- dist_object[[dist_object_name]]

  # if the element is a matrix, perform matrix subsetting (using NULL where appropriate)
  # otherwise, perform vector subsetting (using NULL where appropriate)

  if (is.null(dist_object_element)) {
    dist_sub <- NULL
  } else {
    if (dist_object_name %in% c("distjunc_mat", "mask_mat", "a_mat", "b_mat", "hydro_mat", "w_mat", "euclid_mat")) {
      dist_sub <- dist_object[[dist_object_name]][index, index, drop = FALSE]
    } else if (dist_object_name %in% c("network_index", "pid", "dist_order", "inv_dist_order", ".xcoord", ".ycoord")) {
      dist_sub <- dist_object[[dist_object_name]][index]
    } else {
      stop("Invalid distance object names", call. = FALSE)
    }
  }
  dist_sub
}

#' Find the mask matrix list for each network
#'
#' @param distjunc_list Create list of mask matrices (which zero out tailup
#'   and taildown covariance for observations on separate networks).
#'
#' @noRd
get_mask_matlist <- function(distjunc_list) {
  mask_list <- lapply(distjunc_list, function(x) {
    n_i <- dim(x)[[1]]
    Matrix::Matrix(1, nrow = n_i, ncol = n_i)
  })
}

#' Find the a matrix list for each network
#'
#' @param distjunc_list Create list of a matrices (largest stream distance) for each network.
#'
#' @noRd
get_a_matlist <- function(distjunc_list) {
  # find a matrix list
  a_matlist <- lapply(distjunc_list, function(x) {
    x_matrix <- as.matrix(x)
    Matrix::Matrix(pmax(x_matrix, t(x_matrix)), sparse = TRUE)
  })
}

#' Find the b matrix list for each network
#'
#' @param distjunc_list Create list of a matrices (smallest stream distance) for each network.
#'
#' @noRd
get_b_matlist <- function(distjunc_list) {
  # find a matrix list
  b_matlist <- lapply(distjunc_list, function(x) {
    x_matrix <- as.matrix(x)
    Matrix::Matrix(pmin(x_matrix, t(x_matrix)), sparse = TRUE)
  })
}

#' Find the hydrologic distance matrix list for each network
#'
#' @param distjunc_list Create list of matrices (hydrologic stream distance) for each network.
#'
#' @noRd
get_hydro_matlist <- function(distjunc_list) {
  # find hydro matrix list
  hydro_matlist <- lapply(distjunc_list, function(x) x + t(x))
}

#' Find the w matrix list for each network
#'
#' @param distjunc_list Create list of w matrices (additive function weight) for each network.
#'
#' @noRd
get_w_matlist <- function(ssn.object, additive, network_index, dist_order, b_matlist, mask_matlist, newdata_name = NULL) {
  if (is.null(newdata_name)) {
    additive_val <- as.numeric(ssn.object$obs[[additive]]) # remember a character here
  } else {
    additive_val <- as.numeric(ssn.object$preds[[newdata_name]][[additive]]) # for block kriging
  }

  # order weights by dist_order
  additive_val_order <- additive_val[dist_order]

  # make list
  additive_list <- split(additive_val_order, network_index[dist_order])

  # make additive unmasked
  additive_matlist <- lapply(additive_list, function(x) {
    n_i <- length(x)
    additive_matrix_val <- replicate(n_i, x, simplify = "array")
    w_matrix_val <- pmin(additive_matrix_val, t(additive_matrix_val)) / pmax(additive_matrix_val, t(additive_matrix_val))
    Matrix::Matrix(sqrt(w_matrix_val), sparse = TRUE)
  })

  # apply flow connected and masking requirements
  w_matlist <- mapply(
    FUN = function(additive, b, m) additive * (b == 0) * m, # b == 0 is flow connected
    additive = additive_matlist,
    b = b_matlist,
    m = mask_matlist,
    SIMPLIFY = FALSE
  )
}
