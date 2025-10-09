#' Get the distance matrix oibject
#'
#' @param ssn.object SSN object.
#' @param initial_object Initial value object.
#' @param additive Name of the additive function value column.
#' @param anisotropy Whether there is anisotropy.
#'
#' @return A distance matrix object that contains various distance matrices used in modeling.
#' @noRd
get_dist_object_bigdata <- function(ssn.object, initial_object, additive, anisotropy, local_index, observed_index) {
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

  if (is.null(additive)) {
    additive_val <- NULL
    order_list <- data.frame(
      network_index = network_index,
      pid = pid,
      dist_order = dist_order,
      inv_dist_order = inv_dist_order,
      observed_index
    )
  } else {
    additive_val <- as.numeric(ssn.object$obs[[additive]])
    order_list <- data.frame(
      network_index = network_index,
      pid = pid,
      dist_order = dist_order,
      inv_dist_order = inv_dist_order,
      additive_val = additive_val,
      observed_index
    )
  }



  obs_index <- order_list$observed_index # may need to be sorted
  order_list <- order_list[obs_index, , drop = FALSE]
  local_index_orig <- local_index
  order_bigdata <- order(local_index, order_list$network_index, order_list$pid)
  local_index <- local_index[order_bigdata]
  order_list <- order_list[order_bigdata, , drop = FALSE]
  order_list <- split(order_list, local_index)

  # get list of distance matrices in order of the original data
  dist_matlist <- get_dist_matlist_bigdata(
    ssn.object, initial_object,
    order_list
  )

  # see whether euclid is none to avoid unnecessary computations
  euclid_none <- inherits(initial_object$euclid_initial, "euclid_none")
  if (euclid_none) {
    # not needed if no euclid covariance
    dist_matlist <- lapply(dist_matlist, function(x) c(x, list(euclid_mat = NULL)))
    # dist_matlist$euclid_matix <- NULL does not return anything
  } else {
    # store euclid coordinates (in original data order but then changed)
    obs_coords <- sf::st_coordinates(ssn.object$obs)
    obs_coords <- obs_coords[obs_index, , drop = FALSE]
    obs_coords <- obs_coords[order_bigdata, , drop = FALSE]
    obs_coords_list <- split.data.frame(obs_coords, local_index) # local index should already be sorted

    dist_matlist <- mapply(d = dist_matlist, c = obs_coords_list, function(d, c) {
      if (anisotropy) {
        # store as vectors with drop
        d$.xcoord <- c[, 1, drop = TRUE]
        d$.ycoord <- c[, 2, drop = TRUE]
      } else {
        # compute distance matrix
        d$euclid_mat <- Matrix::Matrix(as.matrix(dist(c)), sparse = TRUE)
      }
      d
    }, SIMPLIFY = FALSE)
  }

  # append lists
  dist_object <- list(dist_matlist = dist_matlist, order_list = order_list)

  # return distance object
  dist_object
}

# a vectorized version of get_dist_object
get_dist_object_oblist_bigdata <- function(dist_object) {

  dist_object_names <- names(dist_object$dist_matlist)
  dist_object <- mapply(d = dist_object$dist_matlist, o = dist_object$order_list, FUN = function(d, o) {
    dist_object_name <- names(d)
    dist_sub <- lapply(dist_object_name, function(x) subset_dist_object_bigdata(x, d, o$observed_index))
    names(dist_sub) <- dist_object_name
    order_sub_name <- names(o)
    order_sub <- lapply(order_sub_name, function(x) subset_dist_object_bigdata(x, o, o$observed_index))
    names(order_sub) <- order_sub_name
    order_sub <- order_sub[-which(names(order_sub) == "observed_index")]
    c(dist_sub, order_sub)
  }, SIMPLIFY = FALSE)
  names(dist_object) <- dist_object_names
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
get_dist_matlist_bigdata <- function(ssn.object, initial_object,
                             order_list) {

  # see whether tailup and taildown are none to avoid unnecssary computations
  tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
  taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")

  # return all NULL if they are both none (no stream distance needed)
  if (tailup_none && taildown_none) {
    dist_matlist <- lapply(order_list, function(x) {
      list(
        distjunc_mat = NULL,
        mask_mat = NULL,
        a_mat = NULL,
        b_mat = NULL,
        hydro_mat = NULL,
        w_mat = NULL
      )
    })
  } else {
    # otherwise

    # get dist junction matrices as a list (for efficiency, do things
    # network by network and then combine so zeroes populate accordingly)
    distjunc_matlist <- lapply(order_list, function(x) get_distjunc_matlist_bigdata(x$network_index, x$pid, ssn.object))


    dist_matlist <- lapply(names(order_list), function(x) {
      dist_matlist <- list(
        distjunc_mat = Matrix::bdiag(distjunc_matlist[[x]]),
        mask_mat = Matrix::bdiag(get_mask_matlist(distjunc_matlist[[x]])),
        a_mat = Matrix::bdiag(get_a_matlist(distjunc_matlist[[x]])),
        b_mat = Matrix::bdiag(get_b_matlist(distjunc_matlist[[x]])),
        hydro_mat = Matrix::bdiag(get_hydro_matlist(distjunc_matlist[[x]]))
      )
    })
    names(dist_matlist) <- names(order_list)

    # # get other matrices as a list
    # dist_matlist <- list(
    #   distjunc_matlist = dist_matlist$distjunc_matlist,
    #   mask_matlist = dist_matlist$mask_matrix,
    #   a_matlist = get_a_matlist(distjunc_matlist),
    #   b_matlist = get_b_matlist(distjunc_matlist),
    #   hydro_matlist = get_hydro_matlist(distjunc_matlist)
    # )
    #
    # dist_matlist <- lapply(names(order_list), function(x) {
    #   dist_matlist <- list(
    #     distjunc_mat = dist_matlist$distjunc_matlist[[x]],
    #     mask_mat = dist_matlist$mask_matlist[[x]],
    #     a_mat = dist_matlist$a_matlist[[x]],
    #     b_mat = dist_matlist$b_matlist[[x]],
    #     hydro_mat = dist_matlist$hydro_matlist[[x]]
    #   )
    # })
    names(dist_matlist) <- names(order_list)

    # if only taildown covariacne, do not need additive matrix
    if (tailup_none) {
      # store additive matrix as NULL
      dist_matlist <- lapply(dist_matlist, function(x) c(x, list(w_mat = NULL)))
    } else {
      # compute additive matrix
      w_mat <- get_w_matlist_bigdata(
        ssn.object,
        order_list,
        dist_matlist
      )
      dist_matlist <- mapply(FUN = function(d, w) {
        c(d, w_mat = w)
      }, d = dist_matlist, w = w_mat, SIMPLIFY = FALSE)
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
get_distjunc_matlist_bigdata <- function(network_index, pid, ssn.object) {

  ext <- "obs"

  # get unique network index vals (from STARS)
  network_index_vals <- sort(as.numeric(as.character(unique(network_index))))

  # find distance junction list by iterating through each network
  distjunc_list <- lapply(network_index_vals, function(x) {
    # on the disk, distance matrices are stored by network
    workspace_name <- paste("dist.net", x, ".bmat", sep = "")
    # path to the distance matrices on disk
    path <- file.path(ssn.object$path, "distance", ext, workspace_name)
    # check to see if the file exists on the disk
    if (!file.exists(path)) {
      stop("Unable to locate required distance matrix", call. = FALSE)
    }
    distjunc_fm <- fm.open(path)
    rownames_val <- rownames(distjunc_fm)
    pid_index <- rownames_val %in% as.character(pid)
    dist_mat <- distjunc_fm[pid_index, pid_index]
    rownames(dist_mat) <- rownames_val[pid_index]
    # get pid order
    # could also get this as dist_order[network_index == x]
    pid_order <- order(as.numeric(rownames(dist_mat)))
    # get distance juncture
    distjunc <- dist_mat[pid_order, pid_order, drop = FALSE]
    close(distjunc_fm)
    distjunc
  })
}


#' Subset relevant objects in the distance object
#'
#' @param dist_object_name Name of the relevant element in the distance object
#' @param dist_object Distance object
#' @param index Index of the distance matrix elements for ordering
#'
#' @noRd
subset_dist_object_bigdata <- function(dist_object_name, dist_object, index) {
  # store the distance object element

  dist_object_element <- dist_object[[dist_object_name]]

  # if the element is a matrix, perform matrix subsetting (using NULL where appropriate)
  # otherwise, perform vector subsetting (using NULL where appropriate)

  if (is.null(dist_object_element)) {
    dist_sub <- NULL
  } else {
    if (dist_object_name %in% c("distjunc_mat", "mask_mat", "a_mat", "b_mat", "hydro_mat", "w_mat", "euclid_mat")) {
      dist_sub <- dist_object[[dist_object_name]][index, index, drop = FALSE]
    } else if (dist_object_name %in% c("network_index", "pid", "dist_order", "inv_dist_order", ".xcoord", ".ycoord", "additive_val", "observed_index")) {
      dist_sub <- dist_object[[dist_object_name]][index]
    } else {
      stop("Invalid distance object names", call. = FALSE)
    }
  }
  dist_sub
}
#' Find the w matrix list for each network
#'
#' @param distjunc_list Create list of w matrices (additive function weight) for each network.
#'
#' @noRd
get_w_matlist_bigdata <- function(ssn.object, order_list, dist_matlist) {

  # make additive unmasked
  additive_matlist <- lapply(order_list, function(x) {
    x <- x$additive_val
    n_i <- length(x)
    additive_matrix_val <- replicate(n_i, x, simplify = "array")
    w_matrix_val <- pmin(additive_matrix_val, t(additive_matrix_val)) / pmax(additive_matrix_val, t(additive_matrix_val))
    Matrix::Matrix(sqrt(w_matrix_val), sparse = TRUE)
  })

  # apply flow connected and masking requirements
  w_matlist <- mapply(
    FUN = function(additive, d) additive * (d$b_mat == 0) * d$mask_mat, # b == 0 is flow connected
    additive = additive_matlist,
    d = dist_matlist,
    SIMPLIFY = FALSE
  )
}
