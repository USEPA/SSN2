#' Get prediction by prediction distance matrix for block Kriging
#'
#' @param object Model object.
#' @param newdata_name Name of prediction data set.
#' @param initial_object Initial value object.
#'
#' @noRd
get_dist_predbk_object <- function(object, newdata_name, initial_object) {
  # get netgeom
  netgeom <- ssn_get_netgeom(object$ssn.object$preds[[newdata_name]], reformat = TRUE)

  # get network index
  network_index <- netgeom$NetworkID

  # get pid
  pid <- netgeom$pid

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
  dist_matlist <- get_dist_predbk_matlist(
    object$ssn.object, newdata_name, initial_object, object$additive,
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
    pred_coords <- sf::st_coordinates(object$ssn.object$preds[[newdata_name]])

    # check for anisotropy and store accordingly
    if (object$anisotropy) {
      # store as vectors with drop
      dist_matlist$.xcoord <- pred_coords[, 1, drop = TRUE]
      dist_matlist$.ycoord <- pred_coords[, 2, drop = TRUE]
    } else {
      # compute distance matrix
      dist_matlist$euclid_mat <- Matrix::Matrix(as.matrix(dist(pred_coords)), sparse = TRUE)
    }
  }

  # append lists
  dist_object <- c(dist_matlist, order_list)

  # return distance object
  dist_object
}

# vectorized version of get_dist_predbk_object
get_dist_predbk_matlist <- function(ssn.object, newdata_name, initial_object, additive,
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
    distjunc_matlist <- get_distjunc_matlist(order_list$network_index, ssn.object, newdata_name)

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
        dist_matlist$mask_matlist,
        newdata_name = newdata_name
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
