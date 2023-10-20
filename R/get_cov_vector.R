get_cov_vector <- function(params_object, dist_pred_object, data, newdata, partition_factor = NULL, anisotropy) {
  tailup_none <- inherits(params_object$tailup, "tailup_none")
  taildown_none <- inherits(params_object$taildown, "taildown_none")
  euclid_none <- inherits(params_object$euclid, "euclid_none")

  if (tailup_none && taildown_none && euclid_none) {
    cov_vector <- Matrix::Matrix(0, nrow = NROW(newdata), ncol = NROW(data), sparse = TRUE)
  } else {
    cov_vector <- cov_vector(params_object$tailup, dist_pred_object) +
      cov_vector(params_object$taildown, dist_pred_object) +
      cov_vector(params_object$euclid, dist_pred_object, anisotropy)
    if (!is.null(params_object$randcov)) {
      cov_vector <- cov_vector + randcov_vector(params_object$randcov, data, newdata)
    }
    if (!is.null(partition_factor)) {
      cov_vector <- cov_vector * partition_vector(partition_factor, data, newdata)
    }
  }
  # Matrix::Matrix(cov_vector, sparse = TRUE)
  cov_vector
}

get_euclid_pred <- function(params, dist_pred_object, anisotropy) {
  if (anisotropy) {
    new_coords_observed <- transform_anis(
      dist_pred_object$.obs_xcoord,
      dist_pred_object$.obs_ycoord,
      rotate = params[["rotate"]],
      scale = params[["scale"]]
    )

    new_coords_pred <- transform_anis(
      dist_pred_object$.pred_xcoord,
      dist_pred_object$.pred_ycoord,
      rotate = params[["rotate"]],
      scale = params[["scale"]]
    )

    dist_vector_x <- outer(X = new_coords_observed$xcoord_val, Y = new_coords_pred$xcoord_val, FUN = function(X, Y) (X - Y)^2)
    dist_vector_y <- outer(X = new_coords_observed$ycoord_val, Y = new_coords_pred$ycoord_val, FUN = function(X, Y) (X - Y)^2)
    dist_vector <- sqrt(dist_vector_x + dist_vector_y)
    euclid_pred <- Matrix::Matrix(dist_vector, sparse = TRUE)
    # transpose to match with final transpose in get_dist_pred_object
    euclid_pred <- t(euclid_pred)
  } else {
    euclid_pred <- dist_pred_object$euclid_pred_mat
  }
  euclid_pred
}
