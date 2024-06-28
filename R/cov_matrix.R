###############################################################################
# GENERIC SETUP
###############################################################################


#' Create a covariance matrix
#'
#' @param params Parameter object.
#' @param dist_object Distance matrix object.
#' @param ... Additional arguments
#'
#' @return A covariance matrix
#' @noRd
cov_matrix <- function(params, dist_object, ...) {
  UseMethod("cov_matrix", params)
}


# m is first for comp reasons (a sparse matrix)

###############################################################################
# TAILUP COVARIANCES
###############################################################################

#' @export
cov_matrix.tailup_linear <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  w <- dist_object$w_mat
  m <- dist_object$mask_mat
  dist_ratio <- h / params[["range"]]
  m * params[["de"]] * (1 - dist_ratio) * (dist_ratio <= 1) * w
}

#' @export
cov_matrix.tailup_spherical <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  w <- dist_object$w_mat
  m <- dist_object$mask_mat
  dist_ratio <- h / params[["range"]]
  m * params[["de"]] * (1 - 1.5 * dist_ratio + 0.5 * dist_ratio^3) * (dist_ratio <= 1) * w
}

#' @export
cov_matrix.tailup_exponential <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  w <- dist_object$w_mat
  m <- dist_object$mask_mat
  dist_ratio <- h / params[["range"]]
  m * params[["de"]] * exp(-dist_ratio) * w
}

#' @export
cov_matrix.tailup_mariah <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  w <- dist_object$w_mat
  m <- dist_object$mask_mat
  dist_ratio <- 90 * h / params[["range"]] # remove 90 for range that is not effective range
  middle <- params[["de"]] * log(dist_ratio + 1) / dist_ratio
  middle[h == 0] <- params[["de"]]
  m * middle * w
}

#' @export
cov_matrix.tailup_epa <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  w <- dist_object$w_mat
  m <- dist_object$mask_mat
  range <- params[["range"]]
  dist_ratio <- h / range
  f_eu <- 16 * range^3 + 17 * range^2 * h - 2 * range * h^2 - h^3
  m * params[["de"]] * (h - range)^2 * f_eu * (dist_ratio <= 1) / (16 * range^5) * w
}

#' @export
cov_matrix.tailup_none <- function(params, dist_object, ...) {
  0
}

###############################################################################
# TAILDOWN COVARIANCES
###############################################################################

#' @export
cov_matrix.taildown_linear <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  a <- dist_object$a_mat
  b <- dist_object$b_mat
  m <- dist_object$mask_mat
  flow_con <- b == 0
  dist_ratio_h <- h / params[["range"]]
  dist_ratio_a <- a / params[["range"]]
  dist_ratio_h_less1 <- dist_ratio_h <= 1
  dist_ratio_a_less1 <- dist_ratio_a <= 1

  h_cor_part <- (1 - dist_ratio_h) * dist_ratio_h_less1 * flow_con
  a_cor_part <- (1 - dist_ratio_a) * dist_ratio_a_less1 * (!flow_con)
  m * params[["de"]] * (h_cor_part + a_cor_part)
}

#' @export
cov_matrix.taildown_spherical <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  a <- dist_object$a_mat
  b <- dist_object$b_mat
  m <- dist_object$mask_mat
  flow_con <- b == 0
  dist_ratio_h <- h / params[["range"]]
  dist_ratio_a <- a / params[["range"]]
  dist_ratio_b <- b / params[["range"]]
  dist_ratio_h_less1 <- dist_ratio_h <= 1
  dist_ratio_a_less1 <- dist_ratio_a <= 1

  h_cor_part <- (1 - 1.5 * dist_ratio_h + 0.5 * dist_ratio_h^3) * dist_ratio_h_less1 * flow_con
  a_cor_part <- (1 - 1.5 * dist_ratio_b + 0.5 * dist_ratio_a) * (1 - dist_ratio_a)^2 * dist_ratio_a_less1 * (!flow_con)
  m * params[["de"]] * (h_cor_part + a_cor_part)
}

#' @export
cov_matrix.taildown_exponential <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  m <- dist_object$mask_mat
  dist_ratio <- h / params[["range"]]
  m * params[["de"]] * exp(-dist_ratio)
}

#' @export
cov_matrix.taildown_mariah <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  a <- dist_object$a_mat
  b <- dist_object$b_mat
  m <- dist_object$mask_mat
  flow_con <- b == 0
  dist_ratio_h <- 90 * h / params[["range"]]
  dist_ratio_a <- 90 * a / params[["range"]]
  dist_ratio_b <- 90 * b / params[["range"]]
  a_eq_b <- which(a == b)


  h_cor_part <- log(dist_ratio_h + 1) / dist_ratio_h * flow_con
  h_cor_part[h == 0] <- 1
  ab_cor_part <- (log(dist_ratio_a + 1) - log(dist_ratio_b + 1)) / (dist_ratio_a - dist_ratio_b)
  ab_cor_part[a_eq_b] <- 1 / (dist_ratio_a[a_eq_b] + 1) # without which a_eq_b to right of <- gets converted to dense format in Matrix?
  ab_cor_part <- ab_cor_part * (!flow_con)

  m * params[["de"]] * (h_cor_part + ab_cor_part)
}

#' @export
cov_matrix.taildown_epa <- function(params, dist_object, ...) {
  h <- dist_object$hydro_mat
  a <- dist_object$a_mat
  b <- dist_object$b_mat
  m <- dist_object$mask_mat
  flow_con <- b == 0
  range <- params[["range"]]
  dist_ratio_h <- h / range
  dist_ratio_a <- a / range
  dist_ratio_h_less1 <- dist_ratio_h <= 1
  dist_ratio_a_less1 <- dist_ratio_a <= 1

  f_eu <- 16 * range^3 + 17 * range^2 * h - 2 * range * h^2 - h^3
  f_ed <- 16 * range^3 + 17 * range^2 * a - 15 * range^2 * b - 20 * range * b^2 - 2 * range * a^2 + 10 * range * a * b + 5 * a^2 * b - a^3 - 10 * a * b^2
  h_cor_part <- (h - range)^2 * f_eu * dist_ratio_h_less1 / (16 * range^5) * flow_con
  a_cor_part <- (a - range)^2 * f_ed * dist_ratio_a_less1 / (16 * range^5) * (!flow_con)
  m * params[["de"]] * (h_cor_part + a_cor_part)
}

#' @export
cov_matrix.taildown_none <- function(params, dist_object, ...) {
  0
}

###############################################################################
# EUCLIDEAN COVARIANCES
###############################################################################

#' @export
cov_matrix.euclid_exponential <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- exp(-dist_ratio)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_spherical <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- (1 - (3 / 2) * dist_ratio + (1 / 2) * dist_ratio^3) * (dist_ratio <= 1)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_gaussian <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- exp(-dist_ratio^2)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_cosine <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  min_val <- pmin(dist_ratio, 1)
  cor_part <- (1 - (2 / pi * (min_val * sqrt(1 - min_val^2) + asin(min_val)))) * (dist_ratio <= 1)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_cubic <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- (1 - (7 / 1 * dist_ratio^2) + (35 / 4 * dist_ratio^3) - (7 / 2 * dist_ratio^5) + (3 / 4 * dist_ratio^7)) * (dist_ratio <= 1)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_pentaspherical <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- (1 - (15 / 8 * dist_ratio) + (5 / 4 * dist_ratio^3) - (3 / 8 * dist_ratio^5)) * (dist_ratio <= 1)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_wave <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- sin(dist_ratio) / (dist_ratio)
  cor_part[h == 0] <- 1
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_jbessel <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_product <- h * params[["range"]]
  cor_part <- besselJ(as.matrix(pmin(dist_product, 100000)), 0)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_gravity <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- (1 + dist_ratio^2)^(-1 / 2)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_rquad <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- (1 + dist_ratio^2)^(-1)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_magnetic <- function(params, dist_object, anisotropy, ...) {
  h <- get_euclid_mat_observed(params, dist_object, anisotropy)
  dist_ratio <- h / params[["range"]]
  cor_part <- (1 + dist_ratio^2)^(-3 / 2)
  params[["de"]] * cor_part
}

#' @export
cov_matrix.euclid_none <- function(params, dist_object, anisotropy, ...) {
  0
}

###############################################################################
# NUGGET COVARIANCES
###############################################################################

#' @export
cov_matrix.nugget_nugget <- function(params, dist_object, de_scale, diagtol = 0, ...) { # de scale for diagonal stability
  n_obs <- length(dist_object$network_index)
  nugget_matrix_observed <- Matrix::Matrix(0, nrow = n_obs, ncol = n_obs)
  diag(nugget_matrix_observed) <- max(params[["nugget"]], 1e-4 * de_scale, diagtol)
  nugget_matrix_observed
}

#' @export
cov_matrix.nugget_none <- function(params, dist_object, de_scale, diagtol = 0, ...) { # de scale for diagonal stability
  # CANT ZERO OUT HERE LIKE TAILUP, TAILDOWN, AND EUCLIDEAN BECAUSE WE NEED DIAGAONAL STABILITY
  n_obs <- length(dist_object$network_index)
  nugget_matrix_observed <- Matrix::Matrix(0, nrow = n_obs, ncol = n_obs)
  diag(nugget_matrix_observed) <- max(0, 1e-4 * de_scale, diagtol)
  nugget_matrix_observed
}

###############################################################################
# EUCILD MATRIX OBSERVED
###############################################################################

get_euclid_mat_observed <- function(params, dist_object, anisotropy) {
  # THIS IS ELEMENT BY ELEMENT

  # compute distance matrix with anisotropic coordinates or return
  # already stored euclidean distance matrix
  if (anisotropy) {
    new_coords_observed <- transform_anis(
      dist_object$.xcoord,
      dist_object$.ycoord,
      rotate = params[["rotate"]],
      scale = params[["scale"]]
    )
    euclid_mat_observed <- Matrix::Matrix(as.matrix(dist(cbind(new_coords_observed$xcoord_val, new_coords_observed$ycoord_val))), sparse = TRUE)
  } else {
    euclid_mat_observed <- dist_object$euclid_mat
  }
  euclid_mat_observed
}
