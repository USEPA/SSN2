#' Perform anisotropy transformation for Euclidean distance
#'
#' @param xcoord_val x-coordinate
#' @param ycoord_val y-coordinate
#' @param rotate Anisotropy rotation parameter
#' @param scale Anisotropy scale parameter
#'
#' @noRd
transform_anis <- function(xcoord_val, ycoord_val, rotate, scale) {
  rotate_clockwise <- matrix(c(cos(rotate), sin(rotate), -sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
  scale_yaxis <- matrix(c(1, 0, 0, 1 / scale), nrow = 2, ncol = 2, byrow = TRUE)
  coords <- rbind(xcoord_val, ycoord_val)
  new_coords <- (scale_yaxis %*% rotate_clockwise) %*% coords
  list(xcoord_val = new_coords[1, ], ycoord_val = new_coords[2, ])
}

#' Perform inverse anisotropy transformation for Euclidean distance
#'
#' @param xcoord_val x-coordinate
#' @param ycoord_val y-coordinate
#' @param rotate Anisotropy rotation parameter
#' @param scale Anisotropy scale parameter
#'
#' @noRd
transform_anis_inv <- function(xcoord_val, ycoord_val, rotate, scale) {
  rotate_clockwise_inv <- matrix(c(cos(rotate), -sin(rotate), sin(rotate), cos(rotate)), nrow = 2, ncol = 2, byrow = TRUE)
  scale_yaxis_inv <- matrix(c(1, 0, 0, scale), nrow = 2, ncol = 2, byrow = TRUE)
  coords <- rbind(xcoord_val, ycoord_val)
  new_coords <- (rotate_clockwise_inv %*% scale_yaxis_inv) %*% coords
  list(xcoord_val = new_coords[1, ], ycoord_val = new_coords[2, ])
}
