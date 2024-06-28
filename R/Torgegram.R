#' Compute the empirical semivariogram
#'
#' @description Compute the empirical semivariogram for varying bin sizes and
#'   cutoff values.
#'
#' @param formula A formula describing the fixed effect structure.
#' @param ssn.object A spatial stream network object with class \code{SSN}.
#' @param type The Torgegram type. A vector with possible values \code{"flowcon"}
#'   for flow-connected distances, \code{"flowuncon"} for flow-unconnected distances,
#'   and \code{"euclid"} for Euclidean distances. The default is to show both
#'   flow-connected and flow-unconnected distances.
#' @param bins The number of equally spaced bins. The default is 15.
#' @param cutoff The maximum distance considered.
#'   The default is half the diagonal of the bounding box from the coordinates.
#' @param partition_factor An optional formula specifying the partition factor.
#'   If specified, semivariances are only computed for observations sharing the
#'   same level of the partition factor.
#'
#' @details The Torgegram is an empirical semivariogram is a tool used to visualize and model
#'   spatial dependence by estimating the semivariance of a process at varying distances
#'   separately for flow-connected, flow-unconnected, and Euclidean distances.
#'   For a constant-mean process, the
#'   semivariance at distance \eqn{h} is denoted \eqn{\gamma(h)} and defined as
#'   \eqn{0.5 * Var(z1  - z2)}. Under second-order stationarity,
#'   \eqn{\gamma(h) = Cov(0) - Cov(h)}, where \eqn{Cov(h)} is the covariance function
#'   at distance \code{h}. Typically the residuals from an ordinary
#'   least squares fit defined by \code{formula} are second-order stationary with
#'   mean zero. These residuals are used to compute the empirical semivariogram.
#'   At a distance \code{h}, the empirical semivariance is
#'   \eqn{1/N(h) \sum (r1 - r2)^2}, where \eqn{N(h)} is the number of (unique)
#'   pairs in the set of observations whose distance separation is \code{h} and
#'   \code{r1} and \code{r2} are residuals corresponding to observations whose
#'   distance separation is \code{h}. In spmodel, these distance bins actually
#'   contain observations whose distance separation is \code{h +- c},
#'   where \code{c} is a constant determined implicitly by \code{bins}. Typically,
#'   only observations whose distance separation is below some cutoff are used
#'   to compute the empirical semivariogram (this cutoff is determined by \code{cutoff}).
#'
#' @return A list with elements correspond to \code{type}. Each element
#'   is data frame with distance bins (\code{bins}), the  average distance
#'   (\code{dist}), the semivariance (\code{gamma}), and the
#'   number of (unique) pairs (\code{np}) for the respective \code{type}.
#'
#' @export
#'
#' @seealso [plot.Torgegram()]
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' tg <- Torgegram(Summer_mn ~ 1, mf04p)
#' plot(tg)
#' @references
#' Zimmerman, D. L., & Ver Hoef, J. M. (2017). The Torgegram for fluvial
#'   variography: characterizing spatial dependence on stream networks.
#'   \emph{Journal of Computational and Graphical Statistics},
#'   \bold{26(2)}, 253--264.
Torgegram <- function(formula, ssn.object,
                      type = c("flowcon", "flowuncon"),
                      bins = 15, cutoff, partition_factor) {
  Torgegram_initial_object <- get_Torgegram_initial_object(type)
  # find distance object
  dist_object <- get_dist_object(ssn.object, Torgegram_initial_object,
    additive = NULL, anisotropy = FALSE
  )

  # find residuals
  lmod <- lm(formula = formula, data = ssn.object$obs)
  residuals <- residuals(lmod)
  residual_mat_sqrt <- as.matrix(dist(residuals))

  # find relevant vectors
  if ("flowcon" %in% type || "flowuncon" %in% type) {
    hydro_mat_mask <- dist_object$hydro_mat * dist_object$mask_mat
    hydro_vector <- as.matrix(hydro_mat_mask)[upper.tri(hydro_mat_mask)] # mat "maybe inefficient" warning
    b_vector <- as.matrix(dist_object$b_mat)[upper.tri(dist_object$b_mat)]
    flowcon_index <- b_vector == 0
    flowcon_vector <- hydro_vector * flowcon_index
    flowuncon_vector <- hydro_vector * !flowcon_index
  }

  if ("euclid" %in% type) {
    euclid_vector <- as.matrix(dist_object$euclid_mat)[upper.tri(dist_object$euclid_mat)]
  }

  residual2_vector <- residual_mat_sqrt[upper.tri(residual_mat_sqrt)]^2

  # handle partition factor
  if (!missing(partition_factor) && !is.null(partition_factor)) {
    # partition_mat_val <- triu(partition_matrix(partition_factor, data = data), k = 1)
    partition_mat_val <- as.matrix(partition_matrix(partition_factor, data = ssn.object$obs))
    partition_vector_val <- partition_mat_val[upper.tri(partition_mat_val)]
    partition_index <- partition_vector_val == 1

    if ("flowcon" %in% type || "flowuncon" %in% type) {
      flowcon_vector <- flowcon_vector[partition_index]
      flowuncon_vector <- flowuncon_vector[partition_index]
    }
    if ("euclid" %in% type) {
      euclid_vector <- euclid_vector[partition_index]
    }
    residual2_vector <- residual2_vector[partition_index]
  }

  # find cutoffs
  if (missing(cutoff) || is.null(cutoff)) {
    cutoff <- NULL
  }

  if (missing(type) || is.null(type)) {
    type <- c("flowcon", "flowuncon")
  }

  esv_list <- list()

  if ("flowcon" %in% type) {
    esv_list$flowcon <- get_esv(flowcon_vector, residual2_vector, bins, cutoff)
  }

  if ("flowuncon" %in% type) {
    esv_list$flowuncon <- get_esv(flowuncon_vector, residual2_vector, bins, cutoff)
  }


  if ("euclid" %in% type) {
    esv_list$euclid <- get_esv(euclid_vector, residual2_vector, bins, cutoff)
  }

  new_esv_list <- structure(esv_list, class = "Torgegram")
  new_esv_list
}

#' Get a single empirical semivariogram for flow-connected, flow-unconnected, or Euclidean distance
#'
#' @param dist_vector Distance vector
#' @param resid2_vector Residual squared vector
#' @param bins Distance bins
#' @param cutoff Distance cutoff
#'
#' @noRd
get_esv <- function(dist_vector, resid2_vector, bins, cutoff) {
  if (is.null(cutoff)) {
    cutoff <- max(dist_vector) * 0.5
  }
  index <- dist_vector > 0 & dist_vector <= cutoff
  dist_vector <- dist_vector[index]
  resid2 <- resid2_vector[index]

  dist_classes <- cut(dist_vector, breaks = seq(0, cutoff, length.out = bins + 1))

  # compute squared differences within each class
  gamma <- tapply(resid2, dist_classes, function(x) mean(x) / 2)

  # compute pairs within each class
  np <- tapply(resid2, dist_classes, length)

  # set as zero if necessary
  np <- ifelse(is.na(np), 0, np)

  # compute average distance within each class
  dist <- tapply(dist_vector, dist_classes, mean)

  # return output
  esv_out <- data.frame(bins = factor(levels(dist_classes), levels = levels(dist_classes)), dist, gamma, np)

  # set row names to NULL
  row.names(esv_out) <- NULL

  # return esv
  esv_out
}
