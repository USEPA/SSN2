#' Compute an empirical stream-network diagnostic
#'
#' @description Compute the empirical semivariogram or autocovariance for
#'   varying bin sizes and cutoff values.
#'
#' @param formula A formula describing the fixed effect structure.
#' @param ssn.object A spatial stream network object with class \code{SSN}.
#' @param type The Torgegram type. A vector with possible values \code{"flowcon"}
#'   for flow-connected distances, \code{"flowuncon"} for flow-unconnected distances,
#'   and \code{"euclid"} for Euclidean distances.
#'   The default is to show both flow-connected and
#'   flow-unconnected distances.
#' @param cloud A logical indicating whether the selected diagnostic should be
#'   summarized by distance class or not. When \code{cloud = FALSE} (the
#'   default), pairwise values are binned and averaged within distance classes.
#'   When \code{cloud} = TRUE, all selected pairwise values and distances are
#'   returned.
#' @param robust A logical indicating whether the robust semivariogram
#' (Cressie and Hawkins, 1980) is used for each \code{type}. The default is \code{FALSE}.
#' @param eacf A logical indicating whether to compute the empirical
#'   autocovariance rather than the empirical semivariogram. The default is
#'   \code{FALSE}. The same pairs, distances, bins, cutoff, and partition
#'   restrictions are used in either case. When \code{TRUE},
#'   \code{robust = TRUE} is not available.
#' @param bins The number of equally spaced bins. The default is 15. Ignored if
#'   \code{cloud = TRUE}.
#' @param cutoff The maximum distance considered. The default is half the
#'   maximum observed distance for each selected pair type.
#' @param partition_factor An optional formula specifying the partition factor.
#'   If specified, diagnostic values are only computed for observations sharing
#'   the same level of the partition factor.
#'
#' @details With \code{eacf = FALSE}, the Torgegram is an empirical semivariogram
#'   used to visualize and model
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
#'   \eqn{1/(2N(h)) \sum (r1 - r2)^2}, where \eqn{N(h)} is the number of (unique)
#'   pairs in the set of observations whose distance separation is \code{h} and
#'   \code{r1} and \code{r2} are residuals corresponding to observations whose
#'   distance separation is \code{h}. The robust version is described by
#'   Cressie and Hawkins (1980).
#'
#'   With \code{eacf = TRUE}, the Torgegram is an empirical autocovariance used to visualize and model
#'   spatial dependence by estimating the autocovariance of a process at varying distances.
#'   For a constant-mean process, the
#'   autocovariance at distance \eqn{h} is denoted \eqn{Cov(h)} and defined as
#'   \eqn{Cov(z1, z2)}. Under second-order stationarity,
#'   \eqn{Cov(h) = Cov(0) - \gamma(h)}, where \eqn{\gamma(h)} is the semivariance function at distance \code{h}. Typically the residuals from an ordinary
#'   least squares fit defined by \code{formula} are second-order stationary with
#'   mean zero. These residuals are used to compute the empirical autocovariance.
#'   At a distance \code{h}, the empirical autocovariance is
#'   \eqn{1/N(h) \sum (r1 \times r2)}, where \eqn{N(h)} is the number of (unique)
#'   pairs in the set of observations whose distance separation is \code{h} and
#'   \code{r1} and \code{r2} are residuals corresponding to observations whose
#'   distance separation is \code{h}.
#'
#'   In \code{SSN2}, the Torgegram distance bins actually
#'   contain observations whose distance separation is \code{h +- c},
#'   where \code{c} is a constant determined implicitly by \code{bins}. Typically,
#'   only observations whose distance separation is below some cutoff are used
#'   to compute either diagnostic (this cutoff is determined by \code{cutoff}).
#'
#' @return With \code{eacf = FALSE}, a list with elements corresponding to \code{type}. Each element
#'   is data frame with distance bins (\code{bins}), the  average distance
#'   (\code{dist}), the semivariance (\code{gamma}), and the
#'   number of (unique) pairs (\code{np}) for the respective \code{type}.
#'   With \code{eacf = TRUE}, the empirical autocovariance (\code{acov})
#'   replaces \code{gamma}; the other columns are unchanged. When
#'   \code{cloud = TRUE}, each element contains \code{dist} and either
#'   \code{gamma} or \code{acov}, with one row for each unique pair.
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
#' Torgegram(Summer_mn ~ 1, mf04p, eacf = TRUE)
#' @references
#' Cressie, N & Hawkins, D.M. 1980. Robust estimation of the variogram.
#'   \emph{Journal of the International Association for Mathematical Geology},
#'   \strong{12}, 115-125.
#' Zimmerman, D. L., & Ver Hoef, J. M. (2017). The Torgegram for fluvial
#'   variography: characterizing spatial dependence on stream networks.
#'   \emph{Journal of Computational and Graphical Statistics},
#'   \bold{26(2)}, 253--264.
Torgegram <- function(formula, ssn.object,
                      type = c("flowcon", "flowuncon"), cloud = FALSE, robust = FALSE,
                      bins = 15, cutoff, partition_factor, eacf = FALSE) {
  if (!is.logical(eacf) || length(eacf) != 1 || is.na(eacf)) {
    stop("eacf must be a single non-missing logical value.", call. = FALSE)
  }

  if (eacf) {
    if (!is.logical(robust) || length(robust) != 1 || is.na(robust)) {
      stop("robust must be a single non-missing logical value.", call. = FALSE)
    }
    if (robust) {
      stop("robust = TRUE is not available when eacf = TRUE.", call. = FALSE)
    }
  }

  Torgegram_initial_object <- get_Torgegram_initial_object(type)
  # find distance object
  dist_object <- get_dist_object(ssn.object, Torgegram_initial_object,
    additive = NULL, anisotropy = FALSE
  )

  # find residuals
  lmod <- lm(formula = formula, data = ssn.object$obs)
  residuals <- residuals(lmod)

  warn_torgegram_zero_distance(ssn.object$obs)

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

  if (eacf) {
    residual_products <- tcrossprod(as.numeric(residuals))
    residual2_vector <- residual_products[upper.tri(residual_products)]
  } else {
    residual_mat_sqrt <- as.matrix(dist(residuals))
    residual_vector <- residual_mat_sqrt[upper.tri(residual_mat_sqrt)]
    residual2_vector <- residual_vector^2
  }

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

  if (cloud) {

    if ("flowcon" %in% type) {
      esv_list$flowcon <- get_esv_cloud(residual2_vector, flowcon_vector, cutoff, eacf)
    }

    if ("flowuncon" %in% type) {
      esv_list$flowuncon <- get_esv_cloud(residual2_vector, flowuncon_vector, cutoff, eacf)
    }

    if ("euclid" %in% type) {
      esv_list$euclid <- get_esv_cloud(residual2_vector, euclid_vector, cutoff, eacf)
    }

  } else {
    if (robust) {

      residual12_vector <- sqrt(residual_vector)
      if ("flowcon" %in% type) {
        esv_list$flowcon <- get_esv_robust(residual12_vector, flowcon_vector, bins, cutoff)
      }

      if ("flowuncon" %in% type) {
        esv_list$flowuncon <- get_esv_robust(residual12_vector, flowuncon_vector, bins, cutoff)
      }

      if ("euclid" %in% type) {
        esv_list$euclid <- get_esv_robust(residual12_vector, euclid_vector, bins, cutoff)
      }
    } else {
      if ("flowcon" %in% type) {
        esv_list$flowcon <- get_esv(residual2_vector, flowcon_vector, bins, cutoff, eacf)
      }

      if ("flowuncon" %in% type) {
        esv_list$flowuncon <- get_esv(residual2_vector, flowuncon_vector, bins, cutoff, eacf)
      }

      if ("euclid" %in% type) {
        esv_list$euclid <- get_esv(residual2_vector, euclid_vector, bins, cutoff, eacf)
      }
    }
  }

  new_esv_list <- structure(esv_list, class = "Torgegram", call = match.call(), cloud = cloud)
  if (eacf) attr(new_esv_list, "eacf") <- TRUE
  new_esv_list
}

#' Warn once when a Torgegram will drop zero-distance (coincident-site) pairs
#'
#' Matches spmodel's \code{esv()} zero-distance warning. A zero hydrologic
#' distance can mean either a genuine coincident-site pair or a masked-out
#' cross-network pair (hydrologic distance is undefined between networks, so
#' \code{hydro_mat} is forced to zero there too), so checking a hydrologic or
#' Euclidean distance vector directly for zeros is ambiguous. Physical
#' coincidence -- the actual condition being warned about -- is unambiguous in
#' the raw observation coordinates directly, regardless of which Torgegram
#' type is being computed, so it's checked there instead.
#'
#' @param obs An sf object of observations (\code{ssn.object$obs}).
#'
#' @return \code{NULL}, invisibly.
#'
#' @noRd
warn_torgegram_zero_distance <- function(obs) {
  if (anyDuplicated(sf::st_coordinates(obs)) > 0) {
    warning("Zero distances observed between at least one pair. Ignoring pairs.", call. = FALSE)
  }
  invisible(NULL)
}

get_esv <- function(resid2_vector, dist_vector, bins, cutoff, eacf = FALSE) {
  if (is.null(cutoff)) {
    cutoff <- max(dist_vector) * 0.5
  }
  index <- dist_vector > 0 & dist_vector <= cutoff
  dist_vector <- dist_vector[index]
  resid2 <- resid2_vector[index]

  dist_classes <- cut(dist_vector, breaks = seq(0, cutoff, length.out = bins + 1))

  gamma <- tapply(resid2, dist_classes, function(x) if (eacf) mean(x) else mean(x) / 2)

  # compute pairs within each class
  np <- tapply(resid2, dist_classes, length)

  # set as zero if necessary
  np <- ifelse(is.na(np), 0, np)

  # compute average distance within each class
  dist <- tapply(dist_vector, dist_classes, mean)

  # return output
  esv_out <- data.frame(bins = factor(levels(dist_classes), levels = levels(dist_classes)), dist, gamma, np)
  if (eacf) names(esv_out)[names(esv_out) == "gamma"] <- "acov"

  # set row names to NULL
  row.names(esv_out) <- NULL

  # return esv
  esv_out
}

get_esv_robust <- function(resid12_vector, dist_vector, bins, cutoff, formula) {
  # Cressie's robust estimator: averages sqrt(|differences|) instead of squared
  # differences (resid12_vector is already on that scale), then raises back
  # to the 4th power with a bias correction -- less sensitive to outlier pairs
  # than the classical estimator above
  if (is.null(cutoff)) {
    cutoff <- max(dist_vector) * 0.5
  }
  index <- dist_vector > 0 & dist_vector <= cutoff
  dist_vector <- dist_vector[index]
  resid12 <- resid12_vector[index]

  # compute semivariogram classes
  dist_classes <- cut(dist_vector, breaks = seq(0, cutoff, length.out = bins + 1))

  # compute squared differences within each class
  gamma <- tapply(resid12, dist_classes, function(x) {
    1 / (0.914 + (0.988 / length(x))) * (mean(x)^4)
  })

  # compute pairs within each class
  np <- tapply(resid12, dist_classes, length)

  # set as zero if necessary
  np <- ifelse(is.na(np), 0, np)

  # compute average distance within each class
  dist <- tapply(dist_vector, dist_classes, mean)

  # return output
  esv_out <- tibble::tibble(
    bins = factor(levels(dist_classes), levels = levels(dist_classes)),
    dist = as.numeric(dist),
    gamma = as.numeric(gamma),
    np = as.numeric(np)
  )

  # set row names to NULL
  # row.names(esv_out) <- NULL

  esv_out
}

get_esv_cloud <- function(residual2_vector, dist_vector, cutoff, eacf = FALSE) {
  # no binning/averaging -- every pair is returned as its own row for
  # plotting, but still restricted to pairs within cutoff, matching the
  # binned and robust paths (get_esv()/get_esv_robust())
  if (is.null(cutoff)) {
    cutoff <- max(dist_vector) * 0.5
  }
  index <- dist_vector > 0 & dist_vector <= cutoff
  dist_vector <- dist_vector[index]
  resid2 <- residual2_vector[index]

  esv_out <- tibble::tibble(dist = dist_vector, gamma = if (eacf) resid2 else resid2 / 2)
  if (eacf) names(esv_out)[names(esv_out) == "gamma"] <- "acov"

  # set row names to NULL
  # row.names(esv_out) <- NULL

  esv_out
}
