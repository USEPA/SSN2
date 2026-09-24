#' Build variance-allocation-by-range starting-value combinations for active covariance components
#'
#' For the active spatial/nugget components, builds one equal-allocation row
#' plus one row per proper subset of components where that subset dominates
#' (90% of variance) and the rest share the remainder (10%). Each allocation
#' row is then crossed with corner (0.25/0.75-quantile-style) range
#' combinations for the active spatial components -- every corner when there
#' are at most two spatial components, or a reduced representative subset
#' (equal, all-spatial, and spatial-pair corners) when all three are active,
#' to avoid a full 8-way cross.
#'
#' @param active A length-4 logical vector (tailup, taildown, euclid,
#'   nugget) indicating which components are active (not type
#'   \code{"none"}).
#'
#' @return A data frame with columns \code{tailup_de}, \code{taildown_de},
#'   \code{euclid_de}, \code{nugget} (variance proportions) and
#'   \code{tailup_range}, \code{taildown_range}, \code{euclid_range}
#'   (proportions of a distance anchor, \code{Inf} for inactive components).
#'
#' @noRd
cov_initial_design <- function(active) {
  components <- which(active)
  spatial <- which(active[1:3])
  m <- length(components)
  subsets <- list(integer())
  if (m > 1L) {
    for (size in seq_len(m - 1L)) {
      subsets <- c(subsets, combn(components, size, simplify = FALSE))
    }
  }
  weights <- matrix(0, length(subsets), 4L)
  if (m) weights[1, components] <- 1 / m
  for (i in seq_along(subsets)[-1]) {
    selected <- subsets[[i]]
    weights[i, components] <- 0.1 / (m - length(selected))
    weights[i, selected] <- 0.9 / length(selected)
  }

  corners <- matrix(Inf, 1L, 3L)
  if (length(spatial)) {
    values <- as.matrix(expand.grid(rep(list(c(0.25, 0.75)), length(spatial))))
    corners <- matrix(Inf, nrow(values), 3L)
    corners[, spatial] <- values
  }
  if (length(spatial) < 3L) {
    allocation <- rep(seq_len(nrow(weights)), each = nrow(corners))
    ranges <- corners[rep(seq_len(nrow(corners)), nrow(weights)), , drop = FALSE]
  } else {
    allocation <- rep(seq_len(nrow(weights)), each = 2L)
    ranges <- corners[rep(c(1L, 8L), nrow(weights)), , drop = FALSE]
    # Spend mixed-range rows on equal, all-spatial, and spatial-pair allocations.
    for (i in seq_along(subsets)) {
      selected <- subsets[[i]]
      if (i == 1L || setequal(selected, spatial)) {
        extra <- corners[2:7, , drop = FALSE]
      } else if (length(selected) == 2L && all(selected %in% spatial)) {
        first <- rep(0.25, 3L)
        first[selected] <- 0.75
        extra <- rbind(first, 1 - first)
      } else {
        next
      }
      allocation <- c(allocation, rep(i, nrow(extra)))
      ranges <- rbind(ranges, extra)
    }
  }
  grid <- as.data.frame(cbind(weights[allocation, , drop = FALSE], ranges))
  names(grid) <- c("tailup_de", "taildown_de", "euclid_de", "nugget",
                   "tailup_range", "taildown_range", "euclid_range")
  grid
}

#' Build the full covariance-parameter starting-value candidate grid
#'
#' Scales \code{\link{cov_initial_design}()}'s variance-allocation/range
#' proportions by an OLS-residual-variance anchor and distance anchors
#' (\code{data_object$tail_max}/\code{euclid_max}), converts range
#' proportions into kernel-specific range parameters via
#' \code{\link{get_cov_initial_range}()}, expands over Euclidean shape
#' parameters (if the active Euclidean kernel has one), GLM dispersion
#' starting values (if \code{is_glm}), and anisotropy rotate/scale, adds
#' random-effect starting allocations via
#' \code{\link{add_randcov_grid}()}, and replaces any supplied/known
#' parameter values (via \code{cov_grid_replace}/\code{cov_grid_replace_glm})
#' before deduplicating.
#'
#' @param initial_NA_object A joint covariance initial-value object with
#'   \code{NA}-filled unknown fields.
#' @param data_object A model data object with \code{s2} (OLS residual
#'   variance), \code{tail_max}, \code{euclid_max} (distance anchors), and
#'   \code{anisotropy}.
#' @param is_glm Whether to also expand over GLM dispersion starting values
#'   and floor variance allocations away from zero.
#'
#' @return A data frame of candidate starting-value combinations, one row per
#'   candidate.
#'
#' @noRd
build_cov_initial_grid <- function(initial_NA_object, data_object, is_glm = FALSE) {
  types <- c("tailup", "taildown", "euclid", "nugget")
  active <- vapply(types, function(type) {
    !inherits(initial_NA_object[[paste0(type, "_initial")]], paste0(type, "_none"))
  }, logical(1))
  grid <- cov_initial_design(active)
  ns2 <- 1.2 * data_object$s2
  variance_names <- c("tailup_de", "taildown_de", "euclid_de", "nugget")
  grid[, variance_names] <- ns2 * grid[, variance_names]
  if (is_glm) {
    for (name in variance_names[active]) grid[[name]] <- pmax(grid[[name]], 0.05)
  }
  grid$tailup_range <- get_cov_initial_range(initial_NA_object$tailup_initial,
                                            grid$tailup_range * data_object$tail_max)
  grid$taildown_range <- get_cov_initial_range(initial_NA_object$taildown_initial,
                                              grid$taildown_range * data_object$tail_max)
  euclid_distance <- grid$euclid_range * data_object$euclid_max
  grid$euclid_range <- get_cov_initial_range(initial_NA_object$euclid_initial, euclid_distance)
  grid$rotate <- 0
  grid$scale <- 1

  euclid <- initial_NA_object$euclid_initial
  euclid_type <- remove_covtype(class(euclid)[1])
  if (euclid_has_extra(euclid_type)) {
    extras <- switch(euclid_type, matern = c(1, 4), cauchy = c(0.5, 2),
                     pexponential = c(0.4, 1.6))
    if (!is.na(euclid$initial[["extra"]])) extras <- euclid$initial[["extra"]]
    grid <- do.call(rbind, lapply(extras, function(extra) {
      candidate <- grid
      starts <- get_euclid_start_values(euclid, euclid_distance,
                                       data_object$euclid_max, extra = extra)
      candidate$euclid_range <- starts$range
      candidate$euclid_extra <- starts$extra
      candidate
    }))
  }
  if (is_glm) {
    grid$dispersion <- 1
    second <- grid
    second$dispersion <- 100
    grid <- rbind(grid, second)
  }
  if (data_object$anisotropy) {
    second <- grid
    second$rotate <- pi / 2
    second$scale <- 0.5
    grid <- rbind(grid, second)
  }
  grid <- add_randcov_grid(grid, initial_NA_object, data_object, ns2)
  replace <- if (is_glm) cov_grid_replace_glm else cov_grid_replace
  grid <- unique(replace(grid, initial_NA_object, data_object))
  rownames(grid) <- NULL
  grid
}

#' Convert a distance anchor into a kernel-specific range starting value
#'
#' Divides the given distance by a family-specific scaling factor (matching
#' spmodel's own automatic range-starting-value convention): 3 for
#' exponential/Matern/powered-exponential, \code{sqrt(3)} for Euclidean
#' Gaussian, the stream Gaussian kernel's own 0.05-correlation-distance
#' divisor for stream (tailup/taildown) Gaussian, and 1 otherwise.
#'
#' @param initial A single covariance component's initial-value object
#'   (used only for its class/type).
#' @param distance A distance value (or vector) to convert.
#'
#' @return \code{distance} divided by the family-specific scaling factor
#'   (\code{Inf} if \code{type} is \code{"none"}).
#'
#' @noRd
get_cov_initial_range <- function(initial, distance) {
  type <- remove_covtype(class(initial)[1])
  if (type == "none") return(rep(Inf, length(distance)))
  if (type == "gaussian" && !inherits(initial, "euclid_gaussian")) {
    # The stream Gaussian kernel differs from the Euclidean Gaussian kernel.
    params <- if (inherits(initial, "tailup_gaussian")) tailup_params else taildown_params
    divisor <- get_effective_range(params("gaussian", de = 1, range = 1))
  } else {
    divisor <- switch(type, exponential = 3, gaussian = sqrt(3),
                      matern = 3, pexponential = 3, 1)
  }
  distance / divisor
}
