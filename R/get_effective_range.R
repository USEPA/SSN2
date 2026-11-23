#' Find the effective (practical) range of a spatial covariance function
#'
#' The "practical range" convention from geostatistics: the distance at
#' which the spatial correlation function drops to (and, for monotone
#' decaying functions, stays below) \code{target} which is 0.05 by default, the
#' standard convention. Compact and monotone functions have closed forms 
#' except for matern, which is numerically solved for. 
#'
#' Two families need special handling because their correlation functions
#' are not monotone decaying:
#' \itemize{
#'   \item Euclidean \code{"wave"} oscillates while decaying, with envelope
#'     \code{1 / (dist / range)}. The value returned
#'     (\code{range / target}) is the distance beyond which the envelope
#'     itself guarantees correlation stays under \code{target} as a
#'     conservative bound, not the first crossing (which happens much
#'     sooner, at \code{pi * range}). Stream wave also returns
#'     \code{range / target}, but this is only a nominal candidate, not a bound:
#'     its correlation decays at order \code{log(h / range) / (h / range)}.
#'     A warning is issued for either component.
#'   \item \code{"jbessel"}'s distance argument is \code{dist * range}, not
#'     \code{dist / range} the way every other family works and so effective
#'     range is \strong{inversely} related to \code{range} here (larger
#'     range means faster decay). The value returned uses the standard
#'     large-argument asymptotic envelope for the Bessel J0 function,
#'     \code{sqrt(2 / (pi * x))}. A warning is issued.
#' }
#'
#' @param params An object from [tailup_params()], [taildown_params()],
#'   [euclid_params()], or [nugget_params()].
#' @param target The correlation threshold defining the effective range. The
#'   default, 0.05, is the standard "practical range" convention. Must be
#'   strictly between zero and one.
#' @return One distance. Tail-up and tail-down results use stream-distance
#'   units; Euclidean results use coordinate-distance units. No spatial
#'   component (`"none"` or nugget) returns zero.
#' @details
#'
#' A single distance therefore does not describe every pair's correlation.
#' Tail-up and tail-down may have different effective ranges, even though
#' both use stream distance. Euclidean ranges are kept separate. With
#' anisotropy, the result uses transformed distance (the major-axis range);
#' the minor-axis range is `scale` times that distance.
#'
#' @noRd
#' @examples
#' get_effective_range(tailup_params("exponential", de = 1, range = 1000))
#' get_effective_range(taildown_params("mariah", de = 1, range = 1000))
#' get_effective_range(euclid_params("gaussian", de = 1, range = 200))
get_effective_range <- function(params, target = 0.05) {
  if (!is.numeric(target) || length(target) != 1L || !is.finite(target) ||
      target <= 0 || target >= 1) {
    stop("target must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  component <- sub("_.*$", "", class(params)[1])
  type <- remove_covtype(class(params)[1])
  if (!is.numeric(params) || !component %in% c("tailup", "taildown", "euclid", "nugget")) {
    stop("params must be an SSN covariance parameter object.", call. = FALSE)
  }
  switch(component,
    tailup = check_tailup_type(type), taildown = check_taildown_type(type),
    euclid = check_euclid_type(type), nugget = check_nugget_type(type)
  )
  if (type == "none" || component == "nugget") return(0)
  range <- unname(params[names(params) == "range"])
  if (length(range) != 1L || !is.finite(range) || range <= 0) {
    stop("The covariance range parameter must be finite and positive.", call. = FALSE)
  }
  if (component == "euclid" && euclid_has_extra(type)) {
    extra <- unname(params[names(params) == "extra"])
    if (length(extra) != 1L || !is.finite(extra) || extra <= 0 ||
        (type == "pexponential" && extra > 2)) {
      stop("The covariance extra parameter is invalid.", call. = FALSE)
    }
  }
  if (type %in% c("linear", "spherical", "epa", "circular", "cubic", "pentaspherical")) return(range)
  if (type == "exponential") return(-log(target) * range)
  if (type == "wave") {
    warning("\"wave\" effective range is not well defined. A potential candidate value is used.", call. = FALSE)
    return(range / target)
  }
  if (component == "euclid") {
    if (type == "gaussian") return(sqrt(-log(target)) * range)
    if (type %in% c("gravity", "rquad", "magnetic", "cauchy")) {
      p <- switch(type, gravity = 0.5, rquad = 1, magnetic = 1.5, cauchy = extra)
      return(range * sqrt(expm1(-log(target) / p)))
    }
    if (type == "pexponential") return((-log(target) * range)^(1 / extra))
    if (type == "jbessel") {
      warning("\"", type, "\" effective range is not well defined. A potential candidate value is used.", call. = FALSE)
      return(2 / (pi * target^2 * range))
    }
  }
  get_effective_range_monotone(params, target)
}

#' Numerically solve for the effective range of a monotone-decreasing kernel
#'
#' Used by \code{\link{get_effective_range}()} for kernels without a closed-form
#' inverse: evaluates the (unit-range) correlation at exponentially growing
#' distances to bracket a root, then solves for the distance at which
#' correlation equals \code{target} via \code{uniroot()}.
#'
#' @param params A single covariance component's parameter object, with
#'   \code{range} rescaled back onto the original scale afterward.
#' @param target The target correlation value (e.g. 0.05).
#'
#' @return The effective range: the distance at which correlation equals
#'   \code{target}.
#'
#' @noRd
get_effective_range_monotone <- function(params, target) {
  range <- params[["range"]]
  params[["range"]] <- params[["de"]] <- 1
  f <- function(distance) {
    distances <- list(hydro_mat = distance, a_mat = 0, b_mat = distance,
                      mask_mat = 1, w_mat = 1, euclid_mat = distance)
    as.numeric(cov_matrix(params, distances, anisotropy = FALSE)) - target
  }
  upper <- 1
  for (i in seq_len(1024L)) {
    value <- f(upper)
    if (!is.finite(value)) stop("The effective-range correlation could not be evaluated.", call. = FALSE)
    if (value <= 0) {
      return(range * uniroot(f, c(0, upper), tol = sqrt(.Machine$double.eps))$root)
    }
    upper <- 2 * upper
  }
  stop("Could not bracket the effective range.", call. = FALSE)
}
