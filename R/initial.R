#' Create a covariance parameter initial object
#'
#' @description Create a covariance parameter initial object that specifies
#'   initial and/or known values to use while estimating specific covariance parameters
#'   with [ssn_lm()] or [ssn_glm()].  See [spmodel::randcov_initial()] for documentation regarding
#'   random effect covariance parameter initial objects.
#'
#' @param tailup_type The tailup covariance function type. Available options
#'   include \code{"linear"}, \code{"spherical"}, \code{"exponential"},
#'   \code{"mariah"}, \code{"epa"}, and \code{"none"}.
#' @param taildown_type The taildown covariance function type. Available options
#'   include \code{"linear"}, \code{"spherical"}, \code{"exponential"},
#'   \code{"mariah"}, \code{"epa"}, and \code{"none"}.
#' @param euclid_type The euclidean covariance function type. Available options
#'   include \code{"spherical"}, \code{"exponential"}, \code{"gaussian"},
#'   \code{"cosine"}, \code{"cubic"}, \code{"pentaspherical"}, \code{"wave"},
#'    \code{"jbessel"}, \code{"gravity"}, \code{"rquad"}, \code{"magnetic"}, and
#'    \code{"none"}.
#' @param nugget_type The nugget covariance function type. Available options
#'   include \code{"nugget"} or \code{"none"}.
#' @param de The spatially dependent (correlated) random error variance. Commonly referred to as
#'   a partial sill.
#' @param range The correlation parameter.
#' @param rotate Anisotropy rotation parameter (from 0 to \eqn{\pi} radians) for
#'   the euclidean portion of the covariance. A value of 0 (the default) implies no rotation.
#' @param scale Anisotropy scale parameter (from 0 to 1) for
#'   the euclidean portion of the covariance. A value of 1 (the default) implies no scaling.
#' @param nugget The spatially independent (not correlated) random error variance. Commonly referred to as
#'   a nugget.
#' @param known A character vector indicating which covariance parameters are to be
#'   assumed known. The value \code{"given"} is shorthand for assuming all
#'   covariance parameters given to \code{*_initial()} are assumed known.
#'
#' @name ssn_initial
#'
#' @export
#'
#'
#' @examples
#' tailup_initial("exponential", de = 1, range = 20, known = "range")
#' tailup_initial("exponential", de = 1, range = 20, known = "given")
#' euclid_initial("spherical", de = 2, range = 4, scale = 0.8, known = c("range", "scale"))
#' dispersion_initial("nbinomial", dispersion = 5)
#'
#' @references
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach
#' to geostatistical modeling in stream networks. \emph{Ecology} \bold{91(3)},
#' 644--651.
#'
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial
#' statistical models of stream networks (with discussion).
#' \emph{Journal of the American Statistical Association} \bold{105}, 6--18.
#' DOI: 10.1198/jasa.2009.ap08248.  Rejoinder pgs. 22--24.
tailup_initial <- function(tailup_type, de, range, known) {
  check_tailup_type(tailup_type)

  # set defaults
  if (missing(de)) de <- NULL
  if (missing(range)) range <- NULL

  tailup_params_given <- c(de = unname(de), range = unname(range))

  is_known <- get_is_known(tailup_params_given, known)

  new_tailup_initial <- structure(list(initial = tailup_params_given, is_known = is_known),
    class = paste("tailup", tailup_type, sep = "_")
  )
  new_tailup_initial
}

#' @rdname ssn_initial
#' @export
taildown_initial <- function(taildown_type, de, range, known) {
  check_taildown_type(taildown_type)

  # set defaults
  if (missing(de)) de <- NULL
  if (missing(range)) range <- NULL

  taildown_params_given <- c(de = unname(de), range = unname(range))

  is_known <- get_is_known(taildown_params_given, known)

  new_taildown_initial <- structure(list(initial = taildown_params_given, is_known = is_known),
    class = paste("taildown", taildown_type, sep = "_")
  )
  new_taildown_initial
}

#' @rdname ssn_initial
#' @export
euclid_initial <- function(euclid_type, de, range, rotate, scale, known) {
  check_euclid_type(euclid_type)

  # set defaults
  if (missing(de)) de <- NULL
  if (missing(range)) range <- NULL
  if (missing(rotate)) rotate <- NULL
  if (missing(scale)) scale <- NULL

  euclid_params_given <- c(
    de = unname(de), range = unname(range),
    rotate = unname(rotate), scale = unname(scale)
  )

  is_known <- get_is_known(euclid_params_given, known)

  new_euclid_initial <- structure(list(initial = euclid_params_given, is_known = is_known),
    class = paste("euclid", euclid_type, sep = "_")
  )
  new_euclid_initial
}

#' @rdname ssn_initial
#' @export
nugget_initial <- function(nugget_type, nugget, known) {
  check_nugget_type(nugget_type)

  if (missing(nugget)) nugget <- NULL
  nugget_params_given <- c(nugget = unname(nugget))

  is_known <- get_is_known(nugget_params_given, known)

  new_nugget_initial <- structure(list(initial = nugget_params_given, is_known = is_known),
    class = paste("nugget", nugget_type, sep = "_")
  )
  new_nugget_initial
}

get_initial_object <- function(tailup_type, taildown_type, euclid_type, nugget_type,
                               tailup_initial, taildown_initial, euclid_initial, nugget_initial) {
  if (is.null(tailup_initial)) tailup_initial <- tailup_initial(tailup_type)
  if (is.null(taildown_initial)) taildown_initial <- taildown_initial(taildown_type)
  if (is.null(euclid_initial)) euclid_initial <- euclid_initial(euclid_type)
  if (is.null(nugget_initial)) nugget_initial <- nugget_initial(nugget_type)

  if (all(c(tailup_type, taildown_type, euclid_type, nugget_type) == "none")) {
    stop("At least one covariance type must be different from \"none\".", call. = FALSE)
  }

  initial_object <- list(
    tailup_initial = tailup_initial,
    taildown_initial = taildown_initial,
    euclid_initial = euclid_initial,
    nugget_initial = nugget_initial
  )

  initial_object
}

get_initial_object_glm <- function(tailup_type, taildown_type, euclid_type, nugget_type,
                                   tailup_initial, taildown_initial, euclid_initial, nugget_initial,
                                   family, dispersion_initial) {
  if (is.null(tailup_initial)) tailup_initial <- tailup_initial(tailup_type)
  if (is.null(taildown_initial)) taildown_initial <- taildown_initial(taildown_type)
  if (is.null(euclid_initial)) euclid_initial <- euclid_initial(euclid_type)
  if (is.null(nugget_initial)) nugget_initial <- nugget_initial(nugget_type)
  if (is.null(dispersion_initial)) dispersion_initial <- dispersion_initial(eval(family))

  initial_object <- list(
    tailup_initial = tailup_initial,
    taildown_initial = taildown_initial,
    euclid_initial = euclid_initial,
    nugget_initial = nugget_initial,
    dispersion_initial = dispersion_initial
  )

  initial_object
}

get_is_known <- function(params, known) {
  if (missing(known)) {
    is_known <- rep(FALSE, length(params))
  } else {
    if (identical(known, "given")) {
      is_known <- rep(TRUE, length(params))
    } else {
      is_known <- names(params) %in% known
    }
  }
  names(is_known) <- names(params)

  # error if NA and known
  params_NA <- which(is.na(params))
  if (any(is_known[params_NA])) {
    stop("params initial values cannot be NA and known.", call. = FALSE)
  }

  is_known
}
