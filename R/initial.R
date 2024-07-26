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
#' @details Create an initial object for use with [ssn_lm()] or [ssn_glm()].
#'   \code{NA} values can be given for \code{ie}, \code{rotate}, and \code{scale}, which lets
#'   these functions find initial values for parameters that are sometimes
#'   otherwise assumed known (e.g., \code{rotate} and \code{scale} with [ssn_lm()] and [ssn_glm()].
#'   Parametric forms for each spatial covariance type are presented below.
#'
#'   \code{tailup_type} Details: Let \eqn{D} be a matrix of hydrologic distances,
#'   \eqn{W} be a diagonal matrix of weights from \code{additive}, \eqn{r = D / range},
#'   and \eqn{I} be
#'   an identity matrix. Then parametric forms for flow-connected
#'   elements of \eqn{R(zu)} are given below:
#'   \itemize{
#'     \item linear: \eqn{(1 - r) * (r <= 1) * W}
#'     \item spherical: \eqn{(1 - 1.5r + 0.5r^3) * (r <= 1) * W}
#'     \item exponential: \eqn{exp(-r) * W}
#'     \item mariah: \eqn{log(90r + 1) / 90r * (D > 0) + 1 * (D = 0) * W}
#'     \item epa: \eqn{(D - range)^2 * F * (r <= 1) * W / 16range^5}
#'     \item none: \eqn{I} * W
#'   }
#'
#'   Details describing the \code{F} matrix in the \code{epa} covariance are given in Garreta et al. (2010).
#'   Flow-unconnected elements of \eqn{R(zu)} are assumed uncorrelated.
#'   Observations on different networks are also assumed uncorrelated.
#'
#'   \code{taildown_type} Details: Let \eqn{D} be a matrix of hydrologic distances,
#'   \eqn{r = D / range},
#'   and \eqn{I} be an identity matrix. Then parametric forms for flow-connected
#'   elements of \eqn{R(zd)} are given below:
#'   \itemize{
#'     \item linear: \eqn{(1 - r) * (r <= 1)}
#'     \item spherical: \eqn{(1 - 1.5r + 0.5r^3) * (r <= 1)}
#'     \item exponential: \eqn{exp(-r)}
#'     \item mariah: \eqn{log(90r + 1) / 90r * (D > 0) + 1 * (D = 0)}
#'     \item epa: \eqn{(D - range)^2 * F1 * (r <= 1) / 16range^5}
#'     \item none: \eqn{I}
#'   }
#'
#'   Now let \eqn{A} be a matrix that contains the shorter of the two distances
#'   between two sites and the common downstream junction, \eqn{r1 = A / range},
#'   \eqn{B} be a matrix that contains the longer of the two distances between two sites and the
#'   common downstream junction, \eqn{r2 = B / range},  and \eqn{I} be an identity matrix.
#'   Then parametric forms for flow-unconnected elements of \eqn{R(zd)} are given below:
#'   \itemize{
#'     \item linear: \eqn{(1 - r2) * (r2 <= 1)}
#'     \item spherical: \eqn{(1 - 1.5r1 + 0.5r2) * (1 - r2)^2 * (r2 <= 1)}
#'     \item exponential: \eqn{exp(-(r1 + r2))}
#'     \item mariah: \eqn{(log(90r1 + 1) - log(90r2 + 1)) / (90r1 - 90r2) * (A =/ B) + (1 / (90r1 + 1)) * (A = B)}
#'     \item epa: \eqn{(B - range)^2 * F2 * (r2 <= 1) / 16range^5}
#'     \item none: \eqn{I}
#'   }
#'
#'   Details describing the \code{F1} and \code{F2} matrices in the \code{epa}
#'   covariance are given in Garreta et al. (2010).
#'   Observations on different networks are assumed uncorrelated.
#'
#'  \code{euclid_type} Details: Let \eqn{D} be a matrix of Euclidean distances,
#'  \eqn{r = D / range}, and \eqn{I} be an identity matrix. Then parametric
#'  forms for elements of \eqn{R(ze)} are given below:
#'   \itemize{
#'     \item exponential: \eqn{exp(- r )}
#'     \item spherical: \eqn{(1 - 1.5r + 0.5r^3) * (r <= 1)}
#'     \item gaussian: \eqn{exp(- r^2 )}
#'     \item cubic: \eqn{(1 - 7r^2 + 8.75r^3 - 3.5r^5 + 0.75r^7) * (r <= 1)}
#'     \item pentaspherical: \eqn{(1 - 1.875r + 1.25r^3 - 0.375r^5) * (r <= 1)}
#'     \item cosine: \eqn{cos(r)}
#'     \item wave: \eqn{sin(r) * (h > 0) / r + (h = 0)}
#'     \item jbessel: \eqn{Bj(h * range)}, Bj is Bessel-J function
#'     \item gravity: \eqn{(1 + r^2)^{-0.5}}
#'     \item rquad: \eqn{(1 + r^2)^{-1}}
#'     \item magnetic: \eqn{(1 + r^2)^{-1.5}}
#'     \item none: \eqn{I}
#'   }
#'
#'   \code{nugget_type} Details: Let \eqn{I} be an identity matrix and \eqn{0}
#'    be the zero matrix. Then parametric
#'    forms for elements the nugget variance are given below:
#'   \itemize{
#'     \item nugget: \eqn{I}
#'     \item none: \eqn{0}
#'   }
#'   In short, the nugget effect is modeled when \code{nugget_type} is \code{"nugget"}
#'   and omitted when \code{nugget_type} is \code{"none"}.
#'
#'   Dispersion and random effect initial objects are specified via
#'   [spmodel::dispersion_initial()] and [spmodel::randcov_initial()], respectively.
#'
#' @return A list with two elements: \code{initial} and \code{is_known}.
#'   \code{initial} is a named numeric vector indicating the spatial covariance parameters
#'   with specified initial and/or known values. \code{is_known} is a named
#'   numeric vector indicating whether the spatial covariance parameters in
#'   \code{initial} are known or not. The class of the list
#'   matches the the relevant spatial covariance type.
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

#' A helper to create an initial value object
#'
#' @param tailup_type Tailup covariance type
#' @param taildown_type Taildown covariance type
#' @param euclid_type Euclidean covariance type
#' @param nugget_type Nugget covariance type
#' @param tailup_initial Tailup initial object (if specified)
#' @param taildown_initial Taildown initial object (if specified)
#' @param euclid_initial Euclidean initial object (if specified)
#' @param nugget_initial Nugget initial object (if specified)
#'
#' @noRd
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

#' A helper to create an initial value glm object
#'
#' @param tailup_type Tailup covariance type
#' @param taildown_type Taildown covariance type
#' @param euclid_type Euclidean covariance type
#' @param nugget_type Nugget covariance type
#' @param tailup_initial Tailup initial object (if specified)
#' @param taildown_initial Taildown initial object (if specified)
#' @param euclid_initial Euclidean initial object (if specified)
#' @param nugget_initial Nugget initial object (if specified)
#' @param family The generalized linear model family
#' @param dispersion_inital Dispersion initial object
#'
#' @noRd
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

#' Create a vector that specifies whether covariance parameters are to be assumed known
#'
#' @param params  A covariance parameter vector
#' @param known Which covariance parameters assumed known
#'
#' @noRd
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
