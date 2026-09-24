#' Create covariance parameter objects
#'
#' @description Create a covariance parameter object for use with other functions.
#'   See [spmodel::randcov_params()] for documentation regarding
#'   random effect covariance parameter objects.
#'
#' @param tailup_type The tailup covariance function type. Available options
#'   include \code{"linear"}, \code{"spherical"}, \code{"exponential"},
#'   \code{"mariah"}, \code{"epa"}, and \code{"none"}.
#' @param taildown_type The taildown covariance function type. Available options
#'   include \code{"linear"}, \code{"spherical"}, \code{"exponential"},
#'   \code{"mariah"}, \code{"epa"}, and \code{"none"}.
#' @param euclid_type The Euclidean covariance function type. Available options
#'   include \code{"spherical"}, \code{"exponential"}, \code{"gaussian"},
#'   \code{"circular"}, \code{"cubic"}, \code{"pentaspherical"}, \code{"wave"},
#'    \code{"jbessel"}, \code{"gravity"}, \code{"rquad"}, \code{"magnetic"},
#'    \code{"matern"}, \code{"cauchy"}, \code{"pexponential"}, and \code{"none"}.
#' @param nugget_type The nugget covariance function type. Available options
#'   include \code{"nugget"} or \code{"none"}.
#' @param de The spatially dependent (correlated) random error variance. Commonly referred to as
#'   a partial sill.
#' @param range The correlation parameter.
#' @param extra An extra covariance parameter used when \code{spcov_type} is
#'   \code{"matern"}, \code{"cauchy"}, \code{"pexponential"}.
#' @param rotate Anisotropy rotation parameter (from 0 to \eqn{\pi} radians) for
#'   the Euclidean portion of the covariance. A value of 0 (the default) implies no rotation.
#' @param scale Anisotropy scale parameter (from 0 to 1) for
#'   the Euclidean portion of the covariance. A value of 1 (the default) implies no scaling.
#' @param nugget The spatially independent (uncorrelated) random error variance. Commonly referred to as
#'   a nugget.
#'
#' @details
#'   Generally, all arguments to \code{tailup_params()}, \code{taildown_params()},
#'   \code{euclid_params()}, and \code{nugget_params()} must be specified, though
#'   default arguments are chosen when the relevant \code{_type} is \code{"none"}.
#'   For full parameterizations of all tailup, taildown, Euclidean, and nugget
#'   covariance functions, see [tailup_initial()], [taildown_initial()],
#'   [euclid_initial()], and [nugget_initial()].
#'
#' @name ssn_params
#'
#' @return A named numeric vector of covariance parameters with class that
#'   matches the relevant \code{type} argument.
#' @export
#'
#' @examples
#' tailup_params("exponential", de = 1, range = 20)
#' taildown_params("exponential", de = 1, range = 20)
#' euclid_params("exponential", de = 1, range = 20, rotate = 0, scale = 1)
#' euclid_params("matern", de = 1, range = 20, extra = 1)
#' nugget_params("nugget", nugget = 1)
#' @references
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach
#' to geostatistical modeling in stream networks. \emph{Ecology} \bold{91(3)},
#' 644--651.
#'
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial
#' statistical models of stream networks (with discussion).
#' \emph{Journal of the American Statistical Association} \bold{105}, 6--18.
#' DOI: 10.1198/jasa.2009.ap08248.  Rejoinder pgs. 22--24.
tailup_params <- function(tailup_type, de, range) {
  check_tailup_type(tailup_type)

  if (tailup_type == "none") {
    de <- 0
    range <- Inf
  } else {
    check_tailup_taildown_parameters(de, range)
  }
  object <- c(de = de, range = range)
  new_object <- structure(object, class = paste("tailup", tailup_type, sep = "_"))
  new_object
}

#' @rdname ssn_params
#' @export
taildown_params <- function(taildown_type, de, range) {
  check_taildown_type(taildown_type)

  if (taildown_type == "none") {
    de <- 0
    range <- Inf
  } else {
    check_tailup_taildown_parameters(de, range)
  }
  object <- c(de = de, range = range)
  new_object <- structure(object, class = paste("taildown", taildown_type, sep = "_"))
  new_object
}

#' @rdname ssn_params
#' @export
euclid_params <- function(euclid_type, de, range, rotate, scale, extra) {
  check_euclid_type(euclid_type)

  if (euclid_type == "none") {
    de <- 0
    range <- Inf
  }

  if (missing(rotate)) {
    rotate <- 0
  }

  if (missing(scale)) {
    scale <- 1
  }
  if (euclid_has_extra(euclid_type)) {
    if (missing(extra)) {
      stop("extra must be specified for this Euclidean covariance type.", call. = FALSE)
    }
  } else if (!missing(extra)) {
    stop("extra is only used by euclid_type = \"matern\", \"cauchy\", or \"pexponential\".", call. = FALSE)
  } else {
    extra <- NULL
  }

  if (euclid_type != "none") {
    check_euclid_extra_parameters(euclid_type, de, range, extra, rotate, scale)
  }

  object <- if (euclid_has_extra(euclid_type)) {
    c(de = de, range = range, extra = extra, rotate = rotate, scale = scale)
  } else {
    c(de = de, range = range, rotate = rotate, scale = scale)
  }
  new_object <- structure(object, class = paste("euclid", euclid_type, sep = "_"))
  new_object
}

#' @rdname ssn_params
#' @export
nugget_params <- function(nugget_type, nugget) {
  check_nugget_type(nugget_type)

  if (nugget_type == "none") {
    nugget <- 0
  } else {
    check_nugget_parameter(nugget)
  }

  object <- c(nugget = nugget)
  new_object <- structure(object, class = paste("nugget", nugget_type, sep = "_"))
  new_object
}

make_euclid_params <- function(euclid_type, de, range, rotate, scale, extra = NULL) {
  if (euclid_has_extra(euclid_type)) {
    euclid_params(euclid_type, de = de, range = range, rotate = rotate, scale = scale, extra = extra)
  } else {
    euclid_params(euclid_type, de = de, range = range, rotate = rotate, scale = scale)
  }
}

get_params_object <- function(classes, cov_orig_val) {
  classes <- remove_covtype(classes)
  tailup_params_val <- tailup_params(
    classes[["tailup"]],
    de = cov_orig_val$orig_ssn[["tailup_de"]],
    range = cov_orig_val$orig_ssn[["tailup_range"]]
  )

  taildown_params_val <- taildown_params(
    classes[["taildown"]],
    de = cov_orig_val$orig_ssn[["taildown_de"]],
    range = cov_orig_val$orig_ssn[["taildown_range"]]
  )

  # class(taildown_params_val) <- classes[["taildown"]]

  euclid_params_val <- make_euclid_params(
    classes[["euclid"]],
    de = cov_orig_val$orig_ssn[["euclid_de"]],
    range = cov_orig_val$orig_ssn[["euclid_range"]],
    rotate = cov_orig_val$orig_ssn[["euclid_rotate"]],
    scale = cov_orig_val$orig_ssn[["euclid_scale"]],
    extra = cov_orig_val$orig_ssn[["euclid_extra"]]
  )

  # class(euclid_params_val) <- classes[["euclid"]]

  nugget_params_val <- nugget_params(
    classes[["nugget"]],
    nugget = cov_orig_val$orig_ssn[["nugget"]]
  )

  # class(nugget_params_val) <- classes[["nugget"]]

  randcov_params_val <- randcov_params(cov_orig_val$orig_randcov)

  params_object <- list(
    tailup = tailup_params_val,
    taildown = taildown_params_val,
    euclid = euclid_params_val,
    nugget = nugget_params_val,
    randcov = randcov_params_val
  )

  params_object
}

#' Get the initial values that are known
#'
#' @param initial_object Initial value object
#'
#' @noRd
get_params_object_known <- function(initial_object) {
  classes <- c(
    tailup = class(initial_object$tailup_initial), taildown = class(initial_object$taildown_initial),
    euclid = class(initial_object$euclid_initial), nugget = class(initial_object$nugget_initial)
  )
  classes <- remove_covtype(classes)

  tailup_params_val <- tailup_params(
    classes[["tailup"]],
    de = initial_object$tailup_initial$initial[["de"]],
    range = initial_object$tailup_initial$initial[["range"]]
  )

  taildown_params_val <- taildown_params(
    classes[["taildown"]],
    de = initial_object$taildown_initial$initial[["de"]],
    range = initial_object$taildown_initial$initial[["range"]]
  )

  euclid_params_val <- make_euclid_params(
    classes[["euclid"]],
    de = initial_object$euclid_initial$initial[["de"]],
    range = initial_object$euclid_initial$initial[["range"]],
    rotate = initial_object$euclid_initial$initial[["rotate"]],
    scale = initial_object$euclid_initial$initial[["scale"]],
    extra = initial_object$euclid_initial$initial[["extra"]]
  )

  nugget_params_val <- nugget_params(
    classes[["nugget"]],
    nugget = initial_object$nugget_initial$initial[["nugget"]]
  )

  randcov_params_val <- randcov_params(initial_object$randcov_initial$initial)

  params_object <- list(
    tailup = tailup_params_val,
    taildown = taildown_params_val,
    euclid = euclid_params_val,
    nugget = nugget_params_val,
    randcov = randcov_params_val
  )

  params_object
}

#' Get the parameter object for glms
#'
#' @param initial_object Initial value object
#' @param cov_orig_val The original covariance parameter values
#'
#' @noRd
get_params_object_glm <- function(classes, cov_orig_val) {
  classes <- remove_covtype(classes)

  tailup_params_val <- tailup_params(
    classes[["tailup"]],
    de = cov_orig_val$orig_ssn[["tailup_de"]],
    range = cov_orig_val$orig_ssn[["tailup_range"]]
  )

  taildown_params_val <- taildown_params(
    classes[["taildown"]],
    de = cov_orig_val$orig_ssn[["taildown_de"]],
    range = cov_orig_val$orig_ssn[["taildown_range"]]
  )

  # class(taildown_params_val) <- classes[["taildown"]]

  euclid_params_val <- make_euclid_params(
    classes[["euclid"]],
    de = cov_orig_val$orig_ssn[["euclid_de"]],
    range = cov_orig_val$orig_ssn[["euclid_range"]],
    rotate = cov_orig_val$orig_ssn[["euclid_rotate"]],
    scale = cov_orig_val$orig_ssn[["euclid_scale"]],
    extra = cov_orig_val$orig_ssn[["euclid_extra"]]
  )

  # class(euclid_params_val) <- classes[["euclid"]]

  nugget_params_val <- nugget_params(
    classes[["nugget"]],
    nugget = cov_orig_val$orig_ssn[["nugget"]]
  )

  # class(nugget_params_val) <- classes[["nugget"]]

  dispersion_params_val <- dispersion_params(
    classes[["dispersion"]],
    dispersion = cov_orig_val$orig_ssn[["dispersion"]]
  )

  randcov_params_val <- randcov_params(cov_orig_val$orig_randcov)

  params_object <- list(
    tailup = tailup_params_val,
    taildown = taildown_params_val,
    euclid = euclid_params_val,
    nugget = nugget_params_val,
    dispersion = dispersion_params_val,
    randcov = randcov_params_val
  )

  params_object
}

#' Get known parameters for glms
#'
#' @param initial_object Initial value object
#'
#' @noRd
get_params_object_glm_known <- function(initial_object) {
  classes <- c(
    tailup = class(initial_object$tailup_initial), taildown = class(initial_object$taildown_initial),
    euclid = class(initial_object$euclid_initial), nugget = class(initial_object$nugget_initial),
    dispersion = class(initial_object$dispersion_initial)
  )
  classes <- remove_covtype(classes)

  tailup_params_val <- tailup_params(
    classes[["tailup"]],
    de = initial_object$tailup_initial$initial[["de"]],
    range = initial_object$tailup_initial$initial[["range"]]
  )

  taildown_params_val <- taildown_params(
    classes[["taildown"]],
    de = initial_object$taildown_initial$initial[["de"]],
    range = initial_object$taildown_initial$initial[["range"]]
  )

  euclid_params_val <- make_euclid_params(
    classes[["euclid"]],
    de = initial_object$euclid_initial$initial[["de"]],
    range = initial_object$euclid_initial$initial[["range"]],
    rotate = initial_object$euclid_initial$initial[["rotate"]],
    scale = initial_object$euclid_initial$initial[["scale"]],
    extra = initial_object$euclid_initial$initial[["extra"]]
  )

  nugget_params_val <- nugget_params(
    classes[["nugget"]],
    nugget = initial_object$nugget_initial$initial[["nugget"]]
  )

  dispersion_params_val <- dispersion_params(
    classes[["dispersion"]],
    dispersion = initial_object$dispersion_initial$initial[["dispersion"]]
  )

  randcov_params_val <- randcov_params(initial_object$randcov_initial$initial)

  params_object <- list(
    tailup = tailup_params_val,
    taildown = taildown_params_val,
    euclid = euclid_params_val,
    nugget = nugget_params_val,
    dispersion = dispersion_params_val,
    randcov = randcov_params_val
  )

  params_object
}

#' Get grid of possible initial values
#'
#' @param cov_grid_vector A configuration of initial covariance parameters
#' @param initial_NA_object The inital NA object (which has information on the known parameters)
#'
#' @noRd
get_params_object_grid <- function(cov_grid_vector, initial_NA_object) {
  classes <- c(
    tailup = class(initial_NA_object$tailup_initial), taildown = class(initial_NA_object$taildown_initial),
    euclid = class(initial_NA_object$euclid_initial), nugget = class(initial_NA_object$nugget_initial)
  )
  classes <- remove_covtype(classes)

  # params object
  tailup_params_val <- tailup_params(
    classes[["tailup"]],
    de = cov_grid_vector[["tailup_de"]],
    range = cov_grid_vector[["tailup_range"]]
  )

  taildown_params_val <- taildown_params(
    classes[["taildown"]],
    de = cov_grid_vector[["taildown_de"]],
    range = cov_grid_vector[["taildown_range"]]
  )

  euclid_params_val <- make_euclid_params(
    classes[["euclid"]],
    de = cov_grid_vector[["euclid_de"]],
    range = cov_grid_vector[["euclid_range"]],
    rotate = cov_grid_vector[["rotate"]],
    scale = cov_grid_vector[["scale"]],
    extra = cov_grid_vector[["euclid_extra"]]
  )

  nugget_params_val <- nugget_params(
    classes[["nugget"]],
    nugget = cov_grid_vector[["nugget"]]
  )

  if (is.null(initial_NA_object$randcov_initial)) {
    randcov_params_val <- NULL
  } else {
    randcov_names <- names(initial_NA_object$randcov_initial$initial)
    randcov_params_val <- randcov_params(cov_grid_vector[randcov_names])
  }

  params_object <- list(
    tailup = tailup_params_val,
    taildown = taildown_params_val,
    euclid = euclid_params_val,
    nugget = nugget_params_val,
    randcov = randcov_params_val
  )
}

#' Get grid of possible initial values for glms
#'
#' @param cov_grid_vector A configuration of initial covariance parameters
#' @param initial_NA_object The inital NA object (which has information on the known parameters)
#'
#' @noRd
get_params_object_grid_glm <- function(cov_grid_vector, initial_NA_object) {
  classes <- c(
    tailup = class(initial_NA_object$tailup_initial), taildown = class(initial_NA_object$taildown_initial),
    euclid = class(initial_NA_object$euclid_initial), nugget = class(initial_NA_object$nugget_initial),
    dispersion = class(initial_NA_object$dispersion_initial)
  )
  classes <- remove_covtype(classes)

  # params object
  tailup_params_val <- tailup_params(
    classes[["tailup"]],
    de = cov_grid_vector[["tailup_de"]],
    range = cov_grid_vector[["tailup_range"]]
  )

  taildown_params_val <- taildown_params(
    classes[["taildown"]],
    de = cov_grid_vector[["taildown_de"]],
    range = cov_grid_vector[["taildown_range"]]
  )

  euclid_params_val <- make_euclid_params(
    classes[["euclid"]],
    de = cov_grid_vector[["euclid_de"]],
    range = cov_grid_vector[["euclid_range"]],
    rotate = cov_grid_vector[["rotate"]],
    scale = cov_grid_vector[["scale"]],
    extra = cov_grid_vector[["euclid_extra"]]
  )

  nugget_params_val <- nugget_params(
    classes[["nugget"]],
    nugget = cov_grid_vector[["nugget"]]
  )

  dispersion_params_val <- dispersion_params(
    classes[["dispersion"]],
    dispersion = cov_grid_vector[["dispersion"]]
  )

  if (is.null(initial_NA_object$randcov_initial)) {
    randcov_params_val <- NULL
  } else {
    randcov_names <- names(initial_NA_object$randcov_initial$initial)
    randcov_params_val <- randcov_params(cov_grid_vector[randcov_names])
  }

  params_object <- list(
    tailup = tailup_params_val,
    taildown = taildown_params_val,
    euclid = euclid_params_val,
    nugget = nugget_params_val,
    dispersion = dispersion_params_val,
    randcov = randcov_params_val
  )
}
