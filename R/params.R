#' Create covariance parameter objects.
#'
#' @description Create a covariance parameter object for us with other functions.
#'   See [spmodel::randcov_params()] for documentation regarding
#'   random effect covariance parameter objects.
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
#'
#' @name ssn_params
#'
#' @return A parameter object with class that matches the relevant \code{type} argument.
#' @export
#'
#' @examples
#' tailup_params("exponential", de = 1, range = 20)
#' taildown_params("exponential", de = 1, range = 20)
#' euclid_params("exponential", de = 1, range = 20, rotate = 0, scale = 1)
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
  }
  object <- c(de = de, range = range)
  new_object <- structure(object, class = paste("taildown", taildown_type, sep = "_"))
  new_object
}

#' @rdname ssn_params
#' @export
euclid_params <- function(euclid_type, de, range, rotate, scale) {
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

  object <- c(de = de, range = range, rotate = rotate, scale = scale)
  new_object <- structure(object, class = paste("euclid", euclid_type, sep = "_"))
  new_object
}

#' @rdname ssn_params
#' @export
nugget_params <- function(nugget_type, nugget) {
  check_nugget_type(nugget_type)

  if (nugget_type == "none") {
    nugget <- 0
  }

  object <- c(nugget = nugget)
  new_object <- structure(object, class = paste("nugget", nugget_type, sep = "_"))
  new_object
}

get_params_object <- function(classes, cov_orig_val) {
  classes <- remove_covtype(classes)
  tailup_params_val <- tailup_params(
    classes[["tailup"]],
    de = cov_orig_val$orig_ssn[["tailup_de"]],
    range = cov_orig_val$orig_ssn[["tailup_range"]]
  )

  # class(tailup_params_val) <- remove_covtype(tailup_params_val)
  # replace class as it has tailup_tailup_exponential structure
  # and don't want to edit user side *_params functions
  # class(tailup_params_val) <- classes[["tailup"]]

  taildown_params_val <- taildown_params(
    classes[["taildown"]],
    de = cov_orig_val$orig_ssn[["taildown_de"]],
    range = cov_orig_val$orig_ssn[["taildown_range"]]
  )

  # class(taildown_params_val) <- classes[["taildown"]]

  euclid_params_val <- euclid_params(
    classes[["euclid"]],
    de = cov_orig_val$orig_ssn[["euclid_de"]],
    range = cov_orig_val$orig_ssn[["euclid_range"]],
    rotate = cov_orig_val$orig_ssn[["euclid_rotate"]],
    scale = cov_orig_val$orig_ssn[["euclid_scale"]]
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

  euclid_params_val <- euclid_params(
    classes[["euclid"]],
    de = initial_object$euclid_initial$initial[["de"]],
    range = initial_object$euclid_initial$initial[["range"]],
    rotate = initial_object$euclid_initial$initial[["rotate"]],
    scale = initial_object$euclid_initial$initial[["scale"]]
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

  # class(tailup_params_val) <- remove_covtype(tailup_params_val)
  # replace class as it has tailup_tailup_exponential structure
  # and don't want to edit user side *_params functions
  # class(tailup_params_val) <- classes[["tailup"]]

  taildown_params_val <- taildown_params(
    classes[["taildown"]],
    de = cov_orig_val$orig_ssn[["taildown_de"]],
    range = cov_orig_val$orig_ssn[["taildown_range"]]
  )

  # class(taildown_params_val) <- classes[["taildown"]]

  euclid_params_val <- euclid_params(
    classes[["euclid"]],
    de = cov_orig_val$orig_ssn[["euclid_de"]],
    range = cov_orig_val$orig_ssn[["euclid_range"]],
    rotate = cov_orig_val$orig_ssn[["euclid_rotate"]],
    scale = cov_orig_val$orig_ssn[["euclid_scale"]]
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

  euclid_params_val <- euclid_params(
    classes[["euclid"]],
    de = initial_object$euclid_initial$initial[["de"]],
    range = initial_object$euclid_initial$initial[["range"]],
    rotate = initial_object$euclid_initial$initial[["rotate"]],
    scale = initial_object$euclid_initial$initial[["scale"]]
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

  euclid_params_val <- euclid_params(
    classes[["euclid"]],
    de = cov_grid_vector[["euclid_de"]],
    range = cov_grid_vector[["euclid_range"]],
    rotate = cov_grid_vector[["rotate"]],
    scale = cov_grid_vector[["scale"]]
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

  euclid_params_val <- euclid_params(
    classes[["euclid"]],
    de = cov_grid_vector[["euclid_de"]],
    range = cov_grid_vector[["euclid_range"]],
    rotate = cov_grid_vector[["rotate"]],
    scale = cov_grid_vector[["scale"]]
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
