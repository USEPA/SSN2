#' Various model checks for ssn_lm
#'
#' @param initial_object Initial value object.
#' @param ssn.object SSN object.
#' @param additive Additive function value name.
#' @param estmethod Estimation method.
#'
#' @noRd
check_ssn_lm <- function(initial_object, ssn.object, additive, estmethod) {
  if (is.null(additive)) {
    if (!grepl("none", class(initial_object$tailup))) {
      stop("Argument additive must be specified.", call. = FALSE)
    }
  } else {
    if (!(additive %in% names(ssn.object$obs))) {
      stop("additive column not found in ssn.object", call. = FALSE)
    }
  }

  if (!estmethod %in% c("reml", "ml")) {
    stop("Estimation method must be \"reml\" or \"ml\".", call. = FALSE)
  }
}

#' Various model checks for ssn_glm
#'
#' @param initial_object Initial value object.
#' @param ssn.object SSN object.
#' @param additive Additive function value name.
#' @param estmethod Estimation method.
#'
#' @noRd
check_ssn_glm <- function(initial_object, ssn.object, additive, estmethod) {
  if (is.null(additive)) {
    if (!grepl("none", class(initial_object$tailup_initial))) {
      stop("Argument additive must be specified.", call. = FALSE)
    }
  } else {
    if (!(additive %in% names(ssn.object$obs))) {
      stop("additive column not found in ssn.object", call. = FALSE)
    }
  }

  if (!estmethod %in% c("reml", "ml")) {
    stop("Estimation method must be \"reml\" or \"ml\".", call. = FALSE)
  }
}


#' Check for valid tailup type
#'
#' @param tailup_type The tailup covariance type.
#'
#' @noRd
check_tailup_type <- function(tailup_type) {
  tailup_valid <- c("linear", "spherical", "exponential", "mariah", "epa", "gaussian", "none")

  if (!(tailup_type %in% tailup_valid)) {
    stop(paste(tailup_type, " is not a valid tailup covariance function."), call. = FALSE)
  }
}

#' Check for valid taildown type
#'
#' @param taildown_type The taildown covariance type.
#'
#' @noRd
check_taildown_type <- function(taildown_type) {
  taildown_valid <- c("linear", "spherical", "exponential", "mariah", "epa", "gaussian", "none")

  if (!(taildown_type %in% taildown_valid)) {
    stop(paste(taildown_type, "is not a valid taildown covariance function."), call. = FALSE)
  }
}

#' Check for valid Euclidean type
#'
#' @param euclid_type The Euclidean covariance type.
#'
#' @noRd
check_euclid_type <- function(euclid_type) {
  euclid_valid <- c(
    "spherical", "exponential", "gaussian", "cosine",
    "cubic", "pentaspherical", "wave", "jbessel", "gravity",
    "rquad", "magnetic", "none"
  )

  if (!(euclid_type %in% euclid_valid)) {
    stop(paste(euclid_type, "is not a valid Euclidean covariance function."), call. = FALSE)
  }
}

#' Check for valid nugget type
#'
#' @param nugget_type The nugget covariance type.
#'
#' @noRd
check_nugget_type <- function(nugget_type) {
  nugget_valid <- c("nugget", "none")
  if (!(nugget_type %in% nugget_valid)) {
    stop(paste(nugget_type, "is not a valid nugget covariance function."), call. = FALSE)
  }
}

#' Check for valid ssn_glm family and dispersion parameter
#'
#' @param family The ssn_glm family.
#' @param dispersion The dispersion parameter.
#'
#' @noRd
check_dispersion <- function(family, dispersion) {
  # family must be a character here
  family_valid <- c("binomial", "poisson", "nbinomial", "Gamma", "inverse.gaussian", "beta")
  if (!(family %in% family_valid)) {
    stop(paste(family, " is not a valid glm family.", sep = ""), call. = FALSE)
  }

  # dispersion can't be missing
  if (!is.null(dispersion) && dispersion != 1 && family %in% c("binomial", "poisson")) {
    stop(paste(family, "dispersion parameter must be fixed at one."), call. = FALSE)
  }
}

#' Various checks on the response variable in ssn_glm
#'
#' @param family The ssn_glm family.
#' @param y The response variable.
#' @param size The number of trials (if family is binomial)
#'
#' @noRd
response_checks_glm <- function(family, y, size) {
  # checks on y
  if (family == "binomial") {
    if (any(size < 1)) {
      stop("All size values must be at least 1.", call. = FALSE)
    }

    if (any(!is.wholenumber(size))) {
      stop("All size values must be a whole number.", call. = FALSE)
    }

    if (any(y < 0)) {
      stop("All response values must be at least 0.", call. = FALSE)
    }

    if (any(!is.wholenumber(y))) {
      stop("All response values must be a whole number.", call. = FALSE)
    }

    if (all(size == 1)) {
      if (!all(y == 0 | y == 1)) {
        stop("All response values must be 0 or 1. 0 indicates a failure and 1 indicates a success.", call. = FALSE)
      }
    }
  } else if (family == "beta") {
    if (any(y <= 0 | y >= 1)) {
      stop("All response values must be greater than 0 and less than 1.", call. = FALSE)
    }
  } else if (family %in% c("poisson", "nbinomial")) {
    if (any(y < 0)) {
      stop("All response values must be at least 0.", call. = FALSE)
    }

    if (any(!is.wholenumber(y))) {
      stop("All response values must be a whole number.", call. = FALSE)
    }
  } else if (family %in% c("Gamma", "inverse.gaussian")) {
    if (any(y <= 0)) {
      stop("All response values must be greater than 0.", call. = FALSE)
    }
  }
}


#' Check if value is a whole number.
#'
#' @param x A vector.
#' @param tol Tolerance to check whether x is a whole number.
#'
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
