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

check_tailup_type <- function(tailup_type) {
  tailup_valid <- c("linear", "spherical", "exponential", "mariah", "epa", "none")

  if (!(tailup_type %in% tailup_valid)) {
    stop(paste(tailup_type, " is not a valid tailup covariance function."), call. = FALSE)
  }
}

check_taildown_type <- function(taildown_type) {
  taildown_valid <- c("linear", "spherical", "exponential", "mariah", "epa", "none")

  if (!(taildown_type %in% taildown_valid)) {
    stop(paste(taildown_type, "is not a valid taildown covariance function."), call. = FALSE)
  }
}

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

check_nugget_type <- function(nugget_type) {
  nugget_valid <- c("nugget", "none")
  if (!(nugget_type %in% nugget_valid)) {
    stop(paste(nugget_type, "is not a valid nugget covariance function."), call. = FALSE)
  }
}

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

# check if whole number
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
