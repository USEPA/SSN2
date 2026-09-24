check_formula_vars_in_data <- function(formula, data, random = NULL, partition_factor = NULL) {
  if ("." %in% all.vars(random)) {
    stop("The `.` shorthand is not supported in random. Explicitly list the desired variable(s).", call. = FALSE)
  }
  if ("." %in% all.vars(partition_factor)) {
    stop("The `.` shorthand is not supported in partition_factor. Explicitly list the desired variable(s).", call. = FALSE)
  }
  formula_vars <- unique(c(all.vars(formula), all.vars(random), all.vars(partition_factor)))
  formula_vars <- setdiff(formula_vars, ".")
  missing_vars <- setdiff(formula_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "Variable(s) ", paste0("\"", missing_vars, "\"", collapse = ", "),
      " used in formula, random, or partition_factor not found in data.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

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
  tailup_valid <- c("linear", "spherical", "exponential", "mariah", "epa", "gaussian", "none")

  if (!(tailup_type %in% tailup_valid)) {
    stop(paste(tailup_type, " is not a valid tailup covariance function."), call. = FALSE)
  }
}

check_taildown_type <- function(taildown_type) {
  taildown_valid <- c("linear", "spherical", "exponential", "mariah", "epa", "gaussian", "none")

  if (!(taildown_type %in% taildown_valid)) {
    stop(paste(taildown_type, "is not a valid taildown covariance function."), call. = FALSE)
  }
}

check_euclid_type <- function(euclid_type) {
  if (identical(euclid_type, "cosine")) {
    stop("'cosine' is no longer a supported Euclidean covariance name. Use 'circular' for the covariance previously named 'cosine' in SSN2.", call. = FALSE)
  }
  euclid_valid <- c(
    "spherical", "exponential", "gaussian", "circular",
    "cubic", "pentaspherical", "wave", "jbessel", "gravity",
    "rquad", "magnetic", "matern", "cauchy", "pexponential", "none"
  )

  if (!(euclid_type %in% euclid_valid)) {
    stop(paste(euclid_type, "is not a valid Euclidean covariance function."), call. = FALSE)
  }
}

euclid_has_extra <- function(euclid_type) {
  euclid_type %in% c("matern", "cauchy", "pexponential")
}

check_euclid_parameter <- function(value, name, lower = -Inf, upper = Inf,
                                   lower_open = FALSE, allow_na = FALSE, allow_inf = FALSE) {
  if (is.null(value)) return(invisible(NULL))
  if (allow_na && length(value) == 1 && is.na(value)) return(invisible(NULL))
  if (!is.numeric(value) || length(value) != 1 || is.na(value)) {
    stop(name, " must be a single finite numeric value.", call. = FALSE)
  }
  if (!allow_inf && !is.finite(value)) {
    stop(name, " must be a single finite numeric value.", call. = FALSE)
  }
  below <- if (lower_open) value <= lower else value < lower
  if (below || value > upper) {
    lower_text <- if (lower_open) "greater than" else "at least"
    if (is.finite(upper)) {
      stop(name, " must be ", lower_text, " ", lower, " and at most ", upper, ".", call. = FALSE)
    }
    stop(name, " must be ", lower_text, " ", lower, ".", call. = FALSE)
  }
  invisible(NULL)
}

check_tailup_taildown_parameters <- function(de, range, allow_na = FALSE) {
  check_euclid_parameter(de, "de", lower = 0, allow_na = allow_na)
  # range = Inf is a legitimate, pervasively-used sentinel for the "none"
  # covariance type (no spatial dependence), not just a user-facing bound
  check_euclid_parameter(range, "range", lower = 0, lower_open = TRUE, allow_na = allow_na, allow_inf = TRUE)
  invisible(NULL)
}

check_nugget_parameter <- function(nugget, allow_na = FALSE) {
  check_euclid_parameter(nugget, "nugget", lower = 0, allow_na = allow_na)
  invisible(NULL)
}

check_euclid_extra_parameters <- function(euclid_type, de, range, extra,
                                          rotate, scale, allow_na = FALSE) {
  check_euclid_parameter(de, "de", lower = 0, allow_na = allow_na)
  # range = Inf is a legitimate, pervasively-used sentinel for the "none"
  # covariance type (no spatial dependence), not just a user-facing bound
  check_euclid_parameter(range, "range", lower = 0, lower_open = TRUE, allow_na = allow_na, allow_inf = TRUE)
  check_euclid_parameter(rotate, "rotate", lower = 0, upper = pi, allow_na = allow_na)
  check_euclid_parameter(scale, "scale", lower = 0, upper = 1, lower_open = TRUE, allow_na = allow_na)

  if (euclid_type == "matern") {
    check_euclid_parameter(extra, "extra", lower = 0.2, upper = 5, allow_na = allow_na)
  } else if (euclid_type == "cauchy") {
    check_euclid_parameter(extra, "extra", lower = 0, lower_open = TRUE, allow_na = allow_na)
  } else if (euclid_type == "pexponential") {
    check_euclid_parameter(extra, "extra", lower = 0, upper = 2, lower_open = TRUE, allow_na = allow_na)
  }
  invisible(NULL)
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

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

# a prediction interval needs a location that has not been observed; without
# newdata, augment() describes the observed data, where only a confidence
# interval (around the fitted mean) is meaningful -- shared by augment.ssn_lm()
# and augment.ssn_glm() so both degrade the combination identically, matching
# spmodel's check_interval_augment() (warn + reset to "none" rather than error)
check_interval_augment <- function(interval, newdata_given) {
  if (!newdata_given && interval == "prediction") {
    warning(
      "interval = \"prediction\" is ignored when newdata is not supplied, because a prediction interval ",
      "requires a location that has not been observed. Supply newdata for prediction intervals, or use ",
      "interval = \"confidence\" for an interval around the fitted mean at the observed locations.",
      call. = FALSE
    )
    interval <- "none"
  }
  interval
}
