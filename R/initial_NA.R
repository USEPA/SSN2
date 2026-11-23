#' Fill in default (NA, unknown) covariance initial values
#'
#' @param initial_object A covariance initial object from \code{get_initial_object()}
#' @param data_object The data object
#'
#' @return An initial object with NA values (to be estimated) filled in for
#'   any tailup, taildown, euclid, nugget, and (if relevant) randcov
#'   parameter not otherwise given an initial or known value
#'
#' @noRd
get_initial_NA_object <- function(initial_object, data_object) {
  # get each initial NA object
  tailup_initial_NA_val <- tailup_initial_NA(initial_object$tailup_initial)
  taildown_initial_NA_val <- taildown_initial_NA(initial_object$taildown_initial)
  euclid_initial_NA_val <- euclid_initial_NA(initial_object$euclid_initial, data_object)
  nugget_initial_NA_val <- nugget_initial_NA(initial_object$nugget_initial)
  randcov_initial_NA_val <- randcov_initial_NA(initial_object$randcov_initial, data_object)

  # put them in a relevant list
  initial_NA_object <- list(
    tailup_initial = tailup_initial_NA_val,
    taildown_initial = taildown_initial_NA_val,
    euclid_initial = euclid_initial_NA_val,
    nugget_initial = nugget_initial_NA_val,
    randcov_initial = randcov_initial_NA_val
  )

  # return all of them
  initial_NA_object
}

#' Fill in default (NA, unknown) tailup covariance initial values
#'
#' @param initial A \code{tailup_initial} object
#'
#' @return A \code{tailup_initial} object with NA values when relevant (which
#'   are replaced later) -- values are NA if they want us to pick those
#'   initial values
#'
#' @noRd
tailup_initial_NA <- function(initial) {
  tailup_names <- c("de", "range")

  if (inherits(initial, "tailup_none")) {
    # set defaults if none covariance
    tailup_val_default <- c(de = 0, range = Inf)
    tailup_known_default <- c(de = TRUE, range = TRUE)
  } else {
    # otherwise we will pick them
    tailup_val_default <- c(de = NA, range = NA)
    tailup_known_default <- c(de = FALSE, range = FALSE)
  }
  # substitute known values
  new_initial <- insert_initial_NA(tailup_names, tailup_val_default, tailup_known_default, initial)
  new_initial
}

#' Fill in default (NA, unknown) taildown covariance initial values
#'
#' @param initial A \code{taildown_initial} object
#'
#' @return A \code{taildown_initial} object with NA values when relevant
#'   (which are replaced later) -- values are NA if they want us to pick
#'   those initial values
#'
#' @noRd
taildown_initial_NA <- function(initial) {
  taildown_names <- c("de", "range")

  if (inherits(initial, "taildown_none")) {
    # set defaults if none covariance
    taildown_val_default <- c(de = 0, range = Inf)
    taildown_known_default <- c(de = TRUE, range = TRUE)
  } else {
    # otherwise we will pick them
    taildown_val_default <- c(de = NA, range = NA)
    taildown_known_default <- c(de = FALSE, range = FALSE)
  }
  # substitute known values
  new_initial <- insert_initial_NA(taildown_names, taildown_val_default, taildown_known_default, initial)
  new_initial
}

#' Fill in default (NA, unknown) euclid covariance initial values
#'
#' @param initial A \code{euclid_initial} object
#' @param data_object The data object, whose \code{anisotropy} element
#'   determines whether \code{rotate} and \code{scale} are estimated or fixed
#'   at their no-op defaults (0 rotation, scale 1)
#'
#' @return A \code{euclid_initial} object with NA values when relevant (which
#'   are replaced later) -- values are NA if they want us to pick those
#'   initial values
#'
#' @noRd
euclid_initial_NA <- function(initial, data_object) {
  has_extra <- inherits(initial, c("euclid_matern", "euclid_cauchy", "euclid_pexponential"))
  euclid_names <- if (has_extra) c("de", "range", "extra", "rotate", "scale") else c("de", "range", "rotate", "scale")

  if (inherits(initial, "euclid_none")) {
    # set defaults if none covariance
    euclid_val_default <- c(de = 0, range = Inf, rotate = 0, scale = 1)
    euclid_known_default <- c(de = TRUE, range = TRUE, rotate = TRUE, scale = TRUE)
  } else {
    if (data_object$anisotropy) {
      # otherwise we will pick them
      euclid_val_default <- if (has_extra) c(de = NA, range = NA, extra = NA, rotate = NA, scale = NA) else c(de = NA, range = NA, rotate = NA, scale = NA)
      euclid_known_default <- if (has_extra) c(de = FALSE, range = FALSE, extra = FALSE, rotate = FALSE, scale = FALSE) else c(de = FALSE, range = FALSE, rotate = FALSE, scale = FALSE)
    } else {
      # otherwise we will pick them (but fix anisotropy parameters)
      euclid_val_default <- if (has_extra) c(de = NA, range = NA, extra = NA, rotate = 0, scale = 1) else c(de = NA, range = NA, rotate = 0, scale = 1)
      euclid_known_default <- if (has_extra) c(de = FALSE, range = FALSE, extra = FALSE, rotate = TRUE, scale = TRUE) else c(de = FALSE, range = FALSE, rotate = TRUE, scale = TRUE)
    }
  }
  # substitute known values
  new_initial <- insert_initial_NA(euclid_names, euclid_val_default, euclid_known_default, initial)
  new_initial
}

#' Fill in default (NA, unknown) nugget covariance initial values
#'
#' @param initial A \code{nugget_initial} object
#'
#' @return A \code{nugget_initial} object with NA values when relevant (which
#'   are replaced later) -- values are NA if they want us to pick those
#'   initial values
#'
#' @noRd
nugget_initial_NA <- function(initial) {
  nugget_names <- c("nugget")

  if (inherits(initial, "nugget_none")) {
    # set defaults if none covariance
    nugget_val_default <- c(nugget = 0)
    nugget_known_default <- c(nugget = TRUE)
  } else {
    # otherwise we will pick them
    nugget_val_default <- c(nugget = NA)
    nugget_known_default <- c(nugget = FALSE)
  }
  # substitute known values
  new_initial <- insert_initial_NA(nugget_names, nugget_val_default, nugget_known_default, initial)
  new_initial
}

#' Substitute default NA values into an initial object for parameters not
#'   otherwise given an initial or known value
#'
#' @param names The names of all parameters for this covariance type
#' @param val_default Default values (NA to estimate, or a fixed value) for
#'   each name in \code{names}
#' @param known_default Default \code{is_known} values for each name in
#'   \code{names}
#' @param initial A partially-specified initial object
#'
#' @return \code{initial} with defaults substituted in for any parameter not
#'   already specified, reordered to match \code{names}
#'
#' @noRd
insert_initial_NA <- function(names, val_default, known_default, initial) {
  # find names with known initial values
  names_replace <- setdiff(names, names(initial$initial))
  # replace other values with NA defaults
  initial$initial[names_replace] <- val_default[names_replace]
  initial$is_known[names_replace] <- known_default[names_replace]

  # reorder names in initial object (with some value for all parameters)
  initial$initial <- initial$initial[names]
  initial$is_known <- initial$is_known[names]

  initial
}

#' Fill random effect parameters with NA's and known FALSE if specified in
#'   the model's \code{random} formula but not given an initial value
#'
#' @param randcov_initial A \code{randcov_initial} object
#' @param data_object The data object, whose \code{randcov_names} gives the
#'   names of the random effects specified in the model
#'
#' @return A \code{randcov_initial} object with appropriate NA's
#'
#' @noRd
randcov_initial_NA <- function(randcov_initial, data_object) {
  if (is.null(randcov_initial)) {
    randcov_initial <- NULL
  } else {
    randcov_names <- data_object$randcov_names
    randcov_val_default <- rep(NA, length = length(randcov_names))
    names(randcov_val_default) <- randcov_names
    randcov_known_default <- rep(FALSE, length = length(randcov_names))
    names(randcov_known_default) <- randcov_names
    # find names not in initial
    randcov_out <- setdiff(randcov_names, names(randcov_initial$initial))
    # put in values not in initial
    randcov_initial$initial[randcov_out] <- randcov_val_default[randcov_out]
    # put in is_known not in initial
    randcov_initial$is_known[randcov_out] <- randcov_known_default[randcov_out]
    # reorder names
    randcov_initial$initial <- randcov_initial$initial[randcov_names]
    randcov_initial$is_known <- randcov_initial$is_known[randcov_names]
  }

  # return randcov_initial
  randcov_initial
}

#' Fill in default (NA, unknown) covariance and dispersion initial values for
#'   GLM-type models
#'
#' @param initial_object A covariance initial object from
#'   \code{get_initial_object_glm()}
#' @param data_object The data object
#'
#' @return An initial object with NA values (to be estimated) filled in for
#'   any tailup, taildown, euclid, nugget, dispersion, and (if relevant)
#'   randcov parameter not otherwise given an initial or known value
#'
#' @noRd
get_initial_NA_object_glm <- function(initial_object, data_object) {
  tailup_initial_NA_val <- tailup_initial_NA(initial_object$tailup_initial)
  taildown_initial_NA_val <- taildown_initial_NA(initial_object$taildown_initial)
  euclid_initial_NA_val <- euclid_initial_NA(initial_object$euclid_initial, data_object)
  nugget_initial_NA_val <- nugget_initial_NA(initial_object$nugget_initial)
  dispersion_initial_NA_val <- dispersion_initial_NA(initial_object$dispersion_initial, data_object)
  randcov_initial_NA_val <- randcov_initial_NA(initial_object$randcov_initial, data_object)

  initial_NA_object <- list(
    tailup_initial = tailup_initial_NA_val,
    taildown_initial = taildown_initial_NA_val,
    euclid_initial = euclid_initial_NA_val,
    nugget_initial = nugget_initial_NA_val,
    dispersion_initial = dispersion_initial_NA_val,
    randcov_initial = randcov_initial_NA_val
  )

  initial_NA_object
}

#' Fill in default (NA, unknown) dispersion initial values
#'
#' @param initial A \code{dispersion_initial} object
#' @param data_object The data object
#'
#' @return A \code{dispersion_initial} object with any missing initial value
#'   and \code{is_known} indicator filled in with defaults (fixed at one for
#'   the binomial and Poisson families, otherwise unknown)
#'
#' @noRd
dispersion_initial_NA <- function(initial, data_object) {
  dispersion_names <- c("dispersion")

  # poisson and binomial dispersion is not identifiable, so it is always
  # fixed at one regardless of what (if anything) the user supplied
  if (data_object$family %in% c("poisson", "binomial")) {
    new_initial <- dispersion_initial(data_object$family, 1, known = "dispersion")
  } else {
    # any dispersion value the user did not specify defaults to NA (to be
    # estimated) and unknown, rather than erroring
    dispersion_val_default <- c(dispersion = NA)
    dispersion_known_default <- c(dispersion = FALSE)
    new_initial <- insert_initial_NA(dispersion_names, dispersion_val_default, dispersion_known_default, initial)
  }
  new_initial
}
