#' An initial value object with NA values that indicate where estimation is required
#'
#' @param initial_object Initial value object
#' @param data_object Data object
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

#' Tailup initial NA object
#'
#' @param initial An initial value specification
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

#' Taildown initial NA object
#'
#' @param initial An initial value specification
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

#' Euclidean initial NA object
#'
#' @param initial An initial value specification
#'
#' @noRd
euclid_initial_NA <- function(initial, data_object) {
  euclid_names <- c("de", "range", "rotate", "scale")

  if (inherits(initial, "euclid_none")) {
    # set defaults if none covariance
    euclid_val_default <- c(de = 0, range = Inf, rotate = 0, scale = 1)
    euclid_known_default <- c(de = TRUE, range = TRUE, rotate = TRUE, scale = TRUE)
  } else {
    if (data_object$anisotropy) {
      # otherwise we will pick them
      euclid_val_default <- c(de = NA, range = NA, rotate = NA, scale = NA)
      euclid_known_default <- c(de = FALSE, range = FALSE, rotate = FALSE, scale = FALSE)
    } else {
      # otherwise we will pick them (but fix anisotropy parameters)
      euclid_val_default <- c(de = NA, range = NA, rotate = 0, scale = 1)
      euclid_known_default <- c(de = FALSE, range = FALSE, rotate = TRUE, scale = TRUE)
    }
  }
  # substitute known values
  new_initial <- insert_initial_NA(euclid_names, euclid_val_default, euclid_known_default, initial)
  new_initial
}

#' Nugget initial NA object
#'
#' @param initial An initial value specification
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

#' Insert NA values when covariance parameters assumed unknown
#'
#' @param names The names to assume unknown
#' @param val_default The default value of values (NA)
#' @param known_default The default value of whether parameters are known
#' @param initial An initial value specification
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

#' Insert NA values when random effect parameters assumed unknown
#'
#' @param randcov_initial Random effect initial object
#' @param data_object Data object
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

#' An initial value object with NA values that indicate where estimation is required for glms
#'
#' @param initial_object Initial value object
#' @param data_object Data object
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

#' Insert NA values when dispersion parameters assumed unknown
#'
#' @param initial Initial value object
#' @param data_object Data object
#'
#' @noRd
dispersion_initial_NA <- function(initial, data_object) {
  dispersion_names <- c("dispersion")

  if (data_object$family %in% c("poisson", "binomial")) {
    new_initial <- dispersion_initial(data_object$family, 1, known = "dispersion")
  } else {
    dispersion_val_default <- c(dispersion = NA)
    dispersion_known_default <- c(dispersion = FALSE)
    new_initial <- insert_initial_NA(dispersion_names, dispersion_val_default, dispersion_known_default, initial)
  }
  new_initial
}
