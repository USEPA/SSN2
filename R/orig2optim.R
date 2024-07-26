#' Transform parameters from original to optim scale
#'
#' @param initial_object Initial value object
#'
#' @noRd
orig2optim <- function(initial_object) {
  # find parameters on optim (transformed) scale

  ## tailup
  tailup_de <- initial_object$tailup_initial$initial[["de"]]
  tailup_de_log <- log(tailup_de)
  tailup_de_is_known <- initial_object$tailup_initial$is_known[["de"]]

  tailup_range <- initial_object$tailup_initial$initial[["range"]]
  tailup_range_log <- log(tailup_range)
  tailup_range_is_known <- initial_object$tailup_initial$is_known[["range"]]

  ## taildown
  taildown_de <- initial_object$taildown_initial$initial[["de"]]
  taildown_de_log <- log(taildown_de)
  taildown_de_is_known <- initial_object$taildown_initial$is_known[["de"]]

  taildown_range <- initial_object$taildown_initial$initial[["range"]]
  taildown_range_log <- log(taildown_range)
  taildown_range_is_known <- initial_object$taildown_initial$is_known[["range"]]

  ## euclid
  euclid_de <- initial_object$euclid_initial$initial[["de"]]
  euclid_de_log <- log(euclid_de)
  euclid_de_is_known <- initial_object$euclid_initial$is_known[["de"]]

  euclid_range <- initial_object$euclid_initial$initial[["range"]]
  euclid_range_log <- log(euclid_range)
  euclid_range_is_known <- initial_object$euclid_initial$is_known[["range"]]

  euclid_rotate <- initial_object$euclid_initial$initial[["rotate"]]
  euclid_rotate_prop <- euclid_rotate / pi
  euclid_rotate_logodds <- logit(euclid_rotate_prop)
  euclid_rotate_is_known <- initial_object$euclid_initial$is_known[["rotate"]]

  euclid_scale <- initial_object$euclid_initial$initial[["scale"]]
  euclid_scale_logodds <- logit(euclid_scale)
  euclid_scale_is_known <- initial_object$euclid_initial$is_known[["scale"]]

  ## nugget
  nugget <- initial_object$nugget_initial$initial[["nugget"]]
  nugget_log <- log(nugget)
  nugget_is_known <- initial_object$nugget_initial$is_known[["nugget"]]

  # make transformed parameter vectors
  orig2optim_val <- c(
    tailup_de_log = tailup_de_log,
    tailup_range_log = tailup_range_log,
    taildown_de_log = taildown_de_log,
    taildown_range_log = taildown_range_log,
    euclid_de_log = euclid_de_log,
    euclid_range_log = euclid_range_log,
    euclid_rotate_logodds = euclid_rotate_logodds,
    euclid_scale_logodds = euclid_scale_logodds,
    nugget_log = nugget_log
  )

  # make transformed is known vector
  orig2optim_is_known <- c(
    tailup_de_is_known = tailup_de_is_known,
    tailup_range_is_known = tailup_range_is_known,
    taildown_de_is_known = taildown_de_is_known,
    taildown_range_is_known = taildown_range_is_known,
    euclid_de_is_known = euclid_de_is_known,
    euclid_range_is_known = euclid_range_is_known,
    euclid_rotate_is_known = euclid_rotate_is_known,
    euclid_scale_is_known = euclid_scale_is_known,
    nugget_is_known = nugget_is_known
  )

  # find the number of unknown parameters
  n_est_ssn <- sum(!orig2optim_is_known)

  # keep values bounded on inverse transformation scale by setting max/mins
  orig2optim_val <- ifelse(orig2optim_val > 50 & !orig2optim_is_known, 50, orig2optim_val)
  orig2optim_val <- ifelse(orig2optim_val < -50 & !orig2optim_is_known, -50, orig2optim_val)


  # handle random effects
  if (!is.null(initial_object$randcov_initial)) {
    randcov_orig2optim_val <- log(initial_object$randcov_initial$initial)
    names(randcov_orig2optim_val) <- paste(names(initial_object$randcov_initial$initial), "log", sep = "_")
    randcov_orig2optim_is_known <- initial_object$randcov_initial$is_known
    names(randcov_orig2optim_is_known) <- paste(names(initial_object$randcov_initial$is_known), "log", sep = "_")
    n_est_rand <- sum(!randcov_orig2optim_is_known)
    randcov_orig2optim_val <- ifelse(randcov_orig2optim_val > 50 & !randcov_orig2optim_is_known, 50, randcov_orig2optim_val)
    randcov_orig2optim_val <- ifelse(randcov_orig2optim_val < -50 & !randcov_orig2optim_is_known, -50, randcov_orig2optim_val)
  } else {
    randcov_orig2optim_val <- NULL
    randcov_orig2optim_is_known <- NULL
    n_est_rand <- 0
  }

  # return the orig2optim initial value (and other) information
  orig2optim_list <- list(
    value = orig2optim_val,
    is_known = orig2optim_is_known,
    n_est_ssn = n_est_ssn,
    randcov_value = randcov_orig2optim_val,
    randcov_is_known = randcov_orig2optim_is_known,
    n_est_rand = n_est_rand,
    n_est = n_est_ssn + n_est_rand,
    classes = c(
      tailup = class(initial_object$tailup_initial),
      taildown = class(initial_object$taildown_initial),
      euclid = class(initial_object$euclid_initial),
      nugget = class(initial_object$nugget_initial)
    )
  )
}

#' Transform parameters from original to optim scale for glms
#'
#' @param initial_object Initial value object
#'
#' @noRd
orig2optim_glm <- function(initial_object) {
  # tailup
  tailup_de <- initial_object$tailup_initial$initial[["de"]]
  tailup_de_log <- log(tailup_de)
  tailup_de_is_known <- initial_object$tailup_initial$is_known[["de"]]

  tailup_range <- initial_object$tailup_initial$initial[["range"]]
  tailup_range_log <- log(tailup_range)
  tailup_range_is_known <- initial_object$tailup_initial$is_known[["range"]]

  # taildown
  taildown_de <- initial_object$taildown_initial$initial[["de"]]
  taildown_de_log <- log(taildown_de)
  taildown_de_is_known <- initial_object$taildown_initial$is_known[["de"]]

  taildown_range <- initial_object$taildown_initial$initial[["range"]]
  taildown_range_log <- log(taildown_range)
  taildown_range_is_known <- initial_object$taildown_initial$is_known[["range"]]

  # euclid
  euclid_de <- initial_object$euclid_initial$initial[["de"]]
  euclid_de_log <- log(euclid_de)
  euclid_de_is_known <- initial_object$euclid_initial$is_known[["de"]]

  euclid_range <- initial_object$euclid_initial$initial[["range"]]
  euclid_range_log <- log(euclid_range)
  euclid_range_is_known <- initial_object$euclid_initial$is_known[["range"]]

  euclid_rotate <- initial_object$euclid_initial$initial[["rotate"]]
  euclid_rotate_prop <- euclid_rotate / pi
  euclid_rotate_logodds <- logit(euclid_rotate_prop)
  euclid_rotate_is_known <- initial_object$euclid_initial$is_known[["rotate"]]

  euclid_scale <- initial_object$euclid_initial$initial[["scale"]]
  euclid_scale_logodds <- logit(euclid_scale)
  euclid_scale_is_known <- initial_object$euclid_initial$is_known[["scale"]]

  # nugget
  nugget <- initial_object$nugget_initial$initial[["nugget"]]
  nugget_log <- log(nugget)
  nugget_is_known <- initial_object$nugget_initial$is_known[["nugget"]]

  # dispersion
  dispersion <- initial_object$dispersion_initial$initial[["dispersion"]]
  dispersion_log <- log(dispersion)
  dispersion_is_known <- initial_object$dispersion_initial$is_known[["dispersion"]]

  # make parameter vectors
  orig2optim_val <- c(
    tailup_de_log = tailup_de_log,
    tailup_range_log = tailup_range_log,
    taildown_de_log = taildown_de_log,
    taildown_range_log = taildown_range_log,
    euclid_de_log = euclid_de_log,
    euclid_range_log = euclid_range_log,
    euclid_rotate_logodds = euclid_rotate_logodds,
    euclid_scale_logodds = euclid_scale_logodds,
    nugget_log = nugget_log,
    dispersion_log = dispersion_log
  )

  orig2optim_is_known <- c(
    tailup_de_is_known = tailup_de_is_known,
    tailup_range_is_known = tailup_range_is_known,
    taildown_de_is_known = taildown_de_is_known,
    taildown_range_is_known = taildown_range_is_known,
    euclid_de_is_known = euclid_de_is_known,
    euclid_range_is_known = euclid_range_is_known,
    euclid_rotate_is_known = euclid_rotate_is_known,
    euclid_scale_is_known = euclid_scale_is_known,
    nugget_is_known = nugget_is_known,
    dispersion_is_known = dispersion_is_known
  )

  n_est_ssn <- sum(!orig2optim_is_known)

  orig2optim_val <- ifelse(orig2optim_val > 50 & !orig2optim_is_known, 50, orig2optim_val)
  orig2optim_val <- ifelse(orig2optim_val < -50 & !orig2optim_is_known, -50, orig2optim_val)


  # random effects
  if (!is.null(initial_object$randcov_initial)) {
    randcov_orig2optim_val <- log(initial_object$randcov_initial$initial)
    names(randcov_orig2optim_val) <- paste(names(initial_object$randcov_initial$initial), "log", sep = "_")
    randcov_orig2optim_is_known <- initial_object$randcov_initial$is_known
    names(randcov_orig2optim_is_known) <- paste(names(initial_object$randcov_initial$is_known), "log", sep = "_")
    n_est_rand <- sum(!randcov_orig2optim_is_known)
    randcov_orig2optim_val <- ifelse(randcov_orig2optim_val > 50 & !randcov_orig2optim_is_known, 50, randcov_orig2optim_val)
    randcov_orig2optim_val <- ifelse(randcov_orig2optim_val < -50 & !randcov_orig2optim_is_known, -50, randcov_orig2optim_val)
  } else {
    randcov_orig2optim_val <- NULL
    randcov_orig2optim_is_known <- NULL
    n_est_rand <- 0
  }


  orig2optim_list <- list(
    value = orig2optim_val,
    is_known = orig2optim_is_known,
    n_est_ssn = n_est_ssn,
    randcov_value = randcov_orig2optim_val,
    randcov_is_known = randcov_orig2optim_is_known,
    n_est_rand = n_est_rand,
    n_est = n_est_ssn + n_est_rand,
    classes = c(
      tailup = class(initial_object$tailup_initial),
      taildown = class(initial_object$taildown_initial),
      euclid = class(initial_object$euclid_initial),
      nugget = class(initial_object$nugget_initial),
      dispersion = class(initial_object$dispersion_initial)
    )
  )
}
