cov_estimate_gloglik <- function(data_object, ssn.object, initial_object, estmethod, optim_dotlist) {
  # find the initial NA object
  initial_NA_object <- get_initial_NA_object(initial_object, data_object)

  # grid search to find initial values
  cov_initial_val <- cov_initial_search(initial_NA_object, ssn.object, data_object, estmethod)

  # find the known parameters for each initial NA object
  tailup_is_known <- initial_NA_object$tailup_initial$is_known
  taildown_is_known <- initial_NA_object$taildown_initial$is_known
  euclid_is_known <- initial_NA_object$euclid_initial$is_known
  nugget_is_known <- initial_NA_object$nugget_initial$is_known
  randcov_is_known <- initial_NA_object$randcov_initial$is_known

  # if all parameters are known, find the likelihood; otherwise optimize
  if (all(tailup_is_known, taildown_is_known, euclid_is_known, nugget_is_known, randcov_is_known)) {
    # find the likelihood
    cov_estimate_val <- use_gloglik_known(cov_initial_val$initial_object, data_object, estmethod)
  } else {
    # optimize
    cov_estimate_val <- use_gloglik(cov_initial_val$initial_object, data_object, estmethod, optim_dotlist = optim_dotlist)
  }
  # return value
  cov_estimate_val
}

cov_estimate_laploglik <- function(data_object, ssn.object, initial_object, estmethod, optim_dotlist) {
  initial_NA_object <- get_initial_NA_object_glm(initial_object, data_object)
  cov_initial_val <- cov_initial_search_glm(initial_NA_object, ssn.object, data_object, estmethod)

  # find is known
  tailup_is_known <- initial_NA_object$tailup_initial$is_known
  taildown_is_known <- initial_NA_object$taildown_initial$is_known
  euclid_is_known <- initial_NA_object$euclid_initial$is_known
  nugget_is_known <- initial_NA_object$nugget_initial$is_known
  dispersion_is_known <- initial_NA_object$dispersion_initial$is_known
  randcov_is_known <- initial_NA_object$randcov_initial$is_known




  if (all(tailup_is_known, taildown_is_known, euclid_is_known, nugget_is_known, dispersion_is_known, randcov_is_known)) {
    cov_estimate_val <- use_laploglik_known(cov_initial_val$initial_object, data_object, estmethod)
  } else {
    cov_estimate_val <- use_laploglik(cov_initial_val$initial_object, data_object, estmethod, optim_dotlist = optim_dotlist)
  }
}
