#' Get covariance parameter output using Laplace log likelihood (for glms)
#'
#' @param initial_object Initial value object
#' @param data_object Data object
#' @param estmethod Estimation method
#' @param optim_dotlist Additional optim arguments
#'
#' @noRd
use_laploglik <- function(initial_object, data_object, estmethod, optim_dotlist) {
  orig2optim_object <- orig2optim_glm(initial_object)

  optim_par <- get_optim_par(orig2optim_object)

  optim_dotlist <- check_optim_method(optim_par, optim_dotlist)

  optim_output <- do.call("optim", c(
    list(
      par = optim_par,
      fn = laploglik,
      orig2optim_object = orig2optim_object,
      data_object = data_object,
      estmethod = estmethod
    ),
    optim_dotlist
  ))


  cov_orig_val <- optim2orig_glm(orig2optim_object, optim_output$par)

  params_object <- get_params_object_glm(orig2optim_object$classes, cov_orig_val)

  if (data_object$anisotropy && !initial_object$euclid_initial$is_known[["rotate"]]) {
    lapll_prods <- laploglik_products(params_object, data_object, estmethod)
    minustwolaploglik <- get_minustwolaploglik(lapll_prods, data_object, estmethod)

    params_object_q2 <- params_object
    params_object_q2$euclid[["rotate"]] <- pi - params_object_q2$euclid[["rotate"]]


    lapll_prods_q2 <- laploglik_products(params_object_q2, data_object, estmethod)


    minustwolaploglik_q2 <- get_minustwolaploglik(lapll_prods_q2, data_object, estmethod)

    ## find appropriate value
    rotate_min <- which.min(c(minustwolaploglik, minustwolaploglik_q2))

    if (rotate_min == 2) {
      params_object <- params_object_q2
    }
  }

  optim_output <- list(
    method = optim_dotlist$method,
    control = optim_dotlist$control, value = optim_output$value,
    counts = optim_output$counts, convergence = optim_output$convergence,
    message = optim_output$message,
    hessian = if (optim_dotlist$hessian) optim_output$hessian else FALSE
  )


  list(
    params_object = params_object, optim_output = optim_output,
    is_known = list(
      tailup = initial_object$tailup_initial$is_known,
      taildown = initial_object$taildown_initial$is_known,
      euclid = initial_object$euclid_initial$is_known,
      nugget = initial_object$nugget_initial$is_known,
      dispersion = initial_object$dispersion$is_known,
      randcov = initial_object$randcov_initial$is_known
    )
  )
}




#' Get covariance parameter output from Laplace log likelihood (for glms) when all parameters are known
#'
#' @param initial_object Initial value object
#' @param data_object Data object
#' @param estmethod Estimation method
#'
#' @noRd
use_laploglik_known <- function(initial_object, data_object, estmethod) {
  params_object <- get_params_object_glm_known(initial_object)

  lapll_prods <- laploglik_products(params_object, data_object, estmethod)


  minustwolaploglik <- get_minustwolaploglik(lapll_prods, data_object, estmethod)

  optim_output <- list(
    method = NA,
    control = NA, value = minustwolaploglik,
    counts = NA, convergence = NA,
    message = NA,
    hessian = NA
  )


  list(
    params_object = params_object, optim_output = optim_output,
    is_known = list(
      tailup = initial_object$tailup_initial$is_known,
      taildown = initial_object$taildown_initial$is_known,
      euclid = initial_object$euclid_initial$is_known,
      nugget = initial_object$nugget_initial$is_known,
      dispersion = initial_object$dispersion$is_known,
      randcov = initial_object$randcov_initial$is_known
    )
  )
}
