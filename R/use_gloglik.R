#' Get covariance parameter output using log likelihood
#'
#' @param initial_object Initial value object
#' @param data_object Data object
#' @param estmethod Estimation method
#' @param optim_dotlist Additional optim arguments
#'
#' @noRd
use_gloglik <- function(initial_object, data_object, estmethod, optim_dotlist) {
  # take original parameter values and put them on optim scale
  orig2optim_object <- orig2optim(initial_object)

  # find relevant parameters to optimize (don't optimize known parameters)
  optim_par <- get_optim_par(orig2optim_object)

  # optim preliminaries
  optim_dotlist <- check_optim_method(optim_par, optim_dotlist)

  # optimize
  optim_output <- do.call("optim", c(
    list(
      par = optim_par,
      fn = gloglik,
      orig2optim_object = orig2optim_object,
      data_object = data_object,
      estmethod = estmethod
    ),
    optim_dotlist
  ))

  # take optim parameter values and put them on original scale
  cov_orig_val <- optim2orig(orig2optim_object, optim_output$par)

  # store as a parameter object
  params_object <- get_params_object(orig2optim_object$classes, cov_orig_val)

  # find appropriate rotation parameter if anisotropy used
  if (data_object$anisotropy && !initial_object$euclid_initial$is_known[["rotate"]]) {
    gll_prods <- gloglik_products(params_object, data_object, estmethod)
    minustwologlik <- get_minustwologlik(gll_prods, data_object, estmethod)

    params_object_q2 <- params_object
    params_object_q2$euclid[["rotate"]] <- pi - params_object_q2$euclid[["rotate"]]


    gll_prods_q2 <- gloglik_products(params_object_q2, data_object, estmethod)


    minustwologlik_q2 <- get_minustwologlik(gll_prods_q2, data_object, estmethod)

    ## find appropriate value
    rotate_min <- which.min(c(minustwologlik, minustwologlik_q2))

    if (rotate_min == 2) {
      params_object <- params_object_q2
    }
  }

  # store optim output
  optim_output <- list(
    method = optim_dotlist$method,
    control = optim_dotlist$control, value = optim_output$value,
    counts = optim_output$counts, convergence = optim_output$convergence,
    message = optim_output$message,
    hessian = if (optim_dotlist$hessian) optim_output$hessian else FALSE
  )


  # return optim output
  list(
    params_object = params_object, optim_output = optim_output,
    is_known = list(
      tailup = initial_object$tailup_initial$is_known,
      taildown = initial_object$taildown_initial$is_known,
      euclid = initial_object$euclid_initial$is_known,
      nugget = initial_object$nugget_initial$is_known,
      randcov = initial_object$randcov_initial$is_known
    )
  )
}



#' Get covariance parameter output from log likelihood when all parameters are known
#'
#' @param initial_object Initial value object
#' @param data_object Data object
#' @param estmethod Estimation method
#'
#' @noRd
use_gloglik_known <- function(initial_object, data_object, estmethod) {
  # all parameters known

  # store as a parameter object
  params_object <- get_params_object_known(initial_object)

  # find -2ll
  gll_prods <- gloglik_products(params_object, data_object, estmethod)
  minustwologlik <- get_minustwologlik(gll_prods, data_object, estmethod)

  # mirror optim output
  optim_output <- list(
    method = NA,
    control = NA, value = minustwologlik,
    counts = NA, convergence = NA,
    message = NA,
    hessian = NA
  )

  # return optim output
  list(
    params_object = params_object, optim_output = optim_output,
    is_known = list(
      tailup = initial_object$tailup_initial$is_known,
      taildown = initial_object$taildown_initial$is_known,
      euclid = initial_object$euclid_initial$is_known,
      nugget = initial_object$nugget_initial$is_known,
      randcov = initial_object$randcov_initial$is_known
    )
  )
}
