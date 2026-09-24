#' Get covariance parameter output using Laplace log likelihood (for glms)
#'
#' @param initial_object Initial value object
#' @param data_object Data object
#' @param estmethod Estimation method
#' @param optim_dotlist Additional optim arguments
#'
#' @noRd
use_laploglik <- function(initial_object, data_object, estmethod, optim_dotlist) {
  orig2optim_object <- orig2optim_glm(initial_object, data_object)

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
    params_object <- resolve_anis_rotation(params_object, data_object, estmethod, is_glm = TRUE)
  }
  params_object <- floor_estimated_nugget(
    params_object, initial_object$nugget_initial$is_known, data_object$diagtol
  )

  optim_output <- trim_optim_output(optim_output, optim_dotlist)

  list(
    params_object = params_object, optim_output = optim_output,
    is_known = c(loglik_is_known_base(initial_object), list(
      dispersion = initial_object$dispersion_initial$is_known,
      randcov = initial_object$randcov_initial$is_known
    ))
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

  optim_output <- known_optim_output_stub(minustwolaploglik)

  list(
    params_object = params_object, optim_output = optim_output,
    is_known = c(loglik_is_known_base(initial_object), list(
      dispersion = initial_object$dispersion_initial$is_known,
      randcov = initial_object$randcov_initial$is_known
    ))
  )
}
