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
  orig2optim_object <- orig2optim(initial_object, data_object)

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
    params_object <- resolve_anis_rotation(params_object, data_object, estmethod, is_glm = FALSE)
  }
  params_object <- floor_estimated_nugget(
    params_object, initial_object$nugget_initial$is_known, data_object$diagtol
  )

  # store optim output
  optim_output <- trim_optim_output(optim_output, optim_dotlist)

  # return optim output
  list(
    params_object = params_object, optim_output = optim_output,
    is_known = c(loglik_is_known_base(initial_object), list(randcov = initial_object$randcov_initial$is_known))
  )
}

#' Closed-form REML/ML for the IID (no spatial dependence, no random effects) case
#'
#' When every spatial covariance component is \code{"none"} and there is no
#' random effect, the covariance matrix is \code{nugget * I}: the GLS
#' estimator reduces to ordinary least squares and the REML/ML-optimal
#' nugget has a direct closed form, so no numerical search over the nugget
#' is needed. Only reached when the nugget itself still needs estimating --
#' a fixed (known) nugget in this same setting already has no free
#' parameters at all and goes through \code{\link{use_gloglik_known}()}
#' before this function is ever called.
#'
#' @param initial_object A joint covariance initial-value object (with
#'   grid-search-resolved starting values, though none are used here except
#'   to build the returned \code{is_known} flags).
#' @param data_object A model data object.
#' @param estmethod The estimation method (\code{"reml"} or \code{"ml"}).
#'
#' @return A list with \code{params_object}, \code{optim_output}, and
#'   \code{is_known}, matching \code{\link{use_gloglik}()}/
#'   \code{\link{use_gloglik_known}()}'s return shape.
#'
#' @noRd
use_gloglik_iid <- function(initial_object, data_object, estmethod) {
  X <- do.call(rbind, data_object$X_list)
  n <- data_object$n
  p <- data_object$p

  # data_object$s2 is already the REML-scale OLS residual variance
  # (sum((y - X %*% betahat)^2) / (n - p)), computed once in get_data_object()
  sse <- data_object$s2 * (n - p)
  nugget <- if (estmethod == "reml") data_object$s2 else sse / n

  Xt_X <- crossprod(X, X)
  # log|X'X| via the Cholesky factor's diagonal is numerically more stable
  # than computing det(Xt_X) directly
  logdet_XtX <- 2 * sum(log(diag(chol(Xt_X))))

  gll_prods <- list(l1 = n * log(nugget), l2 = sse / nugget)
  if (estmethod == "reml") {
    gll_prods$l3 <- logdet_XtX - p * log(nugget)
  }
  minustwologlik <- get_minustwologlik(gll_prods, data_object, estmethod)

  params_object <- get_params_object_known(initial_object)
  params_object$nugget[["nugget"]] <- nugget
  params_object <- floor_estimated_nugget(
    params_object, initial_object$nugget_initial$is_known, data_object$diagtol
  )

  optim_output <- known_optim_output_stub(minustwologlik)

  list(
    params_object = params_object, optim_output = optim_output,
    is_known = c(loglik_is_known_base(initial_object), list(randcov = initial_object$randcov_initial$is_known))
  )
}

use_gloglik_known <- function(initial_object, data_object, estmethod) {
  # all parameters known

  # store as a parameter object
  params_object <- get_params_object_known(initial_object)

  # find -2ll
  gll_prods <- gloglik_products(params_object, data_object, estmethod)
  minustwologlik <- get_minustwologlik(gll_prods, data_object, estmethod)

  # mirror optim output
  optim_output <- known_optim_output_stub(minustwologlik)

  # return optim output
  list(
    params_object = params_object, optim_output = optim_output,
    is_known = c(loglik_is_known_base(initial_object), list(randcov = initial_object$randcov_initial$is_known))
  )
}
