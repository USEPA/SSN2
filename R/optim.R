#' Get optim dotlist
#'
#' Get the dotlist for \code{optim()} within \code{splm()}
#'
#' @param ... Additional arguments to \code{optim()}
#'
#' @return An optim dotlist
#'
#' @noRd
get_optim_dotlist <- function(...) {
  # storing dotlist and setting defaults for optim
  dotlist <- list(...)

  ## l-bfgs-b deafult
  # if (!("method" %in% names(dotlist))) {
  #   dotlist$method <- "L-BFGS-B"
  # }

  # nelder-mead default with lower relative tolerance
  if (!("method" %in% names(dotlist))) {
    dotlist$method <- "Nelder-Mead"
  }

  if (!("hessian" %in% names(dotlist))) {
    dotlist$hessian <- FALSE
  }

  if (!("control" %in% names(dotlist))) {
    dotlist$control <- list()
  }

  if (!("reltol" %in% names(dotlist$control))) {
    dotlist$control$reltol <- 1e-4
  }

  dotlist$lower <- -Inf
  dotlist$upper <- Inf

  # make optim dotlist
  optim_dotlist <- list(gr = NULL, method = dotlist$method, lower = dotlist$lower, upper = dotlist$upper, control = dotlist$control, hessian = dotlist$hessian)
}

#' Get parameters to optimize over in optim (remove known parameters)
#'
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param randcov_orig2optim A \code{randcov_orig2optim} object
#'
#' @return The parameters to optimize over in optim
#'
#' @noRd
get_optim_par <- function(orig2optim_object) {
  if (is.null(orig2optim_object$randcov_value)) {
    par <- orig2optim_object$value[!orig2optim_object$is_known]
  } else {
    ssn_pars <- orig2optim_object$value[!orig2optim_object$is_known]
    randcov_pars <- orig2optim_object$randcov_value[!orig2optim_object$randcov_is_known]
    par <- c(ssn_pars, randcov_pars)
  }
  par
}

#' Fill parameters being optimized with known parameters
#'
#' @param cov_orig2optim A \code{cov_orig2optim} object (parameters on optim scale)
#' @param par The parameters to optimize over in optim
#'
#' @return A covariance parameter vector (on the optim scale)
#'
#' @noRd
fill_optim_par <- function(orig2optim_object, par) {
  orig2optim_object$value[!orig2optim_object$is_known] <- par[seq(1, orig2optim_object$n_est_ssn)]
  par_ssn <- orig2optim_object$value
  if (is.null(orig2optim_object$randcov_value)) {
    par_randcov <- NULL
  } else {
    orig2optim_object$randcov_value[!orig2optim_object$randcov_is_known] <- par[seq(orig2optim_object$n_est_ssn + 1, orig2optim_object$n_est)]
    par_randcov <- orig2optim_object$randcov_value
  }
  list(par_ssn = par_ssn, par_randcov = par_randcov)
}

#' Replace optim method with Brent if only one parameter requires optimization.
#'
#' @param optim_par A vector of parameters to optimize
#' @param optim_dotlist A dotlist with arguments for optim
#'
#' @return An edited method in \code{optim_dotlist} if \code{optim_par} has
#'   length 1. Brent method must have finite lower and upper bounds.
#'
#' @noRd
check_optim_method <- function(optim_par, optim_dotlist) {
  if (length(optim_par) == 1) {
    optim_dotlist$method <- "Brent"
    optim_dotlist$lower <- -50
    optim_dotlist$upper <- 50
  }
  optim_dotlist
}

#' Get parameters to optimize over in optim (remove known parameters) for glms
#'
#' @param spcov_orig2optim A \code{spcov_orig2optim} object
#' @param dispersion_orig2optim A \code{dispersion_orig2optim} object
#' @param randcov_orig2optim A \code{randcov_orig2optim} object
#'
#' @return The parameters to optimize over in optim
#'
#' @noRd
get_optim_par_glm <- function(spcov_orig2optim, dispersion_orig2optim, randcov_orig2optim = NULL) {
  if (is.null(randcov_orig2optim)) {
    par <- spcov_orig2optim$value[!spcov_orig2optim$is_known]
  } else {
    spcov_pars <- spcov_orig2optim$value[!spcov_orig2optim$is_known]
    randcov_pars <- randcov_orig2optim$value[!randcov_orig2optim$is_known]
    par <- c(spcov_pars, randcov_pars)
  }
  dispersion_pars <- dispersion_orig2optim$value[!dispersion_orig2optim$is_known]
  par <- c(par, dispersion_pars)
  par
}
