#' Reject unsupported model/method combinations for Satterthwaite degrees of freedom
#'
#' Validates: \code{object} is an \code{ssn_lm} (not \code{ssn_glm}) fit --
#' the delta-method/numerical-Hessian derivations assume a Gaussian
#' likelihood, so \code{ssn_glm()} fits have no analog; \code{object} was fit
#' without \code{local} -- local fitting partitions the data and
#' approximates the covariance matrix piecewise, so there is no single
#' well-defined likelihood/covariance matrix left to differentiate globally;
#' \code{object$estmethod} is \code{"reml"} or \code{"ml"} -- betahat/theta_hat
#' must come from (RE)ML for these derivations to apply; and \code{method} is
#' exactly \code{"numeric"} -- the only method currently implemented (no
#' closed-form derivatives, unlike spmodel's \code{splm()}/\code{spautor()}).
#'
#' @param object A fitted model object.
#' @param method The requested Satterthwaite method.
#'
#' @return \code{NULL}, invisibly, if every check passes; otherwise an error.
#'
#' @noRd
validate_satterthwaite_scope <- function(object, method) {
  if (inherits(object, "ssn_glm")) {
    stop("Satterthwaite degrees of freedom are only implemented for exact ssn_lm() models fit without 'local'; they are not supported for ssn_glm() models.", call. = FALSE)
  }
  if (!inherits(object, "ssn_lm")) {
    stop("Satterthwaite degrees of freedom are only implemented for ssn_lm() model objects.", call. = FALSE)
  }
  if (!is.null(object$local_index)) {
    stop("Satterthwaite degrees of freedom can only be computed for models fit without 'local'.", call. = FALSE)
  }
  if (!object$estmethod %in% c("reml", "ml")) {
    stop("Satterthwaite degrees of freedom are only defined for estmethod \"reml\" or \"ml\".", call. = FALSE)
  }
  if (missing(method) || is.null(method) || length(method) != 1 || !identical(method, "numeric")) {
    stop("Satterthwaite degrees of freedom currently only support method = \"numeric\". No closed-form or automatic method selection is implemented.", call. = FALSE)
  }
  invisible(NULL)
}

#' Build the context used to differentiate the covariance matrix with respect to its free parameters
#'
#' Bundles every quantity needed to treat the fitted covariance matrix
#' \eqn{\Sigma} as a function of its free parameters \eqn{\theta}
#' (tailup/taildown/euclid/nugget/random-effect variance components), by
#' pinning the fitted (known plus estimated) values as the initial/known
#' values of a fresh initial object. Built once per \code{satterthwaite()}
#' call and reused for every coefficient/contrast, since it does not depend
#' on which one is being tested.
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param method The requested Satterthwaite method (validated via
#'   \code{\link{validate_satterthwaite_scope}()}).
#'
#' @return A list with \code{orig2optim_object}, \code{eta} (the fitted
#'   optimization-scale free parameters), \code{cov_names_free}/
#'   \code{cov_names_free_orig} (optimization- and original-scale free
#'   parameter names), \code{cov_val_free} (their fitted original-scale
#'   values), \code{data_object}, \code{estmethod}, \code{params_object},
#'   \code{is_known}, and \code{object}.
#'
#' @noRd
get_satterthwaite_context <- function(object, method) {
  validate_satterthwaite_scope(object, method)

  params_object <- object$coefficients$params_object
  is_known <- object$is_known

  tailup_type <- remove_covtype(class(params_object$tailup))
  taildown_type <- remove_covtype(class(params_object$taildown))
  euclid_type <- remove_covtype(class(params_object$euclid))
  nugget_type <- remove_covtype(class(params_object$nugget))

  pinned_tailup <- tailup_initial(tailup_type,
    de = params_object$tailup[["de"]], range = params_object$tailup[["range"]],
    known = names(which(is_known$tailup))
  )
  pinned_taildown <- taildown_initial(taildown_type,
    de = params_object$taildown[["de"]], range = params_object$taildown[["range"]],
    known = names(which(is_known$taildown))
  )
  if (euclid_has_extra(euclid_type)) {
    pinned_euclid <- euclid_initial(euclid_type,
      de = params_object$euclid[["de"]], range = params_object$euclid[["range"]],
      rotate = params_object$euclid[["rotate"]], scale = params_object$euclid[["scale"]],
      known = names(which(is_known$euclid)), extra = params_object$euclid[["extra"]]
    )
  } else {
    pinned_euclid <- euclid_initial(euclid_type,
      de = params_object$euclid[["de"]], range = params_object$euclid[["range"]],
      rotate = params_object$euclid[["rotate"]], scale = params_object$euclid[["scale"]],
      known = names(which(is_known$euclid))
    )
  }
  pinned_nugget <- nugget_initial(nugget_type,
    nugget = params_object$nugget[["nugget"]],
    known = names(which(is_known$nugget))
  )

  pinned_initial_object <- get_initial_object(
    tailup_type = tailup_type, taildown_type = taildown_type,
    euclid_type = euclid_type, nugget_type = nugget_type,
    tailup_initial = pinned_tailup, taildown_initial = pinned_taildown,
    euclid_initial = pinned_euclid, nugget_initial = pinned_nugget
  )

  if (!is.null(object$random)) {
    pinned_randcov_initial <- list(initial = params_object$randcov, is_known = is_known$randcov)
    pinned_initial_object$randcov_initial <- pinned_randcov_initial
  } else {
    pinned_randcov_initial <- NULL
  }

  # rebuilt fresh (not reused from the fit) with range_constrain forced off:
  # constraining is purely an optimizer-search-space device, and once a point
  # estimate is pinned here, the plain log-scale Jacobian is what the delta
  # method below wants -- logit-odds-scale gradients near a boundary are
  # worse-conditioned for no benefit (matches spmodel's own gradient-context
  # rebuild, which does the same).
  data_object <- get_data_object(
    formula = object$formula, ssn.object = object$ssn.object, additive = object$additive,
    anisotropy = object$anisotropy, initial_object = pinned_initial_object,
    random = object$random, randcov_initial = pinned_randcov_initial,
    partition_factor = object$partition_factor, local = NULL, range_constrain = FALSE
  )

  orig2optim_object <- orig2optim(pinned_initial_object, data_object)
  eta <- get_optim_par(orig2optim_object)
  cov_names_free <- names(eta)

  if (length(cov_names_free) == 0) {
    stop("All covariance parameters are known; Satterthwaite degrees of freedom are not applicable.", call. = FALSE)
  }

  cov_names_free_orig <- sub("_logodds$", "", sub("_log$", "", cov_names_free))

  euclid_orig_named <- c(
    euclid_de = params_object$euclid[["de"]], euclid_range = params_object$euclid[["range"]]
  )
  if (euclid_has_extra(euclid_type)) {
    euclid_orig_named <- c(euclid_orig_named, euclid_extra = params_object$euclid[["extra"]])
  }
  full_orig_named <- c(
    tailup_de = params_object$tailup[["de"]], tailup_range = params_object$tailup[["range"]],
    taildown_de = params_object$taildown[["de"]], taildown_range = params_object$taildown[["range"]],
    euclid_orig_named,
    euclid_rotate = params_object$euclid[["rotate"]], euclid_scale = params_object$euclid[["scale"]],
    nugget = params_object$nugget[["nugget"]]
  )
  if (!is.null(params_object$randcov)) {
    full_orig_named <- c(full_orig_named, params_object$randcov)
  }
  cov_val_free <- full_orig_named[cov_names_free_orig]
  names(cov_val_free) <- cov_names_free_orig

  list(
    orig2optim_object = orig2optim_object,
    eta = eta,
    cov_names_free = cov_names_free,
    cov_names_free_orig = cov_names_free_orig,
    cov_val_free = cov_val_free,
    data_object = data_object,
    estmethod = object$estmethod,
    params_object = params_object,
    is_known = is_known,
    object = object
  )
}

#' Get (reusing a fit-time cache when available) the covariance-parameter covariance matrix
#'
#' Reuses \code{object$vcov$cov} (the numerical Hessian/Jacobian result
#' cached at fit time when \code{ssn_lm(..., ddf = "satterthwaite")}
#' succeeded) instead of paying for the expensive
#' \eqn{\mathrm{Cov}(\hat\theta)} computation a second time; the (cheap)
#' context is always rebuilt fresh since it is not itself stored on the
#' object.
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param method The requested Satterthwaite method.
#'
#' @return A list with \code{context} (from
#'   \code{\link{get_satterthwaite_context}()}) and \code{vcov_theta} (the
#'   covariance-parameter covariance matrix, or \code{NULL} if it could not
#'   be computed).
#'
#' @noRd
get_satterthwaite_cached <- function(object, method) {
  context <- get_satterthwaite_context(object, method)

  if (!is.null(object$ddf)) {
    return(list(context = context, vcov_theta = object$vcov$cov))
  }

  vt <- get_vcov_theta_numeric(context)
  list(context = context, vcov_theta = vt$vcov_theta)
}

#' Determine a \code{ddf} argument's method
#'
#' Determine whether to use satterthwaite or asymptotic degrees of freedom.
#'
#' @param ddf The \code{ddf} argument as passed by the caller (already
#'   resolved from \code{missing()} to \code{NULL})
#' @param n The fitted model's sample size
#'
#' @return Either \code{"asymptotic"} or \code{"satterthwaite"}
#'
#' @noRd
determine_ddf <- function(ddf, n) {
  if (is.null(ddf)) ddf <- if (n <= 500) "satterthwaite" else "asymptotic"
  if (!is.character(ddf) || length(ddf) != 1L || is.na(ddf) ||
      !ddf %in% c("asymptotic", "satterthwaite")) {
    stop("ddf must be \"asymptotic\" or \"satterthwaite\".", call. = FALSE)
  }
  ddf
}

#' Compute (or skip) denominator degrees of freedom for a fitted model
#'
#' Implements \code{ssn_lm()}'s \code{ddf} argument (see \code{determine_ddf()}
#' for how it resolves the sample size to \code{"asymptotic"} or
#' \code{"satterthwaite"}): \code{"asymptotic"} returns every element
#' \code{NULL} (unchanged from before this argument existed);
#' \code{"satterthwaite"} attempts \code{compute_satterthwaite_fit_time()} on
#' the already-fitted object and returns everything it produces, or every
#' element \code{NULL} if the computation errors, produces invalid degrees of
#' freedom, or a non-positive-definite covariance matrix -- a failure here
#' must never prevent \code{ssn_lm()} from returning a fitted model. An
#' explicit \code{ddf} request (as opposed to the size-based default) still
#' surfaces a diagnostic warning on failure. \code{vcov_cov} (the estimated
#' covariance matrix of the free covariance parameters) is a byproduct of
#' computing \code{ddf} already needed here, so both are returned together
#' rather than computing \code{vcov_theta} a second time via
#' \code{satterthwaite()} -- see [vcov.SSN2()], which exposes it as
#' \code{vcov(object, type = "cov")}.
#'
#' @param object The fitted model object.
#' @param ddf The \code{ddf} method.
#'
#' @return A list with elements \code{ddf} (a named numeric vector of
#'   denominator df, see [satterthwaite.SSN2()], or \code{NULL}) and
#'   \code{vcov_cov} (the corresponding covariance-parameter covariance
#'   matrix, or \code{NULL})
#'
#' @noRd
get_fit_ddf <- function(object, ddf) {
  none <- list(ddf = NULL, vcov_cov = NULL)
  if (determine_ddf(ddf, object$n) == "asymptotic") return(none)
  out <- withCallingHandlers(
    tryCatch(compute_satterthwaite_fit_time(object), error = function(e) {
      if (!is.null(ddf)) {
        warning("Satterthwaite degrees of freedom could not be computed (", conditionMessage(e),
                "); using asymptotic inference.", call. = FALSE)
      }
      NULL
    }),
    warning = function(w) {
      if (is.null(ddf) && grepl("not numerically positive definite", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (is.null(out) || is.null(out$vcov_cov) ||
      !is.numeric(out$ddf) || length(out$ddf) != object$p ||
      any(!is.finite(out$ddf) | out$ddf <= 0)) return(none)
  out
}

#' Compute Satterthwaite degrees of freedom for every coefficient at fit time
#'
#' Stores covariance-parameter uncertainty in \code{object$vcov$cov} so
#' subsequent inference can reuse it without recomputing the Hessian.
#'
#' @param object A fitted \code{ssn_lm} model object.
#'
#' @return A list with \code{ddf} (a named numeric vector, one entry per
#'   fixed-effect coefficient) and \code{vcov_cov} (the covariance-parameter
#'   covariance matrix, or \code{NULL} if it could not be computed).
#'
#' @noRd
compute_satterthwaite_fit_time <- function(object) {
  context <- get_satterthwaite_context(object, "numeric")
  vt <- get_vcov_theta_numeric(context)

  p <- object$p
  coef_names <- rownames(object$vcov$fixed)
  df <- stats::setNames(rep(NA_real_, p), coef_names)

  if (is.null(vt$vcov_theta)) {
    return(list(ddf = df, vcov_cov = NULL))
  }

  ident <- diag(p)
  for (k in seq_len(p)) {
    Li <- as.numeric(ident[k, ])
    df[k] <- satterthwaite_df_for_contrast(
      Li, object, context, vt$vcov_theta,
      contrast_label = coef_names[k]
    )
  }

  list(ddf = df, vcov_cov = vt$vcov_theta)
}

#' Numerically compute the sampling covariance matrix of the free covariance parameters
#'
#' Computes \eqn{\mathrm{Cov}(\hat\theta)}, used by
#' \code{\link{satterthwaite_df_for_contrast}()}: standard MLE asymptotics,
#' \eqn{\mathrm{Cov}(\hat\eta) \approx [\mathrm{Fisher\ info}]^{-1}},
#' approximated by numerically differentiating the fitting log-likelihood
#' itself at the fitted optimizer-scale value (\code{eta}), then delta-method
#' mapped back to the original (\eqn{\theta}) parameterization. Unlike
#' spmodel's \code{splm()}/\code{spautor()}, no closed-form alternative is
#' implemented, so this numeric approach is always used.
#'
#' @param context A Satterthwaite context from
#'   \code{\link{get_satterthwaite_context}()}.
#' @param numderiv_args Optional \code{method.args} passed to
#'   \code{numDeriv::hessian()}/\code{numDeriv::jacobian()}.
#'
#' @return A list with \code{vcov_theta} (the covariance matrix on the
#'   original parameter scale, or \code{NULL} if it is not numerically
#'   positive definite), \code{vcov_eta} (the optimizer-scale covariance
#'   matrix), \code{J} (the delta-method Jacobian), and \code{H_eta} (the
#'   optimizer-scale Hessian of \eqn{-2} times the log-likelihood).
#'
#' @noRd
get_vcov_theta_numeric <- function(context, numderiv_args = NULL) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Install the numDeriv package before using method = \"numeric\" for Satterthwaite degrees of freedom.", call. = FALSE)
  }

  obj_eta <- function(eta) {
    gloglik(eta, context$orig2optim_object, context$data_object, context$estmethod)
  }

  # obj_eta() returns -2*log-likelihood (see gloglik()), so its Hessian is
  # -2 * Hessian(log-lik); dividing by 2 converts that into the observed
  # Fisher information (the negative log-lik Hessian), whose inverse is the
  # usual MLE asymptotic covariance estimate.
  #
  # numDeriv perturbs the covariance parameters at many nearby points to
  # build this Hessian via finite differences, and those perturbations
  # routinely wander into numerically awkward (but perfectly legitimate)
  # regions -- e.g. matern's smoothness parameter pushing base R's besselK()
  # to warn "value out of range in 'bessel_k'", or a near-zero variance
  # component briefly going slightly negative. These warnings are incidental
  # to the finite-differencing process itself, not a sign that the resulting
  # Cov(theta_hat) is wrong (checked separately via the positive-definiteness
  # check just below), and a single ssn_lm() fit can otherwise emit thousands
  # of them. Genuine failures still surface: this suppression only affects
  # warning(), not error(), and stop()'d likelihood evaluations propagate
  # normally.
  hessian_args <- list(func = obj_eta, x = context$eta)
  if (!is.null(numderiv_args)) hessian_args$method.args <- numderiv_args
  H_eta <- suppressWarnings(do.call(numDeriv::hessian, hessian_args))
  H_eta <- as.matrix(Matrix::forceSymmetric(H_eta))

  vcov_eta <- tryCatch(chol2inv(chol(H_eta / 2)), error = function(e) NULL)
  if (is.null(vcov_eta)) {
    warning("The covariance matrix of the covariance parameters is not numerically positive definite.", call. = FALSE)
    return(list(vcov_theta = NULL, vcov_eta = NULL, J = NULL, H_eta = H_eta))
  }

  # theta(eta): unpack the optim-scale vector back to the original scale and
  # keep only the free entries -- computed via a full numerical Jacobian
  # because it accommodates any transformation (e.g., exp vs expit)
  delta_transform_eta <- function(eta) {
    fill <- fill_optim_par(context$orig2optim_object, eta)
    euclid_type <- remove_covtype(context$orig2optim_object$classes[["euclid"]])
    orig_ssn <- optim2orig_ssn_components(
      fill$par_ssn, euclid_type, context$orig2optim_object$range_constrain_value
    )
    orig_randcov <- optim2orig_randcov_components(fill$par_randcov)
    full <- c(orig_ssn, orig_randcov)
    unname(full[context$cov_names_free_orig])
  }

  # delta method: Cov(theta_hat) ~= J Cov(eta_hat) J', where J is the
  # Jacobian of the (nonlinear, e.g. log/logit) optim-scale-to-original-scale
  # transform evaluated at the fitted value -- maps the Hessian-based
  # covariance above off of the unconstrained optimizer scale and onto the
  # original, interpretable covariance parameter scale used everywhere else
  #
  # suppressWarnings(): see the H_eta call above -- same finite-difference
  # perturbation, same incidental warnings
  jacobian_args <- list(func = delta_transform_eta, x = context$eta)
  if (!is.null(numderiv_args)) jacobian_args$method.args <- numderiv_args
  J <- suppressWarnings(do.call(numDeriv::jacobian, jacobian_args))

  vcov_theta <- J %*% tcrossprod(vcov_eta, J)
  dimnames(vcov_theta) <- list(context$cov_names_free_orig, context$cov_names_free_orig)

  list(vcov_theta = vcov_theta, vcov_eta = vcov_eta, J = J, H_eta = H_eta)
}

#' List every spatial covariance field name in canonical order
#'
#' @return A character vector of the 10 possible spatial covariance field
#'   names (\code{tailup_de}, \code{tailup_range}, \code{taildown_de},
#'   \code{taildown_range}, \code{euclid_de}, \code{euclid_range},
#'   \code{euclid_extra}, \code{euclid_rotate}, \code{euclid_scale},
#'   \code{nugget}).
#'
#' @noRd
get_spcov_field_names <- function() {
  c(
    "tailup_de", "tailup_range", "taildown_de", "taildown_range",
    "euclid_de", "euclid_range", "euclid_extra", "euclid_rotate", "euclid_scale", "nugget"
  )
}

#' Build a covariance parameter object with free fields perturbed to given values
#'
#' Starts from \code{context$params_object} (the fitted values) and
#' overwrites only the free (not-known) fields named in \code{cov_val_free},
#' used to evaluate the log-likelihood/fixed-effect covariance at nearby
#' covariance-parameter values during numerical differentiation.
#'
#' @param cov_val_free A named numeric vector of free covariance-parameter
#'   values (names matching \code{context$cov_names_free_orig}).
#' @param context A Satterthwaite context from
#'   \code{\link{get_satterthwaite_context}()}.
#'
#' @return A joint covariance parameter object with free fields set to
#'   \code{cov_val_free} and known fields unchanged.
#'
#' @noRd
fill_perturbed_params_object <- function(cov_val_free, context) {
  params_object <- context$params_object
  cov_val_free <- as.numeric(cov_val_free)
  names(cov_val_free) <- context$cov_names_free_orig

  set_if_present <- function(full_name) {
    if (full_name %in% names(cov_val_free)) unname(cov_val_free[[full_name]]) else NULL
  }

  val <- set_if_present("tailup_de")
  if (!is.null(val)) params_object$tailup[["de"]] <- val
  val <- set_if_present("tailup_range")
  if (!is.null(val)) params_object$tailup[["range"]] <- val

  val <- set_if_present("taildown_de")
  if (!is.null(val)) params_object$taildown[["de"]] <- val
  val <- set_if_present("taildown_range")
  if (!is.null(val)) params_object$taildown[["range"]] <- val

  val <- set_if_present("euclid_de")
  if (!is.null(val)) params_object$euclid[["de"]] <- val
  val <- set_if_present("euclid_range")
  if (!is.null(val)) params_object$euclid[["range"]] <- val
  val <- set_if_present("euclid_extra")
  if (!is.null(val)) params_object$euclid[["extra"]] <- val
  val <- set_if_present("euclid_rotate")
  if (!is.null(val)) params_object$euclid[["rotate"]] <- val
  val <- set_if_present("euclid_scale")
  if (!is.null(val)) params_object$euclid[["scale"]] <- val

  val <- set_if_present("nugget")
  if (!is.null(val)) params_object$nugget[["nugget"]] <- val

  if (!is.null(params_object$randcov)) {
    randcov_free <- intersect(names(params_object$randcov), names(cov_val_free))
    if (length(randcov_free) > 0) {
      params_object$randcov[randcov_free] <- cov_val_free[randcov_free]
    }
  }

  params_object
}

#' Evaluate \eqn{g(\theta) = L_i' V_\beta(\theta) L_i} at a covariance-parameter value
#'
#' The function numerically differentiated by
#' \code{\link{satterthwaite_df_for_contrast}()} (via \code{numDeriv::grad()}):
#' rebuilds the fixed-effect covariance matrix at a (possibly perturbed)
#' covariance-parameter value and evaluates one linear contrast's variance
#' under it.
#'
#' @param cov_val_free A named numeric vector of free covariance-parameter
#'   values to evaluate at.
#' @param Li A single linear contrast vector (length \code{p}).
#' @param context A Satterthwaite context from
#'   \code{\link{get_satterthwaite_context}()}.
#'
#' @return The scalar \eqn{g(\theta) = L_i' V_\beta(\theta) L_i}.
#'
#' @noRd
get_grad_gi <- function(cov_val_free, Li, context) {
  params_object <- fill_perturbed_params_object(cov_val_free, context)

  cov_matrix_list <- get_cov_matrix_list(params_object, context$data_object)
  eigenprods <- get_eigenprods(
    cov_matrix_list[[1]], context$data_object$X_list[[1]],
    context$data_object$y_list[[1]], context$data_object$ones_list[[1]]
  )
  invcov_betahat <- crossprod(eigenprods$SqrtSigInv_X)
  cov_betahat_theta <- chol2inv(chol(Matrix::forceSymmetric(invcov_betahat)))

  as.numeric(crossprod(Li, as.matrix(cov_betahat_theta) %*% Li))
}
