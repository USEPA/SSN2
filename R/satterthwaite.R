#' Compute Satterthwaite denominator degrees of freedom
#'
#' @description Compute Satterthwaite denominator degrees of freedom
#'   \eqn{t}-based (rather than asymptotic \eqn{z}-based)
#'   fixed effect inference in small samples.
#'
#' @param object A fitted model object from [ssn_lm()] fit without \code{local}.
#'   Not supported for [ssn_glm()] objects.
#' @param method The method by which to compute gradients. Must be supplied
#'   as \code{"numeric"} for numerical differentiation, the only supported method.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Satterthwaite degrees of freedom are generally more appropriate than
#'   asymptotic degrees of freedom for small samples. They can be computationally costly
#'   for sample sizes exceeding 500; however, for sample sizes this large, they Satterthwaite
#'   and asymptotic degrees of freedom should yield very similar inferences.
#'
#' @return A named numeric vector of Satterthwaite degrees of freedom for each
#'   fixed effect.
#'
#' @seealso [ssn_lm()] [anova.SSN2()]
#'
#' @name satterthwaite.SSN2
#' @method satterthwaite ssn_lm
#' @order 1
#' @export
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' satterthwaite(ssn_mod, method = "numeric")
#'
#' @references
#'   Rencher, Alvin C. and Schaalje, G. Bruce (2008). Linear Models in
#'   Statistics, Second Edition. John Wiley & Sons.
satterthwaite.ssn_lm <- function(object, method, ...) {
  # Reuses object$vcov$cov (the numerical Hessian/Jacobian result) when the
  # model was fit with ssn_lm(..., ddf = "satterthwaite"), skipping the
  # expensive recomputation; always rebuilds the (cheap) context fresh.
  sw <- get_satterthwaite_cached(object, method)
  context <- sw$context
  vcov_theta <- sw$vcov_theta

  p <- object$p
  coef_names <- rownames(object$vcov$fixed)
  df <- stats::setNames(rep(NA_real_, p), coef_names)

  if (is.null(vcov_theta)) {
    return(df)
  }

  ident <- diag(p)
  for (k in seq_len(p)) {
    Li <- as.numeric(ident[k, ])
    df[k] <- satterthwaite_df_for_contrast(
      Li, object, context, vcov_theta,
      contrast_label = coef_names[k]
    )
  }

  df
}

#' @rdname satterthwaite.SSN2
#' @method satterthwaite ssn_glm
#' @order 2
#' @export
satterthwaite.ssn_glm <- function(object, ...) {
  stop("Satterthwaite degrees of freedom are only implemented for exact ssn_lm() models fit without 'local'; they are not supported for ssn_glm() models.", call. = FALSE)
}

#' Compute Satterthwaite degrees of freedom for one linear contrast
#'
#' Computes the Satterthwaite approximation for a linear contrast:
#' \eqn{df = 2 g(\hat\theta)^2 / (\nabla g(\hat\theta)' \mathrm{Cov}(\hat\theta)
#' \nabla g(\hat\theta))}, where \eqn{g(\theta) = L_i' V_\beta(\theta) L_i}.
#' Returns \code{NA} (with a warning) if the denominator is non-positive or
#' non-finite (an aliased or otherwise non-estimable contrast).
#'
#' @param Li A single linear contrast vector (length \code{p}).
#' @param object A fitted \code{ssn_lm} model object.
#' @param context A Satterthwaite context from
#'   \code{\link{get_satterthwaite_context}()}.
#' @param vcov_theta The numeric variance-covariance matrix of the free
#'   covariance parameters.
#' @param contrast_label An optional label for the contrast, used in the
#'   degeneracy warning message.
#'
#' @return The contrast's Satterthwaite degrees of freedom, or \code{NA_real_}
#'   if degenerate.
#'
#' @noRd
satterthwaite_df_for_contrast <- function(Li, object, context, vcov_theta, contrast_label = NULL) {
  g <- as.numeric(crossprod(Li, object$vcov$fixed %*% Li))

  grad_g <- suppressWarnings(numDeriv::grad(
    function(theta_free) get_grad_gi(theta_free, Li, context),
    context$cov_val_free
  ))

  denom <- as.numeric(crossprod(grad_g, vcov_theta) %*% grad_g)

  if (!is.finite(denom) || denom <= 0) {
    label <- if (is.null(contrast_label)) "a contrast" else contrast_label
    warning("Satterthwaite degrees of freedom for \"", label, "\" are degenerate (the contrast may be aliased or non-estimable); returning NA.", call. = FALSE)
    return(NA_real_)
  }

  2 * g^2 / denom
}

#' Compute Fai-Cornelius approximate denominator degrees of freedom for a joint contrast
#'
#' Diagonalizes the joint contrast's covariance, computes each eigen-contrast's
#' own Satterthwaite degrees of freedom via
#' \code{\link{satterthwaite_df_for_contrast}()}, and combines them via the
#' Fai-Cornelius approximation used for multi-row (joint, e.g. \code{anova()})
#' contrasts.
#'
#' @param L A contrast matrix (or vector, coerced to a one-row matrix).
#' @param object A fitted \code{ssn_lm} model object.
#' @param context A Satterthwaite context from
#'   \code{\link{get_satterthwaite_context}()}.
#' @param vcov_theta The numeric variance-covariance matrix of the free
#'   covariance parameters, or \code{NULL}.
#'
#' @return The joint contrast's approximate denominator degrees of freedom,
#'   or \code{NA_real_} if \code{vcov_theta} is \code{NULL} or the
#'   approximation is degenerate.
#'
#' @noRd
fai_cornelius <- function(L, object, context, vcov_theta) {
  if (!is.matrix(L)) {
    L <- matrix(L, nrow = 1)
  }
  q <- nrow(L)

  if (is.null(vcov_theta)) {
    return(NA_real_)
  }

  C <- L %*% tcrossprod(object$vcov$fixed, L)
  C <- as.matrix(Matrix::forceSymmetric(C))
  eig <- eigen(C, symmetric = TRUE)
  U <- eig$vectors
  eta <- crossprod(U, L)

  nu <- vapply(seq_len(q), function(i) {
    satterthwaite_df_for_contrast(
      eta[i, ], object, context, vcov_theta,
      contrast_label = paste0("joint contrast row ", i)
    )
  }, numeric(1))

  if (any(!is.finite(nu)) || any(nu <= 2)) {
    return(NA_real_)
  }

  E <- sum(nu / (nu - 2))
  if (E <= q) {
    return(NA_real_)
  }

  2 * E / (E - q)
}
