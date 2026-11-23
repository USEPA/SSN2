#' Recover data for \code{emmeans} support
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()]
#' @param frame The model frame
#' @param ... Additional arguments passed to \code{emmeans::recover_data()}
#'
#' @return The recovered data, for use by the \code{emmeans} package
#'
#' @details Registered dynamically for the \code{emmeans} generic in
#'   \code{.onLoad()} rather than exported directly, since \code{emmeans}
#'   is only in Suggests.
#'
#' @noRd
recover_data.ssn_lm <- function(object, frame = model.frame(object), ...) {
  # check to see if emmeans installed
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Install the emmeans package before using", call. = FALSE)
  }
  # recover data (using emmeans code)
  fcall = object$call
  # recognize that lm objects have a $model element that is model.frame(object)
  emmeans::recover_data(fcall, delete.response(terms(object)), frame = frame, na.action = NULL, ...)
}

recover_data.ssn_glm <- recover_data.ssn_lm

#' Build the \code{emmeans} basis for a fitted model
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()]
#' @param trms Model terms
#' @param xlev Factor levels
#' @param grid A reference grid
#' @param ... Additional arguments passed to \code{emmeans} helpers
#'
#' @return A list describing the linear model basis (design matrix, coefficients,
#'   covariance matrix, etc.), for use by the \code{emmeans} package
#'
#' @details Registered dynamically for the \code{emmeans} generic in
#'   \code{.onLoad()} rather than exported directly, since \code{emmeans}
#'   is only in Suggests.
#'
#' @noRd
emm_basis.ssn_lm <- function(object, trms, xlev, grid, ...) {
  # emm_basis
  # check to see if emmeans installed
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Install the emmeans package before using", call. = FALSE)
  }
  bhat = coef(object)
  nm = if (is.null(names(bhat)))
    row.names(bhat)
  else names(bhat)
  m = suppressWarnings(model.frame(trms, grid, na.action = na.pass,
                                   xlev = xlev))
  X = model.matrix(trms, m, contrasts.arg = object$contrasts)
  assign = attr(X, "assign")
  # reorder/select columns of X to match the order of the fitted coefficient
  # names so X %*% bhat lines up correctly
  X = X[, nm, drop = FALSE]
  bhat = as.numeric(bhat)
  V = emmeans::.my.vcov(object, ...)
  nbasis = estimability::all.estble # returns a 1x1 NA which says all functions estimable
  misc = list()
  df_setup <- get_emmeans_dffun(object)
  dffun <- df_setup$dffun
  dfargs <- df_setup$dfargs

  mm <- model.matrix(object)
  mm = emmeans::.cmpMM(mm, assign = attr(mm, "assign"))
  if (inherits(object, c("ssn_glm"))) {
    # SSN2 doesn't store the link name directly on glm objects, so look it up
    # from the family (binomial/beta use logit, everything else uses log)
    famdat <- data.frame(
      family = c("poisson", "nbinomial", "binomial", "beta", "Gamma", "inverse.gaussian")
    )
    famdat$link <- ifelse(famdat$family %in% c("binomial", "beta"), "logit", "log")
    fam = famdat[match(object$family, famdat$family), ]
    fam = list(family = fam$family, link = fam$link)
    misc = emmeans::.std.link.labels(fam, misc)
  }
  list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun,
       dfargs = dfargs, misc = misc, model.matrix = mm)
}

emm_basis.ssn_glm <- emm_basis.ssn_lm

#' Build the degrees-of-freedom callback \code{emm_basis.ssn_lm()} passes to \code{emmeans}
#'
#' Satterthwaite denominator degrees of freedom (matching \code{object$ddf})
#' are used for exact \code{ssn_lm()} fits when available; \code{ssn_glm()}
#' fits (and any other case) fall back to the asymptotic (\code{Inf}) degrees
#' of freedom. See [satterthwaite.SSN2()].
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#'
#' @return A list with \code{dffun} (a function of \code{k}/\code{dfargs}
#'   emmeans calls to get a contrast's degrees of freedom, with a
#'   \code{"mesg"} attribute of \code{"asymptotic"} or \code{"satterthwaite"})
#'   and \code{dfargs} (the values that function needs).
#'
#' @noRd
get_emmeans_dffun <- function(object) {
  dfargs <- list()
  dffun <- function(k, dfargs) Inf
  attr(dffun, "mesg") <- "asymptotic"

  if (is.null(object$ddf)) {
    return(list(dffun = dffun, dfargs = dfargs))
  }

  sw <- tryCatch(
    get_satterthwaite_cached(object, "numeric"),
    error = function(e) {
      warning("Satterthwaite degrees of freedom could not be computed for emmeans (", conditionMessage(e), "); falling back to asymptotic (Inf) degrees of freedom.", call. = FALSE)
      NULL
    }
  )

  if (!is.null(sw) && is.null(sw$vcov_theta)) {
    warning("Falling back to asymptotic (Inf) degrees of freedom for emmeans.", call. = FALSE)
    sw <- NULL
  }

  if (!is.null(sw)) {
    # Every function/object dffun needs must be captured as a VALUE in
    # dfargs, not referenced lexically -- emmeans::ref_grid() strips the
    # closure's environment to baseenv().
    dfargs <- list(
      object = object, context = sw$context, vcov_theta = sw$vcov_theta,
      get_df_fn = satterthwaite_df_for_contrast
    )
    dffun <- function(k, dfargs) {
      dfargs$get_df_fn(as.numeric(k), dfargs$object, dfargs$context, dfargs$vcov_theta)
    }
    attr(dffun, "mesg") <- "satterthwaite"
  }

  list(dffun = dffun, dfargs = dfargs)
}
