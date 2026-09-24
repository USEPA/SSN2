#' Perform leave-one-out cross validation
#'
#' @description Perform leave-one-out cross validation with options for computationally
#'   efficient approximations for big data.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param cv_predict A logical indicating whether the leave-one-out fitted values
#'   should be returned. Defaults to \code{FALSE}. If \code{object} is from [ssn_glm()],
#'   the fitted values returned are on the link scale.
#' @param se.fit A logical indicating whether the leave-one-out
#'   prediction standard errors should be returned. Defaults to \code{FALSE}.
#'   If \code{object} is from [ssn_glm()],
#'   the standard errors correspond to the fitted values returned on the link scale.
#' @param local A list or logical. If a list, specific list elements described
#'   in [predict.SSN2()] control the big data approximation behavior.
#'   If a logical, \code{TRUE} chooses default list elements for the list version
#'   of \code{local} as specified in [predict.SSN2()]. Defaults to \code{FALSE},
#'   which performs exact computations.
#' @param interval Whether to report coverage at a requested prediction interval
#'   level. \code{"none"} preserves the legacy coverage columns; \code{"prediction"}
#'   additionally reports coverage at \code{level}. Only available for
#'   \code{ssn_lm()} objects.
#' @param level The requested prediction interval coverage level. Defaults to
#'   \code{0.95}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Each observation is held-out from the data set and the remaining data
#'   are used to make a prediction for the held-out observation. This is compared
#'   to the true value of the observation and several fit statistics are (sometimes optionally) computed:
#'   bias, mean-squared-prediction error (MSPE), root-mean-squared-prediction
#'   error (RMSPE), and the squared correlation (cor2) between the observed data
#'   and leave-one-out predictions (regarded as a prediction version of r-squared
#'   appropriate for comparing across spatial and nonspatial models), and 
#'   prediction interval coverage (cover.XX). Generally,
#'   bias should be near zero and prediction interval coverage at the
#'   intended level for well-fitting models. The lower the MSPE and RMSPE,
#'   the better the model fit (according to the leave-out-out criterion).
#'   The higher the cor2, the better the model fit (according to the leave-out-out
#'   criterion). cor2 and cover.XX are not returned when \code{object} was fit using
#'   \code{ssn_glm()} because we do not observe the underlying latent mean.
#'
#' @return If \code{cv_predict = FALSE} and \code{se.fit = FALSE},
#'   a fit statistics tibble (with bias, MSPE, RMSPE, and cor2; see Details).
#'   If \code{cv_predict = TRUE} or \code{se.fit = TRUE},
#'   a list with elements: \code{stats}, a fit statistics tibble
#'   (with bias, MSPE, RMSPE, and cor2; see Details); \code{cv_predict}, a numeric vector
#'   with leave-one-out predictions for each observation (if \code{cv_predict = TRUE});
#'   and \code{se.fit}, a numeric vector with leave-one-out prediction standard
#'   errors for each observation (if \code{se.fit = TRUE}). When \code{object} is from
#'   \code{splm()} or \code{spautor()} and \code{interval = "prediction"}, the fit
#'   statistics tibble also has \code{cover.XX} columns (e.g. \code{cover.95}
#'   for \code{level = 0.95}; see Details).
#'
#' @name loocv.SSN2
#' @method loocv ssn_lm
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
#' loocv(ssn_mod)
loocv.ssn_lm <- function(object, cv_predict = FALSE, se.fit = FALSE,
                         local, interval = c("none", "prediction"), level = 0.95, ...) {
  if (missing(local)) local <- NULL
  interval <- match.arg(interval)
  local <- resolve_cv_auto_local(local, object$n, "loocv")
  local <- resolve_loocv_local(local, iid = is_loocv_iid(object))
  validate_loocv_level(level)
  loocv_val <- get_loocv.ssn_lm(object, cv_predict = TRUE, se.fit = TRUE, local = local)
  response_val <- model.response(model.frame(object))
  error_val <- response_val - loocv_val$cv_predict
  se_val <- loocv_val$se.fit
  validate_loocv_se(se_val)

  bias <- mean(error_val)
  # Standardized bias uses the prediction standard error itself.
  std.bias <- mean(error_val / se_val)
  MSPE <- loocv_val$mspe
  RMSPE <- sqrt(loocv_val$mspe)
  std.MSPE <- mean(error_val^2 / se_val^2)
  RAV <- sqrt(mean(se_val^2))
  cor2 <- cor(loocv_val$cv_predict, response_val)^2
  zstat <- abs(error_val / se_val)
  cover.80 <- mean(zstat < qnorm(0.90))
  cover.90 <- mean(zstat < qnorm(0.95))
  cover.95 <- mean(zstat < qnorm(0.975))

  loocv_stats <- tibble(
    bias = bias,
    std.bias = std.bias,
    MSPE = MSPE,
    RMSPE = RMSPE,
    std.MSPE = std.MSPE,
    RAV = RAV,
    cor2 = cor2,
    cover.80 = cover.80,
    cover.90 = cover.90,
    cover.95 = cover.95
  )
  if (interval == "prediction") {
    cover_name <- loocv_coverage_name(level)
    if (!cover_name %in% names(loocv_stats)) {
      loocv_stats[[cover_name]] <- mean(zstat < qnorm(1 - (1 - level) / 2))
    }
  }

  if (!cv_predict && !se.fit) {
    return(loocv_stats)
  } else {
    loocv_out <- list()
    loocv_out$stats <- loocv_stats

    if (cv_predict) {
      loocv_out$cv_predict <- loocv_val$cv_predict
    }

    if (se.fit) {
      loocv_out$se.fit <- loocv_val$se.fit
    }
    return(loocv_out)
  }
}

validate_loocv_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1) {
    stop("level must be a single number strictly between 0 and 1.", call. = FALSE)
  }
}

is_loocv_iid <- function(object) {
  is.null(object$random) && inherits(coef(object, "nugget"), "nugget_nugget") &&
    all(vapply(c("tailup", "taildown", "euclid"), function(type) {
      inherits(coef(object, type), paste0(type, "_none"))
    }, logical(1)))
}

resolve_loocv_local <- function(local, iid = FALSE) {
  if (is.logical(local) && (length(local) != 1 || is.na(local))) {
    stop("local must be TRUE, FALSE, or a local-control list.", call. = FALSE)
  }
  if (!is.logical(local) && !is.list(local)) {
    stop("local must be TRUE, FALSE, or a local-control list.", call. = FALSE)
  }
  local_list <- get_local_list_prediction(local)
  if (local_list$method == "covariance" &&
      (!is.numeric(local_list$size) || length(local_list$size) != 1 ||
       !is.finite(local_list$size) || local_list$size < 1 || local_list$size != as.integer(local_list$size))) {
    stop("local$size must be a positive integer.", call. = FALSE)
  }
  local_list
}

validate_loocv_se <- function(se.fit) {
  if (any(!is.finite(se.fit) | se.fit <= 0)) {
    stop("Gaussian LOOCV standardized statistics require finite, positive prediction standard errors.", call. = FALSE)
  }
}

loocv_coverage_name <- function(level) {
  paste0("cover.", sub("^0\\.", "", as.character(level)))
}

get_loocv.ssn_lm <- function(object, cv_predict = FALSE, se.fit = FALSE, local = FALSE, ...) {
  iid <- is_loocv_iid(object)
  local_list <- resolve_loocv_local(local, iid = iid)

  if (iid) {
    model_frame <- model.frame(object)
    X <- model.matrix(object)
    y <- model.response(model_frame)
    model_offset <- model.offset(model_frame)
    y_krige <- if (is.null(model_offset)) y else y - as.vector(model_offset)
    qr_X <- qr(X)
    leverage <- rowSums(qr.Q(qr_X)^2)
    # PRESS updates coefficients after deletion without an n-by-n covariance.
    cv_predict_val <- y - qr.resid(qr_X, y_krige) / (1 - leverage)
    if (se.fit) {
      total_var <- max(coef(object, "nugget")[["nugget"]], object$diagtol, 0)
      cv_predict_se <- sqrt(total_var / (1 - leverage))
    }
  } else if (local_list$method == "all") {
    cov_matrix_val <- covmatrix(object)

    # actually need inverse because of HW blocking
    cov_matrixInv_val <- chol2inv(chol(forceSymmetric(cov_matrix_val)))
    model_frame <- model.frame(object)
    X <- model.matrix(object)
    y <- model.response(model_frame)
    # The kriging update models the offset-free response. Add each held-out
    # row's known offset back after the update, as in predict.ssn_lm().
    model_offset <- model.offset(model_frame)
    y_krige <- if (is.null(model_offset)) y else y - as.vector(model_offset)
    yX <- cbind(y_krige, X)
    SigInv_yX <- cov_matrixInv_val %*% yX

    cv_predict_val_list <- run_pred_dispatch(
      get_loocv, seq_len(object$n), local_list,
      Sig = cov_matrix_val,
      SigInv = cov_matrixInv_val, Xmat = X, y = y_krige, yX = yX,
      SigInv_yX = SigInv_yX, se.fit = se.fit
    )
    # cv_predict_val <- unlist(cv_predict_val_list)
    cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    if (!is.null(model_offset)) {
      cv_predict_val <- cv_predict_val + as.vector(model_offset)
    }
    if (se.fit) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  } else {
    cov_matrix_val <- covmatrix(object)
    model_frame <- model.frame(object)
    X <- model.matrix(object)
    y <- model.response(model_frame)
    model_offset <- model.offset(model_frame)
    y_krige <- if (is.null(model_offset)) y else y - as.vector(model_offset)
    cv_predict_val_list <- run_pred_dispatch(
      get_loocv_local_lm, seq_len(object$n), local_list,
      Sig = cov_matrix_val, Xmat = X, y = y_krige,
      size = local_list$size, se.fit = se.fit,
      betahat = coef(object), cov_betahat = vcov(object)
    )
    cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    if (!is.null(model_offset)) {
      cv_predict_val <- cv_predict_val + as.vector(model_offset)
    }
    if (se.fit) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  }
  if (cv_predict) {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), cv_predict = as.vector(cv_predict_val), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), cv_predict = as.vector(cv_predict_val))
    }
  } else {
    if (se.fit) {
      cv_output <- list(mspe = mean((cv_predict_val - y)^2), se.fit = as.vector(cv_predict_se))
    } else {
      cv_output <- mean((cv_predict_val - y)^2)
    }
  }
  cv_output
}

get_loocv_local_index <- function(obs, Sig, size) {
  retain <- seq_len(NROW(Sig))[-obs]
  n <- length(retain)
  retain[order(abs(as.numeric(Sig[obs, retain])))[seq.int(n, max(1L, n - size + 1L))]]
}

get_loocv_local_lm <- function(obs, Sig, Xmat, y, size, se.fit, betahat, cov_betahat) {
  retain <- get_loocv_local_index(obs, Sig, size)
  new_Sig <- Sig[retain, retain, drop = FALSE]
  new_SigInv <- chol2inv(chol(forceSymmetric(new_Sig)))
  new_X <- Xmat[retain, , drop = FALSE]
  new_y <- y[retain]
  obs_c <- Sig[obs, retain, drop = FALSE]
  new_pred <- Xmat[obs, , drop = FALSE] %*% betahat +
    obs_c %*% new_SigInv %*% (new_y - new_X %*% betahat)
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c %*% new_SigInv %*% new_X
    var_fit <- Sig[obs, obs] - obs_c %*% new_SigInv %*% t(obs_c) +
      Q %*% cov_betahat %*% t(Q)
    se_fit <- sqrt(as.numeric(var_fit))
  } else {
    se_fit <- NULL
  }
  list(pred = as.numeric(new_pred), se.fit = se_fit)
}

#' Get loocv residual
#'
#' @param obs An observation to leave out
#' @param Sig The full covariance matrix
#' @param SigInv The full inverse covariance matrix
#' @param Xmat Model matrix
#' @param y response vector
#'
#' @return A loocv residual
#'
#' @noRd
get_loocv <- function(obs, Sig, SigInv, Xmat, y, yX, SigInv_yX, se.fit) {
  # Rather than refit the model n times with each observation dropped (which
  # would mean inverting a new (n-1)x(n-1) matrix each time), this uses a
  # partitioned-matrix (Sherman-Morrison-type) update: SigInv for the (n-1)
  # remaining observations can be recovered algebraically from the full-data
  # SigInv and the row/column for the dropped observation "obs".
  SigInv_mm <- SigInv[obs, obs] # a constant
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  newX <- Xmat[-obs, , drop = FALSE]
  newyX <- yX[-obs, , drop = FALSE]

  # SigInv %*% yX for the data with "obs" removed, obtained via the partitioned
  # inverse update rather than recomputing SigInv from scratch
  new_SigInv_oo_newyX <- SigInv_yX[-obs, , drop = FALSE] - SigInv_om %*% yX[obs, , drop = FALSE]
  newSigInv_newyX <- new_SigInv_oo_newyX - SigInv_om %*% (crossprod(SigInv_om, newyX) / SigInv_mm)

  newSigInv_newX <- newSigInv_newyX[, -1, drop = FALSE]
  newSigInv_newy <- newSigInv_newyX[, 1, drop = FALSE]
  # refit betahat using only the n-1 remaining observations
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, newSigInv_newX))))
  new_betahat <- new_covbetahat %*% crossprod(newX, newSigInv_newy)
  obs_c <- Sig[obs, -obs, drop = FALSE]
  # kriging predictor at "obs" using the refit betahat and the covariance
  # between "obs" and the remaining observations (universal kriging equation)
  new_pred <- Xmat[obs, , drop = FALSE] %*% new_betahat + obs_c %*% (newSigInv_newy - newSigInv_newX %*% new_betahat)

  # var
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c %*% newSigInv_newX
    new_SigInv_oo_obs_c <- tcrossprod(SigInv[-obs, -obs], obs_c) - SigInv_om %*% (crossprod(SigInv_om, t(obs_c)) / SigInv_mm)
    # kriging prediction variance: marginal variance minus variance explained
    # by the observed data, plus extra uncertainty from estimating betahat
    var_fit <- Sig[obs, obs] - obs_c %*% new_SigInv_oo_obs_c + Q %*% tcrossprod(new_covbetahat, Q)
    se_fit <- sqrt(var_fit)
  } else {
    se_fit <- NULL
  }

  # return
  list(pred = as.numeric(new_pred), se.fit = as.numeric(se_fit))
}
