#' Perform leave-one-out cross validation
#'
#' @description Perform leave-one-out cross validation with options for computationally
#'   efficient approximations for big data.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param cv_predict A logical indicating whether the leave-one-out fitted values
#'   should be returned. Defaults to \code{FALSE}.
#' @param se.fit A logical indicating whether the leave-one-out
#'   prediction standard errors should be returned. Defaults to \code{FALSE}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Each observation is held-out from the data set and the remaining data
#'   are used to make a prediction for the held-out observation. This is compared
#'   to the true value of the observation and several model-fit statistics are computed
#'   across all observations.
#'
#' @return If \code{cv_predict = FALSE} and \code{se.fit = FALSE},
#'   a tibble indicating several
#'   leave-one-out cross validation error metrics. If \code{cv_predict = TRUE} or \code{se.fit = TRUE},
#'   a list with elements: \code{stats}, a tibble indicating several
#'   leave-one-out cross validation metrics; \code{cv_predict}, a numeric vector
#'   with leave-one-out predictions for each observation (if \code{cv_predict = TRUE});
#'   and \code{se.fit}, a numeric vector with leave-one-out prediction standard
#'   errors for each observation (if \code{se.fit = TRUE}).
#'
#'   If an \code{ssn_lm} object, the cross validation error metrics are:
#'   \itemize{
#'     \item bias: The average difference between the predicted value and true value
#'     \item std.bias: The average standardized difference between the predicted value and true value
#'     \item MSPE: The average squared difference between the predicted value and true value
#'     \item RMSPE: The root average squared difference between the predicted value and true value
#'     \item std.MSPE: The average standardized squared difference between the predicted value and true value
#'     \item RAV: The root of the average estimated variance of the predicted value
#'     \item cor2: The squared correlation between the predicted and true values
#'     \item cover.80: Coverage rates of 80% prediction intervals built for the true values
#'     \item cover.90: Coverage rates of 90% prediction intervals built for the true values
#'     \item cover.95: Coverage rates of 95% prediction intervals built for the true values
#'   }
#'
#'   If an \code{ssn_glm} object, the cross validation error metrics are:
#'   \itemize{
#'     \item bias: The average difference between the predicted value and true value
#'     \item MSPE: The average squared difference between the predicted value and true value
#'     \item RMSPE: The root average squared difference between the predicted value and true value
#'     \item RAV: The root of the average estimated variance of the predicted value (on the link scale)
#'   }
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
loocv.ssn_lm <- function(object, cv_predict = FALSE, se.fit = FALSE, ...) {
  loocv_val <- get_loocv.ssn_lm(object, cv_predict = TRUE, se.fit = TRUE)
  response_val <- model.response(model.frame(object))
  error_val <- loocv_val$cv_predict - response_val
  se_val <- loocv_val$se.fit

  bias <- mean(error_val)
  std.bias <- mean(error_val / sqrt(se_val))
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

#' A helper to get the leave one out predictions
#'
#' @param object Model object
#' @param cv_predict Whether cross validation predictions should be returned
#' @param se.fit Whether the standard error should be returned
#' @param ... Additional arguments
#'
#' @noRd
get_loocv.ssn_lm <- function(object, cv_predict = FALSE, se.fit = FALSE, ...) {
  # local stuff (save for later)
  # if (missing(local)) local <- NULL
  # if (is.null(local)) {
  #   if (object$n > 5000) {
  #     local <- TRUE
  #     message("Because the sample size exceeds 5000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun loocv() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
  #   } else {
  #     local <- FALSE
  #   }
  # }
  local <- FALSE
  local_list <- get_local_list_prediction(local)

  if (local_list$method == "all") {
    cov_matrix_val <- covmatrix(object)

    # actually need inverse because of HW blocking
    cov_matrixInv_val <- chol2inv(chol(forceSymmetric(cov_matrix_val)))
    model_frame <- model.frame(object)
    X <- model.matrix(object)
    y <- model.response(model_frame)
    yX <- cbind(y, X)
    SigInv_yX <- cov_matrixInv_val %*% yX

    # parallel stuff
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      cv_predict_val_list <- parallel::parLapply(cl, seq_len(object$n), get_loocv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX, se.fit = se.fit
      )
      cl <- parallel::stopCluster(cl)
    } else {
      cv_predict_val_list <- lapply(seq_len(object$n), get_loocv,
        Sig = cov_matrix_val,
        SigInv = cov_matrixInv_val, Xmat = X, y = y, yX = yX,
        SigInv_yX = SigInv_yX, se.fit = se.fit
      )
    }
    # cv_predict_val <- unlist(cv_predict_val_list)
    cv_predict_val <- vapply(cv_predict_val_list, function(x) x$pred, numeric(1))
    if (se.fit) {
      cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    }
  } else {
    # model_frame <- model.frame(object)
    # X <- model.matrix(object)
    # y <- model.response(model_frame)
    # betahat <- coef(object)
    # cov_betahat <- vcov(object)
    #
    # cov_matrix_val <- covmatrix(object)
    # total_var <- cov_matrix_val[1, 1]
    #
    # if (local_list$parallel) {
    #   # turn of parallel as it is used different in predict
    #   local_list$parallel <- FALSE
    #   cl <- parallel::makeCluster(local_list$ncores)
    #   cv_predict_val_list <- parallel::parLapply(
    #     cl, seq_len(object$n), loocv_local, se.fit, cov_matrix_val,
    #     total_var, X, y, betahat, cov_betahat, local_list
    #   )
    #   cl <- parallel::stopCluster(cl)
    # } else {
    #   cv_predict_val_list <- lapply(
    #     seq_len(object$n), loocv_local, se.fit, cov_matrix_val,
    #     total_var, X, y, betahat, cov_betahat, local_list
    #   )
    # }
    # if (se.fit) {
    #   cv_predict_val <- vapply(cv_predict_val_list, function(x) x$fit, numeric(1))
    #   cv_predict_se <- vapply(cv_predict_val_list, function(x) x$se.fit, numeric(1))
    # } else {
    #   cv_predict_val <- unlist(cv_predict_val_list)
    # }
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

# not currently used (but will be later when local prediction implemented)
loocv_local <- function(row, se.fit, cov_matrix_val, total_var, Xmat, y, betahat, cov_betahat, local) {
  # new_cov_matrix_val <- cov_matrix_val[-row, -row, drop = FALSE]
  # new_cov_vector_val <- cov_matrix_val[row, -row, drop = FALSE]
  # new_Xmat <- Xmat[-row, , drop = FALSE]
  # new_y <- y[-row]
  #
  # if (local$method == "covariance") { # always on covariance as all gets routed to full cholesky and this is not called
  #   n <- length(new_cov_vector_val)
  #   cov_index <- order(as.numeric(new_cov_vector_val))[seq(from = n, to = max(1, n - local$size + 1))] # use abs() here?
  #   new_cov_vector_val <- new_cov_vector_val[cov_index]
  #   new_cov_matrix_val <- new_cov_matrix_val[cov_index, cov_index, drop = FALSE]
  #   cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(new_cov_matrix_val)))
  #   new_Xmat <- new_Xmat[cov_index, , drop = FALSE]
  #   new_y <- new_y[cov_index]
  # }
  #
  # c0 <- as.numeric(new_cov_vector_val)
  # SqrtSigInv_X <- forwardsolve(cov_lowchol, new_Xmat)
  # SqrtSigInv_y <- forwardsolve(cov_lowchol, new_y)
  # residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
  # SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
  # x0 <- Xmat[row, , drop = FALSE]
  #
  # fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
  # H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  # if (se.fit) {
  #   total_var <- total_var
  #   var <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
  #   pred_list <- list(fit = fit, se.fit = sqrt(var))
  # } else {
  #   pred_list <- list(fit = fit)
  # }
  # pred_list
}

#' Get leave one out cross validation value
#'
#' @param obs Index of observed data
#' @param Sig Covariance matrix
#' @param SigInv Inverse covariance matrix
#' @param Xmat Model matrix
#' @param y response
#' @param yX Matrix containing model matrix and response
#' @param SigInv_yX Product of inverse covariance matrix with model matrix and response
#' @param se.fit Whether to return standard errors of leave one out predictions
#'
#' @noRd
get_loocv <- function(obs, Sig, SigInv, Xmat, y, yX, SigInv_yX, se.fit) {
  SigInv_mm <- SigInv[obs, obs] # a constant
  SigInv_om <- SigInv[-obs, obs, drop = FALSE]

  newX <- Xmat[-obs, , drop = FALSE]
  newyX <- yX[-obs, , drop = FALSE]

  new_SigInv_oo_newyX <- SigInv_yX[-obs, , drop = FALSE] - SigInv_om %*% yX[obs, , drop = FALSE]
  newSigInv_newyX <- new_SigInv_oo_newyX - SigInv_om %*% (crossprod(SigInv_om, newyX) / SigInv_mm)

  newSigInv_newX <- newSigInv_newyX[, -1, drop = FALSE]
  newSigInv_newy <- newSigInv_newyX[, 1, drop = FALSE]
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(newX, newSigInv_newX))))
  new_betahat <- new_covbetahat %*% crossprod(newX, newSigInv_newy)
  obs_c <- Sig[obs, -obs, drop = FALSE]
  new_pred <- Xmat[obs, , drop = FALSE] %*% new_betahat + obs_c %*% (newSigInv_newy - newSigInv_newX %*% new_betahat)

  # var
  if (se.fit) {
    Q <- Xmat[obs, , drop = FALSE] - obs_c %*% newSigInv_newX
    new_SigInv_oo_obs_c <- tcrossprod(SigInv[-obs, -obs], obs_c) - SigInv_om %*% (crossprod(SigInv_om, t(obs_c)) / SigInv_mm)
    var_fit <- Sig[obs, obs] - obs_c %*% new_SigInv_oo_obs_c + Q %*% tcrossprod(new_covbetahat, Q)
    se_fit <- sqrt(var_fit)
  } else {
    se_fit <- NULL
  }

  # return
  list(pred = as.numeric(new_pred), se.fit = as.numeric(se_fit))
}
