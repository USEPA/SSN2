#' @param type The prediction type, either on the response scale, link scale (only for
#'   \code{ssn_glm()} model objects), terms scale,
#'   or prediction (i.e., Kriging) weight scale.
#' @param newdata_size The \code{size} value for each observation in \code{newdata}
#'   used when predicting for the binomial family, with a default value of 1.
#' @param var_correct A logical indicating whether to return the corrected prediction
#'   variances when predicting via models fit using \code{ssn_glm()}. The default is
#'   \code{TRUE}.
#' @param dispersion The dispersion assumed when computing the prediction standard errors
#'   for \code{ssn_glm()} model objects when \code{family}
#'   is \code{"nbinomial"}, \code{"beta"}, \code{"Gamma"}, or \code{"inverse.gaussian"}.
#'   If omitted, the model object dispersion parameter is used.
#' @param delta A logical indicating whether to return delta method standard errors
#' on the response scale when \code{se.fit = TRUE} and \code{type = "response"}. The default is \code{FALSE}.
#' @rdname predict.SSN2
#' @method predict ssn_glm
#' @export
predict.ssn_glm <- function(object, newdata, type = c("link", "response", "terms", "weight"), se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                            level = 0.95, block = FALSE, dispersion = NULL, terms = NULL, local, var_correct = TRUE, delta = FALSE, newdata_size, na.action = na.fail, ...) {
  # match type argument so the two display
  type <- match.arg(type)

  # match interval argument so the three display
  interval <- match.arg(interval)

  if (!is.logical(block) || length(block) != 1L || is.na(block)) {
    stop("block must be TRUE or FALSE.", call. = FALSE)
  }
  if (block) {
    stop("Block prediction is not supported for ssn_glm() models.", call. = FALSE)
  }

  if (type == "weight") {
    se.fit <- FALSE
    interval <- "none"
  }

  if (!is.logical(delta)) {
    stop("delta must be TRUE or FALSE", call. = FALSE)
  }

  # new data name
  if (missing(newdata)) newdata <- "all"
  newdata_name <- newdata

  # deal with newdata_size
  if (missing(newdata_size)) newdata_size <- NULL

  # handle dispersion argument if provided
  if (!is.null(dispersion)) {
    if (object$family %in% c("binomial", "poisson") && dispersion != 1) {
      stop("dispersion is fixed at one for binomial and poisson families.", call. = FALSE)
    }
    object$coefficients$params_object$dispersion[1] <- dispersion
  }

  if (missing(local)) local <- NULL
  # deal with local
  if (is.null(local)) {
    if (newdata != "all") {
      if (object$n > 10000 || (object$n * NROW(object$ssn.object$preds[[newdata]]) > 1e8)) {
        # if (object$n > 5000 || NROW(newdata) > 5000) {
        local <- TRUE
        message("Because the large sample size or number of predictions, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun predict() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
      } else {
        local <- FALSE
      }
    }
  }
  # make local list
  local_list <- get_local_list_prediction(local)
  local <- local_list

  # iterate through prediction names
  newdata_name <- resolve_newdata_name(object, newdata_name)
  if (length(newdata_name) > 1) {
    if (!is.null(newdata_size)) {
      stop("newdata_size cannot be used when predicting for multiple newdata sets (newdata = \"all\" or omitted); call predict() separately for the single relevant dataset.", call. = FALSE)
    }
    pred_list <- lapply(newdata_name, function(x) {
      predict(object, x,
        type = type, se.fit = se.fit, interval = interval, level = level,
        terms = terms, var_correct = var_correct, delta = delta, na.action = na.action, local = local, ...
      )
    })
    names(pred_list) <- newdata_name
    return(pred_list)
  }

  pn <- get_prediction_newdata(object, newdata_name)
  obdata <- pn$obdata
  newdata <- pn$newdata
  add_newdata_rows <- pn$add_newdata_rows

  # stop if zero rows
  if (NROW(newdata) == 0) {
    return(NULL)
  }

  dispersion_params_val <- as.vector(coef(object, type = "dispersion"))

  nm <- get_newdata_model_matrix(object, newdata)
  newdata <- nm$newdata
  newdata_model <- nm$newdata_model
  newdata_offset <- nm$offset

  # call terms if needed
  if (type == "terms") {
    # glm supports standard errors for terms objects but not intervals (no interval argument)
    # scale df not used for glms
    return(predict_terms(object, newdata_model, se.fit, scale = NULL, df = Inf, interval, level, add_newdata_rows, terms, ...))
  }

  # confidence intervals for the mean only need the fixed-effect design and
  # coefficient covariance, so return before reading any prediction-to-observed
  # distance/covariance data
  if (interval == "confidence") {
    # finding fitted values of the mean parameters
    fit <- as.numeric(newdata_model %*% coef(object))
    if (!is.null(newdata_offset)) {
      fit <- fit + newdata_offset
    }
    vars <- get_diag_XVXt(newdata_model, vcov(object))
    se <- sqrt(vars)
    # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
    tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    if (type == "response") {
      fit <- invlink(fit, object$family, newdata_size)
      lwr <- invlink(lwr, object$family, newdata_size)
      upr <- invlink(upr, object$family, newdata_size)
    }
    return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows, object$missing_index))
  }

  # make covariance object, in bounded row-chunks over newdata (mirroring
  # block prediction's chunking) so the full n_obs x n_pred covariance is
  # never materialized at once; ctx also bundles the marginal-variance and
  # random-effect setup below so it is built once per call, not once per row
  ctx <- get_point_pred_context(object, newdata_name, newdata, newdata_model, local_list)
  cov_vector_list <- ctx$cov_vector_list
  newdata_list <- ctx$newdata_list
  cov_matrix_val <- ctx$cov_matrix_val
  spatial_nugget_var <- ctx$spatial_nugget_var
  randcov_params <- ctx$randcov_params
  cov_lowchol <- ctx$cov_lowchol
  randcov_context <- ctx$randcov_context
  Xmat <- ctx$Xmat
  y <- ctx$y
  offset <- ctx$offset

  # only "none" and "prediction" reach this point (match.arg() at the top
  # restricts interval to one of "none"/"confidence"/"prediction", and
  # "confidence" already returned above)

  if (local_list$method == "all") {
      predvar_adjust_ind <- FALSE
      predvar_adjust_all <- TRUE
    } else {
      predvar_adjust_ind <- TRUE
      predvar_adjust_all <- FALSE
    }

    # change predvar adjust based on var correct
    if (!var_correct) {
      predvar_adjust_ind <- FALSE
      predvar_adjust_all <- FALSE
    }

    w <- fitted(object, type = "link")
    size <- object$size

    pred_val <- run_pred_dispatch(get_pred_glm, newdata_list, local_list = local_list,
      se.fit = se.fit, interval = interval, formula = object$formula,
      obdata = obdata, cov_matrix_val = cov_matrix_val,
      spatial_nugget_var = spatial_nugget_var, randcov_params = randcov_params, cov_lowchol = cov_lowchol,
      randcov_context = randcov_context,
      Xmat = Xmat, y = y, offset = offset,
      betahat = coefficients(object), cov_betahat = vcov(object, var_correct = FALSE),
      contrasts = object$contrasts, local = local_list,
      family = object$family, w = w,
      size = size, dispersion = dispersion_params_val,
      predvar_adjust_ind = predvar_adjust_ind, xlevels = object$xlevels, type = type
    )

    if (type == "weight") {
      fit <- do.call("rbind", lapply(pred_val, function(x) x$fit))
      fit <- as.matrix(fit)
      colnames(fit) <- object$observed_index
      if (add_newdata_rows) {
        rownames(fit) <- object$missing_index
      }
      return(fit)
    }

    if (interval == "none") {
      fit <- vapply(pred_val, function(x) x$fit, numeric(1))
      if (!is.null(newdata_offset)) {
        fit <- fit + newdata_offset
      }
      se <- NULL
      if (se.fit) {
        vars <- vapply(pred_val, function(x) x$var, numeric(1))
        if (predvar_adjust_all) {
          # predvar_adjust is for the local function so FALSE there is TRUE
          # here
          vars_adj <- get_wts_varw(
            family = object$family,
            Xmat = model.matrix(object),
            y = model.response(model.frame(object)),
            w = fitted(object, type = "link"),
            size = object$size,
            dispersion = dispersion_params_val,
            cov_lowchol = cov_lowchol,
            x0 = newdata_model,
            c0 = do.call(rbind, cov_vector_list)
          )
          vars <- vars_adj + vars
        }
        se <- sqrt(vars)
        if (type == "response" && delta) {
          se <- get_delta_se(fit, se, object$family, newdata_size)
        }
      }
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
      }
      return(finalize_interval_none(fit, se, add_newdata_rows, object$missing_index))
    }

    if (interval == "prediction") {
      fit <- vapply(pred_val, function(x) x$fit, numeric(1))
      if (!is.null(newdata_offset)) {
        fit <- fit + newdata_offset
      }
      vars <- vapply(pred_val, function(x) x$var, numeric(1))
      if (predvar_adjust_all) {
        vars_adj <- get_wts_varw(
          family = object$family,
          Xmat = model.matrix(object),
          y = model.response(model.frame(object)),
          w = fitted(object, type = "link"),
          size = object$size,
          dispersion = dispersion_params_val,
          cov_lowchol = cov_lowchol,
          x0 = newdata_model,
          c0 = do.call(rbind, cov_vector_list)
        )
        vars <- vars_adj + vars
      }
      se <- sqrt(vars)
      # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
      tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      if (type == "response" && se.fit && delta) {
        se <- get_delta_se(fit, se, object$family, newdata_size)
      }
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
        lwr <- invlink(lwr, object$family, newdata_size)
        upr <- invlink(upr, object$family, newdata_size)
      }
      return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows, object$missing_index))
    }
}






#' Get a prediction (and its standard error) for glms
#'
#' @param newdata_list A row of prediction data
#' @param se.fit Whether standard errors should be returned
#' @param interval The interval type
#' @param formula Model formula
#' @param obdata Observed data
#' @param cov_matrix_val Covariance matrix
#' @param total_var Total variance in the process
#' @param cov_lowchol Lower triangular of Cholesky decomposition matrix
#' @param Xmat Model matrix
#' @param y Response variable
#' @param betahat Fixed effect estimates
#' @param cov_betahat Covariance of fixed effects
#' @param contrasts Possible contrasts
#' @param local Local neighborhood options (not yet implemented)
#' @param family glm family
#' @param w Latent effects
#' @param size Number of binomial trials
#' @param dispersion Dispersion parameter
#' @param predvar_adjust_ind Whether prediction variance should be adjusted for uncertainty in w
#' @param xlevels Levels of explanatory variables
#' @param type type scale
#' @param randcov_context random effect context
#'
#' @noRd
get_pred_glm <- function(newdata_list, se.fit, interval,
                         formula, obdata, cov_matrix_val, spatial_nugget_var, randcov_params, cov_lowchol,
                         Xmat, y, offset, betahat, cov_betahat, contrasts, local,
                         family, w, size, dispersion, predvar_adjust_ind, xlevels, type = "link", randcov_context = NULL) {
  cov_vector_val <- newdata_list$c0
  n <- length(cov_vector_val)
  if (local$method == "covariance") {
    keep <- order(abs(as.numeric(cov_vector_val)))[seq.int(n, max(1L, n - local$size + 1L))]
    cov_vector_val <- cov_vector_val[keep]
    cov_lowchol <- t(chol(cov_matrix_val[keep, keep, drop = FALSE]))
    Xmat <- Xmat[keep, , drop = FALSE]
    y <- y[keep]
    w <- w[keep]
    if (!is.null(offset)) offset <- offset[keep]
    if (!is.null(size)) size <- size[keep]
  }

  c0 <- as.numeric(cov_vector_val)
  SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
  x0 <- newdata_list$x0

  if (type == "weight") {
    Xt_SigInv <- t(backsolve(t(cov_lowchol), SqrtSigInv_X))
    betahat_wt <- cov_betahat %*% Xt_SigInv
    residuals_weight <- -1 * Xmat %*% betahat_wt
    diag(residuals_weight) <- diag(residuals_weight) + 1
    fit <- x0 %*% betahat_wt + Matrix::crossprod(SqrtSigInv_c0, forwardsolve(cov_lowchol, residuals_weight))
    if (local$method == "covariance") {
      weights <- matrix(0, 1L, n)
      weights[, keep] <- fit
      fit <- weights
    }
    return(list(fit = fit))
  }

  w_free <- if (!is.null(offset)) w - offset else w
  SqrtSigInv_w <- forwardsolve(cov_lowchol, w_free)
  residuals_pearson <- SqrtSigInv_w - SqrtSigInv_X %*% betahat

  fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
  H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  if (se.fit || interval == "prediction") {
    total_var <- spatial_nugget_var + randcov_newvar(randcov_params, newdata_list$row, context = randcov_context)
    var <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
    if (predvar_adjust_ind) {
      var_adj <- get_wts_varw(family, Xmat, y, w, size, dispersion, cov_lowchol, x0, c0)
      var <- var_adj + var
    }
    pred_list <- list(fit = fit, var = var)
  } else {
    pred_list <- list(fit = fit)
  }
  pred_list
}

#' Get weights by which to adjust variances involving w
#'
#' @param family glm family
#' @param Xmat Model matrix
#' @param y Response variable
#' @param w Latent effects
#' @param size Number of binomial trials
#' @param dispersion Dispersion parameter
#' @param cov_lowchol Lower triangular of Cholesky decomposition matrix
#' @param x0 Explanatory variable values for newdata
#' @param c0 Covariance between observed and newdata
#'
#' @noRd
get_wts_varw <- function(family, Xmat, y, w, size, dispersion, cov_lowchol, x0, c0) {


  SigInv <- chol2inv(t(cov_lowchol)) # works on upchol

  SigInv_X <- SigInv %*% Xmat
  cov_betahat <- chol2inv(chol(crossprod(Xmat, SigInv_X)))
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  Ptheta <- SigInv - SigInv_X %*% wts_beta

  D <- get_D(family, w, y, size, dispersion)
  H <- D - Ptheta
  mHInv <- solve(-H) # chol2inv(chol(Matrix::forceSymmetric(-H))) # solve(-H)

  if (is.vector(x0)) { # for length-one predicts result x0 c0 are vectors (how splm pred operates)
    wts_pred <- x0 %*% wts_beta + c0 %*% SigInv - (c0 %*% SigInv_X) %*% wts_beta
    var_adj <- as.numeric(wts_pred %*% tcrossprod(mHInv, wts_pred))
  } else { # this is to handle the matrix arguments for non-local predict calls with spglm
    if (NROW(x0) == 1) {
      wts_pred <- x0 %*% wts_beta + c0 %*% SigInv - (c0 %*% SigInv_X) %*% wts_beta
      var_adj <- as.numeric(wts_pred %*% tcrossprod(mHInv, wts_pred))
    } else {
      var_adj <- vapply(seq_len(NROW(x0)), function(x) { # this is so that only the diagonal of these products is returned
        x0_new <- x0[x, , drop = FALSE]
        c0_new <- c0[x, , drop = FALSE]
        wts_pred <- x0_new %*% wts_beta + c0_new %*% SigInv - (c0_new %*% SigInv_X) %*% wts_beta
        as.numeric(wts_pred %*% tcrossprod(mHInv, wts_pred))
      }, numeric(1))
    }
  }
  var_adj
}
