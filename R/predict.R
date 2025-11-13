#' Model predictions (Kriging)
#'
#' @description Predicted values and intervals based on a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param newdata A character vector that indicates the name of the prediction data set
#'   for which predictions are desired (accessible via \code{object$ssn.object$preds}).
#'   Note that the prediction data must be in the original SSN object used to fit the model.
#'   If \code{newdata} is omitted, predictions
#'   for all prediction data sets are returned. Note that the name \code{".missing"}
#'   indicates the prediction data set that contains the missing observations in the data used
#'   to fit the model.
#' @param se.fit A logical indicating if standard errors are returned.
#'   The default is \code{FALSE}.
#' @param scale A numeric constant by which to scale the regular standard errors and intervals.
#'   Similar to but slightly different than \code{scale} for [stats::predict.lm()], because
#'   predictions form a spatial model may have different residual variances for each
#'   observation in \code{newdata}. The default is \code{NULL}, which returns
#'   the regular standard errors and intervals.
#' @param df Degrees of freedom to use for confidence or prediction intervals
#'   (ignored if \code{scale} is not specified). The default is \code{Inf}.
#' @param interval Type of interval calculation. The default is \code{"none"}.
#'   Other options are \code{"confidence"} (for confidence intervals) and
#'   \code{"prediction"} (for prediction intervals).
#' @param level Tolerance/confidence level. The default is \code{0.95}.
#' @param terms If \code{type} is \code{"terms"}, the type of terms to be returned,
#'   specified via either numeric position or name. The default is all terms are included.
#' @param block A logical indicating whether a block prediction over the entire
#'  region in \code{newdata} should be returned. The default is \code{FALSE}, which returns point
#'  predictions for each location in \code{newdata}. Currently only available for
#'  model fit using \code{ssn_lm()} or models fit using \code{ssn_glm()} where
#'  \code{family} is \code{"gaussian"}.
#' @param local A optional logical or list controlling the big data approximation. If omitted, \code{local}
#'   is set to \code{TRUE} or \code{FALSE} based on the observed data sample size (i.e., sample size of the fitted
#'   model object) -- if the sample size exceeds 10,000, \code{local} is
#'   set to \code{TRUE}, otherwise it is set to \code{FALSE}. This default behavior
#'   occurs because main computational
#'   burden of the big data approximation depends almost exclusively on the
#'   observed data sample size, not the number of predictions desired
#'   (which we feel is not intuitive at first glance).
#'   If \code{local} is \code{FALSE}, no big data approximation
#'   is implemented. If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{method}: The big data approximation method. If \code{method = "all"},
#'       all observations are used and \code{size} is ignored.
#'       If \code{method = "covariance"}, the \code{size} data observations
#'       having the average highest covariance with the prediction locations are used.
#'       The default
#'       is \code{"covariance"}. Only used with models fit using [ssn_lm()].
#'     \item \code{size}: The number of data observations to use when \code{method}
#'       is \code{"distance"} or \code{"covariance"}. The default is 4000. Only used
#'       with models fit using [ssn_lm()].
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. This can significantly speed
#'       up computations even when \code{method = "all"} (i.e., no big data
#'       approximation is used), as predictions
#'       are spread out over multiple cores. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 4000, method = "covariance", parallel = FALSE)}.
#' @param terms If \code{type} is \code{"terms"}, the type of terms to be returned,
#'   specified via either numeric position or name. The default is all terms are included.
#' @param na.action Missing (\code{NA}) values in \code{newdata} will return an error and should
#'   be removed before proceeding.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details The (empirical) best linear unbiased predictions (i.e., Kriging
#'   predictions) at each site are returned when \code{interval} is \code{"none"}
#'   or \code{"prediction"} alongside standard errors. Prediction intervals
#'   are also returned if \code{interval} is \code{"prediction"}. When
#'   \code{interval} is \code{"confidence"}, the estimated mean is returned
#'   alongside standard errors and confidence intervals for the mean.
#'
#' @return If \code{se.fit} is \code{FALSE}, \code{predict.ssn()} returns
#'   a vector of predictions or a matrix of predictions with column names
#'   \code{fit}, \code{lwr}, and \code{upr} if \code{interval} is \code{"confidence"}
#'   or \code{"prediction"}. If \code{se.fit} is \code{TRUE}, a list with the following components is returned:
#'   \itemize{
#'     \item \code{fit}: vector or matrix as above
#'     \item \code{se.fit:} standard error of each fit
#'   }
#'
#' @name predict.SSN2
#' @method predict ssn_lm
#' @export
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, predpts = "pred1km", overwrite = TRUE)
#'
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' predict(ssn_mod, "pred1km")
predict.ssn_lm <- function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf, interval = c("none", "confidence", "prediction"),
                           level = 0.95, type = c("response", "terms"), block = FALSE, local, terms = NULL, na.action = na.fail, ...) {


  # match interval argument so the three display
  interval <- match.arg(interval)
  type <- match.arg(type)

  # call predict_block if necessary
  # if (block) { # gives ::: warning so no rd exported
  #   call_val <- match.call()
  #   call_val[[1]] <- as.symbol("predict_block")
  #   call_list <- as.list(call_val)
  #   call_list <- call_list[-which(names(call_list) %in% c("block"))]
  #   # call_list[[1]] <- quote(SSN2:::predict_block)
  #   call_val <- as.call(call_list)
  #   object <- eval(call_val, envir = parent.frame())
  #   return(object)
  # }

  if (missing(newdata)) newdata <- "all"
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

  # if (is.list(local) || local && (object$n > 15000 || NROW(object$ssn.object$preds[[newdata]]) * object$n > 15000^2)) {
  #   object <- predict_bigdata_ssn_lm(object, newdata, se.fit, interval, level, block, local, ...)
  #   return(object)
  # }

  #  safter but potentially passes block
  if (block) {
    object <- predict_block(object, newdata, se.fit, interval, level, local, ...)
    return(object)
  }

  # deal with local (omitted for now)
  # if (missing(local)) local <- NULL
  # if (is.null(local)) local <- FALSE

  # new data name
  if (missing(newdata)) newdata <- NULL
  newdata_name <- newdata
  # iterate through prediction names
  if (is.null(newdata_name) || newdata_name == "all") {
    newdata_name <- names(object$ssn.object$preds)
  }
  if (length(newdata_name) > 1) {
    pred_list <- lapply(newdata_name, function(x) predict(object, x, se.fit = se.fit, interval = interval, level = level, local = local, ...))
    names(pred_list) <- newdata_name
    return(pred_list)
  }

  if (newdata_name == ".missing") {
    add_newdata_rows <- TRUE
  } else {
    add_newdata_rows <- FALSE
  }

  # rename relevant quantities
  obdata <- object$ssn.object$obs

  # newdata and newdata name
  newdata <- object$ssn.object$preds[[newdata_name]]

  # stop if zero rows
  if (NROW(newdata) == 0) {
    return(NULL)
  }

  # get params object
  params_object <- object$coefficients$params_object

  # make covariance object
  cov_vector <- covmatrix(object, newdata_name)
  if (local_list$method == "covariance") {
    cov_vector_means <- colMeans(cov_vector)
    cov_index <- order(as.numeric(cov_vector_means))[seq(from = object$n, to = max(1, object$n - local$size + 1))]
    cov_vector <- cov_vector[, cov_index, drop = FALSE]
  } else {
    cov_index <- NULL
  }
  cov_vector_list <- split(cov_vector, seq_len(NROW(cov_vector)))

  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  if (any(grepl("nmatrix.", attributes(formula_newdata)$dataClasses, fixed = TRUE)) && NROW(newdata) == 1) {
    newdata <- newdata[c(1, 1), , drop = FALSE]
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
    newdata_model <- newdata_model[1, , drop = FALSE]
    # find offset
    offset <- model.offset(newdata_model_frame)
    if (!is.null(offset)) {
      offset <- offset[1]
    }
    newdata <- newdata[1, , drop = FALSE]
  } else {
    newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
    # assumes that predicted observations are not outside the factor levels
    newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
    # find offset
    offset <- model.offset(newdata_model_frame)
  }

  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # call terms if needed
  if (type == "terms") {
    return(predict_terms(object, newdata_model, se.fit, scale, df, interval, level, add_newdata_rows, terms, ...))
  }

  # storing newdata as a list
  newdata_rows_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_list <- mapply(
    x = newdata_rows_list, y = newdata_model_list, c = cov_vector_list,
    FUN = function(x, y, c) list(row = x, x0 = y, c0 = c), SIMPLIFY = FALSE
  )

  # storing cov matrix
  cov_matrix_val <- covmatrix(object)



  # total var (could do sums params object)
  total_var <- cov_matrix_val[1, 1]

  if (interval %in% c("none", "prediction")) {
    # local_list <- get_local_list_prediction(local)

    if (local_list$method == "all") {
      cov_lowchol <- t(chol(cov_matrix_val))
    } else {
      cov_matrix_val <- cov_matrix_val[cov_index, cov_index, drop = FALSE]
      cov_lowchol <- t(chol(cov_matrix_val))
    }

    Xmat <- model.matrix(object)
    y <- model.response(model.frame(object))
    offset <- model.offset(model.frame(object))

    if (!is.null(cov_index)) {
      Xmat <- Xmat[cov_index, , drop = FALSE]
      y <- y[cov_index]
      if (!is.null(offset)) {
        offset <- offset[cov_index]
      }
    }
    # until big data back
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      pred_val <- parallel::parLapply(cl, newdata_list, get_pred,
                                      se.fit = se.fit, interval = interval, formula = object$formula,
                                      obdata = obdata, cov_matrix_val = cov_matrix_val, total_var = total_var, cov_lowchol = cov_lowchol,
                                      Xmat = Xmat, y = y,
                                      offset = offset,
                                      betahat = coefficients(object), cov_betahat = vcov(object),
                                      contrasts = object$contrasts, local = local_list,
                                      xlevels = object$xlevels, cov_index)
      cl <- parallel::stopCluster(cl)
    } else {
      pred_val <- lapply(newdata_list, get_pred,
                         se.fit = se.fit, interval = interval, formula = object$formula,
                         obdata = obdata, cov_matrix_val = cov_matrix_val, total_var = total_var, cov_lowchol = cov_lowchol,
                         Xmat = Xmat, y = y,
                         offset = offset,
                         betahat = coefficients(object), cov_betahat = vcov(object),
                         contrasts = object$contrasts, local = local_list,
                         xlevels = object$xlevels, cov_index)
    }

    # pred_val <- lapply(newdata_list, get_pred,
    #   se.fit = se.fit, interval = interval, formula = object$formula,
    #   obdata = obdata, cov_matrix_val = cov_matrix_val, total_var = total_var, cov_lowchol = cov_lowchol,
    #   Xmat = model.matrix(object), y = model.response(model.frame(object)),
    #   offset = model.offset(model.frame(object)),
    #   betahat = coefficients(object), cov_betahat = vcov(object),
    #   contrasts = object$contrasts, local = local_list,
    #   xlevels = object$xlevels
    # )


    if (interval == "none") {
      fit <- vapply(pred_val, function(x) x$fit, numeric(1))
      if (se.fit) {
        vars <- vapply(pred_val, function(x) x$var, numeric(1))
        se <- sqrt(vars)
        if (!is.null(scale)) {
          se <- se * scale
        }
        if (add_newdata_rows) {
          names(fit) <- object$missing_index
          names(se) <- object$missing_index
        }
        return(list(fit = fit, se.fit = se))
      } else {
        if (add_newdata_rows) {
          names(fit) <- object$missing_index
        }
        return(fit)
      }
    }

    if (interval == "prediction") {
      fit <- vapply(pred_val, function(x) x$fit, numeric(1))
      vars <- vapply(pred_val, function(x) x$var, numeric(1))
      se <- sqrt(vars)
      if (!is.null(scale)) {
        se <- se * scale
        df <- df
      } else {
        df <- Inf
      }
      tstar <- qt(1 - (1 - level) / 2, df = df)
      # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
      # tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      fit <- cbind(fit, lwr, upr)
      row.names(fit) <- seq_len(NROW(fit))
      if (se.fit) {
        if (add_newdata_rows) {
          row.names(fit) <- object$missing_index
          names(se) <- object$missing_index
        }
        return(list(fit = fit, se.fit = se))
      } else {
        if (add_newdata_rows) {
          row.names(fit) <- object$missing_index
        }
        return(fit)
      }
    }
  } else if (interval == "confidence") {
    # finding fitted values of the mean parameters
    fit <- as.numeric(newdata_model %*% coef(object))
    vars <- as.numeric(vapply(newdata_model_list, function(x) crossprod(x, vcov(object) %*% x), numeric(1)))
    se <- sqrt(vars)
    if (!is.null(scale)) {
      se <- se * scale
      df <- df
    } else {
      df <- Inf
    }
    tstar <- qt(1 - (1 - level) / 2, df = df)
    # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
    # tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    fit <- cbind(fit, lwr, upr)
    row.names(fit) <- seq_len(NROW(fit))
    if (se.fit) {
      if (add_newdata_rows) {
        row.names(fit) <- object$missing_index
        names(se) <- object$missing_index
      }
      return(list(fit = fit, se.fit = se))
    } else {
      if (add_newdata_rows) {
        row.names(fit) <- object$missing_index
      }
      return(fit)
    }
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}

  ' Title
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
  #' @param offset A possible offset
  #' @param betahat Fixed effect estimates
  #' @param cov_betahat Covariance of fixed effects
  #' @param contrasts Possible contrasts
  #' @param local Local neighborhood options (not yet implemented)
  #' @param xlevels Levels of explanatory variables
  #'
  #' @noRd
  get_pred <- function(newdata_list, se.fit, interval, formula, obdata, cov_matrix_val, total_var, cov_lowchol,
                       Xmat, y, offset, betahat, cov_betahat, contrasts, local, xlevels, cov_index) {



    cov_vector_val <- newdata_list$c0

    # moved indexing of relevant quantities outside the function (cov_index no longer needed)
    # if (!is.null(cov_index)) {
    #   obdata <- obdata[cov_index, , drop = FALSE]
    #   model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass, xlev = xlevels)
    #   Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
    #   y <- model.response(model_frame)
    #   offset <- model.offset(model_frame)
    #   # Xmat <- Xmat[cov_index, , drop = FALSE]
    #   # y <- y[cov_index]
    #   # if (!is.null(offset)) {
    #   #   offset <- offset[cov_index]
    #   # }
    # }

    # if (local$method == "covariance") {
    #   n <- length(cov_vector_val)
    #   cov_index <- order(as.numeric(cov_vector_val))[seq(from = n, to = max(1, n - local$size + 1))] # use abs() here?
    #   obdata <- obdata[cov_index, , drop = FALSE]
    #   cov_vector_val <- cov_vector_val[cov_index]
    #   cov_matrix_val <- cov_matrix_val[cov_index, cov_index, drop = FALSE]
    #   cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(cov_matrix_val)))
    #   model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass, xlev = xlevels)
    #   Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
    #   y <- model.response(model_frame)
    #   offset <- model.offset(model_frame)
    # }

    # handle offset
    if (!is.null(offset)) {
      y <- y - offset
    }



    c0 <- as.numeric(cov_vector_val)
    SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
    SqrtSigInv_y <- forwardsolve(cov_lowchol, y)
    residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat
    SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
    x0 <- newdata_list$x0

    fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
    H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)

    if (se.fit || interval == "prediction") {
      total_var <- total_var
      var <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
      pred_list <- list(fit = fit, var = var)
    } else {
      pred_list <- list(fit = fit)
    }
    pred_list
}

predict_block <- function(object, newdata, se.fit, interval, level, local, ...) {
  # deal with local (omitted for now)
  # if (missing(local)) local <- NULL
  # if (is.null(local)) local <- FALSE

  # new data name
  if (missing(newdata)) newdata <- NULL
  newdata_name <- newdata
  # iterate through prediction names
  if (is.null(newdata_name) || newdata_name == "all") {
    newdata_name <- names(object$ssn.object$preds)
  }
  if (length(newdata_name) > 1) {
    pred_list <- lapply(newdata_name, function(x) predict_block(object, x, se.fit = se.fit, interval = interval, level = level, local = local, ...))
    names(pred_list) <- newdata_name
    return(pred_list)
  }

  # rename relevant quantities
  obdata <- object$ssn.object$obs

  # newdata and newdata name
  newdata <- object$ssn.object$preds[[newdata_name]]

  # stop if zero rows
  if (NROW(newdata) == 0) {
    return(NULL)
  }

  # get params object
  params_object <- object$coefficients$params_object

  # make covariance object
  cov_vector <- covmatrix(object, newdata_name)
  # adjustment one for local
  if (local$method == "covariance") {
    cov_vector_means <- colMeans(cov_vector)
    cov_index <- order(as.numeric(cov_vector_means))[seq(from = object$n, to = max(1, object$n - local$size + 1))]
    cov_vector <- cov_vector[, cov_index, drop = FALSE]
  } else {
    cov_index <- NULL
  }
  formula_newdata <- delete.response(terms(object))
  newdata_model_frame <- model.frame(formula_newdata, newdata, drop.unused.levels = FALSE, na.action = na.pass, xlev = object$xlevels)
  # assumes that predicted observations are not outside the factor levels
  newdata_model <- model.matrix(formula_newdata, newdata_model_frame, contrasts = object$contrasts)
  # find offset
  offset <- model.offset(newdata_model_frame)

  # newdata model stuff
  attr_assign <- attr(newdata_model, "assign")
  attr_contrasts <- attr(newdata_model, "contrasts")
  keep_cols <- which(colnames(newdata_model) %in% colnames(model.matrix(object)))
  newdata_model <- newdata_model[, keep_cols, drop = FALSE]
  attr(newdata_model, "assign") <- attr_assign[keep_cols]
  attr(newdata_model, "contrasts") <- attr_contrasts

  # block Kriging stuff

  Xmat <- model.matrix(object)
  x0 <- colMeans(newdata_model)
  c0 <- colMeans(cov_vector)
  y <- model.response(model.frame(object))
  fitted_val <- fitted(object)
  if (!is.null(offset)) {
    y <- y - offset
  }

  # storing cov matrix
  cov_matrix_val <- covmatrix(object)

  # adjustment two for local
  if (local$method == "covariance") {
   Xmat <- Xmat[cov_index, , drop = FALSE]
   y <- y[cov_index]
   cov_matrix_val <- cov_matrix_val[cov_index, cov_index, drop = FALSE]
   fitted_val <- fitted_val[cov_index]
  }



  # cholesky
  cov_lowchol <- t(chol(cov_matrix_val))

  # total var (could do sums params object)
  cov_matrix_preds <- covmatrix(object, newdata = newdata_name, cov_type = "pred.pred")
  total_var <- mean(cov_matrix_preds)

  # betahat
  betahat <- coefficients(object)
  cov_betahat <- vcov(object)

  # block Kriging prediction
  if (interval %in% c("none", "prediction")) {
    # should use cholesky matrix not eigendecomposition matrix
    # ie can't call residuals(object, type = "pearson") directly
    residuals_pearson <- forwardsolve(cov_lowchol, y - fitted_val)
    SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
    fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
    names(fit) <- "1"

    if (interval == "none") {
      if (se.fit) {
        SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
        H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
        vars <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
        se <- sqrt(vars)
        names(se) <- "1"
        return(list(fit = fit, se.fit = se))
      } else {
        return(fit)
      }
    } else if (interval == "prediction") {
      SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
      H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
      vars <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
      se <- sqrt(vars)
      tstar <- qnorm(1 - (1 - level) / 2)
      lwr <- fit - tstar * se
      upr <- fit + tstar * se
      fit <- cbind(fit, lwr, upr)
      row.names(fit) <- "1"
      if (se.fit) {
        names(se) <- "1"
        return(list(fit = fit, se.fit = se))
      } else {
        return(fit)
      }
    }
  } else if (interval == "confidence") {
    fit <- as.numeric(x0 %*% betahat)
    vars <- as.numeric(crossprod(x0, cov_betahat %*% x0))
    se <- sqrt(vars)
    # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
    tstar <- qnorm(1 - (1 - level) / 2)
    lwr <- fit - tstar * se
    upr <- fit + tstar * se
    fit <- cbind(fit, lwr, upr)
    row.names(fit) <- "1"
    if (se.fit) {
      names(se) <- "1"
      return(list(fit = fit, se.fit = se))
    } else {
      return(fit)
    }
  } else {
    stop("Interval must be none, confidence, or prediction")
  }
}
