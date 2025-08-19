#' @param type The scale (\code{response} or \code{link}) of predictions obtained
#'   using \code{ssn_glm} objects.
#' @param newdata_size The \code{size} value for each observation in \code{newdata}
#'   used when predicting for the binomial family.
#' @param var_correct A logical indicating whether to return the corrected prediction
#'   variances when predicting via models fit using \code{ssn_glm}. The default is
#'   \code{TRUE}.
#' @param dispersion The dispersion of assumed when computing the prediction standard errors
#'   for \code{ssn_glm()} model objects when \code{family}
#'   is \code{"nbinomial"}, \code{"beta"}, \code{"Gamma"}, or \code{"inverse.gaussian"}.
#'   If omitted, the model object dispersion parameter is used.
#' @rdname predict.SSN2
#' @method predict ssn_glm
#' @export
predict.ssn_glm <- function(object, newdata, type = c("link", "response", "terms"), se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                            level = 0.95, dispersion = NULL, terms = NULL, local, var_correct = TRUE, newdata_size, na.action = na.fail, ...) {
  # match type argument so the two display
  type <- match.arg(type)

  # match interval argument so the three display
  interval <- match.arg(interval)

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

  # handle local for now
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
  dispersion_params_val <- as.vector(coef(object, type = "dispersion"))

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
    # glm supports standard errors for terms objects but not intervals (no interval argument)
    # scale df not used for glms
    return(predict_terms(object, newdata_model, se.fit, scale = NULL, df = Inf, interval, level, add_newdata_rows, terms, ...))
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
    local_list <- get_local_list_prediction(local)

    if (local_list$method == "all") {
      cov_lowchol <- t(chol(cov_matrix_val))
      predvar_adjust_ind <- FALSE
      predvar_adjust_all <- TRUE
    } else {
      cov_matrix_val <- cov_matrix_val[cov_index, cov_index, drop = FALSE]
      cov_lowchol <- t(chol(cov_matrix_val))
      # predvar_adjust_ind <- TRUE
      # predvar_adjust_all <- FALSE
      predvar_adjust_ind <- FALSE
      predvar_adjust_all <- TRUE # changed to this with averaging subsetting
    }

    # # matrix cholesky
    # if (local_list$method == "all") {
    #   cov_matrix_val <- covmatrix(object)
    #   cov_lowchol <- t(chol(cov_matrix_val))
    #   predvar_adjust_ind <- FALSE
    #   predvar_adjust_all <- TRUE
    # } else {
    #   cov_lowchol <- NULL
    #   predvar_adjust_ind <- TRUE
    #   predvar_adjust_all <- FALSE
    # }

    # change predvar adjust based on var correct
    if (!var_correct) {
      predvar_adjust_ind <- FALSE
      predvar_adjust_all <- FALSE
    }

    Xmat <- model.matrix(object)
    y <- model.response(model.frame(object))
    offset <- model.offset(model.frame(object))
    w <- fitted(object, type = "link")
    size <- object$size
    if (!is.null(cov_index)) {
      Xmat <- Xmat[cov_index, , drop = FALSE]
      y <- y[cov_index]
      w <- w[cov_index]
      if (!is.null(offset)) {
        offset <- offset[cov_index]
      }
      if (!is.null(size)) {
        size <- size[cov_index]
      }
    }

    # until big data back
    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      pred_val <- parallel::parLapply(cl, newdata_list, get_pred_glm,
                                      se.fit = se.fit, interval = interval, formula = object$formula,
                       obdata = obdata, cov_matrix_val = cov_matrix_val, total_var = total_var, cov_lowchol = cov_lowchol,
                       Xmat = Xmat, y = y,
                       betahat = coefficients(object), cov_betahat = vcov(object, var_correct = FALSE),
                       contrasts = object$contrasts, local = local_list,
                       family = object$family, w = w,
                       size = size, dispersion = dispersion_params_val,
                       predvar_adjust_ind = predvar_adjust_ind, xlevels = object$xlevels, cov_index = cov_index)
      cl <- parallel::stopCluster(cl)
    } else {
      pred_val <- lapply(newdata_list, get_pred_glm,
                       se.fit = se.fit, interval = interval, formula = object$formula,
                       obdata = obdata, cov_matrix_val = cov_matrix_val, total_var = total_var, cov_lowchol = cov_lowchol,
                       Xmat = Xmat, y = y,
                       betahat = coefficients(object), cov_betahat = vcov(object, var_correct = FALSE),
                       contrasts = object$contrasts, local = local_list,
                       family = object$family, w = w,
                       size = size, dispersion = dispersion_params_val,
                       predvar_adjust_ind = predvar_adjust_ind, xlevels = object$xlevels, cov_index = cov_index)

    }

    # pred_val <- lapply(newdata_list, get_pred_glm,
    #   se.fit = se.fit, interval = interval, formula = object$formula,
    #   obdata = obdata, cov_matrix_val = cov_matrix_val, total_var = total_var, cov_lowchol = cov_lowchol,
    #   Xmat = model.matrix(object), y = model.response(model.frame(object)),
    #   betahat = coefficients(object), cov_betahat = vcov(object, var_correct = FALSE),
    #   contrasts = object$contrasts, local = local_list,
    #   family = object$family, w = fitted(object, type = "link"),
    #   size = object$size, dispersion = dispersion_params_val,
    #   predvar_adjust_ind = predvar_adjust_ind, xlevels = object$xlevels, cov_index
    # )






    if (interval == "none") {
      fit <- vapply(pred_val, function(x) x$fit, numeric(1))
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
      }
      if (type == "response") {
        fit <- invlink(fit, object$family, newdata_size)
      }
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
            c0 = cov_vector,
            cov_index = cov_index,
            cov_betahat = vcov(object, var_correct = FALSE)
          )
          vars <- vars_adj + vars
        }
        se <- sqrt(vars)
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
      # apply offset
      if (!is.null(offset)) {
        fit <- fit + offset
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
          c0 = cov_vector,
          cov_index = cov_index,
          cov_betahat = vcov(object, var_correct = FALSE)
        )
        vars <- vars_adj + vars
      }
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
    # apply offset
    if (!is.null(offset)) {
      fit <- fit + offset
    }
    newdata_model_list <- split(newdata_model, seq_len(NROW(newdata_model)))
    vars <- as.numeric(vapply(newdata_model_list, function(x) crossprod(x, vcov(object) %*% x), numeric(1)))
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
#'
#' @noRd
get_pred_glm <- function(newdata_list, se.fit, interval,
                         formula, obdata, cov_matrix_val, total_var, cov_lowchol,
                         Xmat, y, betahat, cov_betahat, contrasts, local,
                         family, w, size, dispersion, predvar_adjust_ind, xlevels, cov_index) {
  cov_vector_val <- newdata_list$c0

  # moved indexing of relevant quantities outside the function (cov_index no longer needed)
  # if (!is.null(cov_index)) {
  #   obdata <- obdata[cov_index, , drop = FALSE]
  #   model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass, xlev = xlevels)
  #   Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
  #   # Xmat <- Xmat[cov_index, , drop = FALSE]
  #   y <- model.response(model_frame)
  #   # y <- y[cov_index]
  #   offset <- model.offset(model_frame)
  #   # if (!is.null(offset)) {
  #   #   offset <- offset[cov_index]
  #   # }
  #   w <- w[cov_index]
  #   if (!is.null(size)) {
  #     size <- size[cov_index]
  #   }
  # }

  # if (local$method == "covariance") {
  #   # n <- length(cov_vector_val)
  #   # cov_index <- order(as.numeric(cov_vector_val))[seq(from = n, to = max(1, n - local$size + 1))] # use abs() here?
  #   # obdata <- obdata[cov_index, , drop = FALSE]
  #   # cov_vector_val <- cov_vector_val[cov_index]
  #   # cov_matrix_val <- cov_matrix_val[cov_index, cov_index, drop = FALSE]
  #   # w <- w[cov_index]
  #   # y <- y[cov_index]
  #   # if (!is.null(size)) {
  #   #   size <- size[cov_index]
  #   # }
  #   # cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(cov_matrix_val)))
  #   # model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass, xlev = xlevels)
  #   # Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
  # }

  c0 <- as.numeric(cov_vector_val)
  SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
  SqrtSigInv_w <- forwardsolve(cov_lowchol, w)
  residuals_pearson <- SqrtSigInv_w - SqrtSigInv_X %*% betahat
  SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
  x0 <- newdata_list$x0

  fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
  H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)
  if (se.fit || interval == "prediction") {
    var <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
    if (predvar_adjust_ind) {
      var_adj <- get_wts_varw(family, Xmat, y, w, size, dispersion, cov_lowchol, x0, c0, cov_index = cov_index, cov_betahat = cov_betahat)
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
get_wts_varw <- function(family, Xmat, y, w, size, dispersion, cov_lowchol, x0, c0, cov_index, cov_betahat) {


  SigInv <- chol2inv(t(cov_lowchol)) # works on upchol

  if (!is.null(cov_index)) {
    Xmat <- Xmat[cov_index, , drop = FALSE]
    y <- y[cov_index]
    w <- w[cov_index]
    if (!is.null(size)) {
      size <- size[cov_index]
    }
  }
  SigInv_X <- SigInv %*% Xmat
  # cov_betahat <- chol2inv(chol(Matrix::forceSymmetric(crossprod(Xmat, SigInv_X)))) # invertibility issues big data
  cov_betahat <- cov_betahat
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  Ptheta <- SigInv - SigInv_X %*% wts_beta

  d <- get_d(family, w, y, size, dispersion)
  # and then the gradient vector
  # g <-  d - Ptheta %*% w
  # Next, compute H
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
