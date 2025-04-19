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
#' @param interval Type of interval calculation. The default is \code{"none"}.
#'   Other options are \code{"confidence"} (for confidence intervals) and
#'   \code{"prediction"} (for prediction intervals).
#' @param level Tolerance/confidence level. The default is \code{0.95}.
#' @param block A logical indicating whether a block prediction over the entire
#'  region in \code{newdata} should be returned. The default is \code{FALSE}, which returns point
#'  predictions for each location in \code{newdata}. Currently only available for
#'  model fit using \code{ssn_lm()} or models fit using \code{ssn_glm()} where
#'  \code{family} is \code{"gaussian"}.
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
predict_bigdata_ssn_lm <- function(object, newdata, se.fit = FALSE, interval = c("none", "confidence", "prediction"),
                           level = 0.95, block = FALSE, local, ...) {
  # match interval argument so the three display
  interval <- match.arg(interval)


  #  safter but potentially passes block
  # if (block) {
  #   object <- predict_block(object, newdata, se.fit, interval, level, ...)
  #   return(object)
  # }

  # new data name
  if (missing(newdata)) newdata <- NULL
  newdata_name <- newdata
  # iterate through prediction names
  if (is.null(newdata_name) || newdata_name == "all") {
    newdata_name <- names(object$ssn.object$preds)
  }
  if (length(newdata_name) > 1) {
    pred_list <- lapply(newdata_name, function(x) predict_bigdata_ssn_lm(object, x, se.fit = se.fit, interval = interval, level = level, local = local, ...))
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
  obdata_netgeom <- ssn_get_netgeom(obdata)
  network_index_obs <- obdata_netgeom$NetworkID
  pid_obs <- obdata_netgeom$pid

  # newdata and newdata name
  newdata <- object$ssn.object$preds[[newdata_name]]

  # stop if zero rows
  if (NROW(newdata) == 0) {
    return(NULL)
  }

  # get params object
  params_object <- object$coefficients$params_object
  initial_object <-


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

  # storing newdata as a list
  newdata_rows_list <- split(newdata, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))

  # storing newdata as a list
  newdata_list <- mapply(
    x = newdata_rows_list, y = newdata_model_list,
    FUN = function(x, y) list(row = x, x0 = y), SIMPLIFY = FALSE
  )


  if (interval %in% c("none", "prediction")) {

    # local prediction list
    local_list <- get_local_list_prediction(local)

    dotlist <- list(...)
    dotlist_names <- names(dotlist)


    if ("extra_randcov_list" %in% dotlist_names && !is.null(dotlist[["extra_randcov_list"]])) {
      extra_randcov_list <- dotlist$extra_randcov_list
    } else {
      extra_randcov_list <- get_extra_randcov_list(object, obdata, newdata)
    }
    reform_bar2_list <- extra_randcov_list$reform_bar2_list
    Z_index_obdata_list <- extra_randcov_list$Z_index_obdata_list
    reform_bar1_list <- extra_randcov_list$reform_bar1_list
    Z_val_obdata_list <- extra_randcov_list$Z_val_obdata_list

    if ("extra_partition_list" %in% dotlist_names && !is.null(dotlist[["extra_partition_list"]])) {
      extra_partition_list <- dotlist$extra_partition_list
    } else {
      extra_partition_list <- get_extra_partition_list(object, obdata, newdata)
    }
    reform_bar2 <- extra_partition_list$reform_bar2
    partition_index_obdata <- extra_partition_list$partition_index_obdata

    if (local_list$method == "all") {
      cov_lowchol <- t(chol(covmatrix(object)))
    } else {
      cov_lowchol <- NULL
    }

    if (local_list$parallel) {
      cl <- parallel::makeCluster(local_list$ncores)
      pred_val <- parallel::parLapply(cl, newdata_list, get_pred_bigdata_ssn_lm,
                                      se.fit = se.fit, interval = interval, formula = object$formula,
                                      obdata = obdata, params_object = params_object,
                                      random = object$random, reform_bar2_list,
                                      Z_index_obdata_list, reform_bar1_list, Z_val_obdata_list, object$partition_factor,
                                      reform_bar2, partition_index_obdata, cov_lowchol = cov_lowchol,
                                      Xmat = model.matrix(object), y = model.response(model.frame(object)),
                                      offset = model.offset(model.frame(object)),
                                      betahat = coefficients(object), cov_betahat = vcov(object),
                                      contrasts = object$contrasts, local = local_list,
                                      xlevels = object$xlevels, network_index_obs, pid_obs, object$ssn.object, newdata_name, object$additive, object$anisotropy)
      cl <- parallel::stopCluster(cl)
    } else {
      pred_val <- lapply(newdata_list, get_pred_bigdata_ssn_lm,
                         se.fit = se.fit, interval = interval, formula = object$formula,
                         obdata = obdata, params_object = params_object,
                         random = object$random, reform_bar2_list,
                         Z_index_obdata_list, reform_bar1_list, Z_val_obdata_list,  object$partition_factor,
                         reform_bar2, partition_index_obdata, cov_lowchol = cov_lowchol,
                         Xmat = model.matrix(object), y = model.response(model.frame(object)),
                         offset = model.offset(model.frame(object)),
                         betahat = coefficients(object), cov_betahat = vcov(object),
                         contrasts = object$contrasts, local = local_list,
                         xlevels = object$xlevels, network_index_obs, pid_obs, object$ssn.object, newdata_name, object$additive, object$anisotropy)
    }


    if (interval == "none") {
      fit <- vapply(pred_val, function(x) x$fit, numeric(1))
      if (se.fit) {
        vars <- vapply(pred_val, function(x) x$var, numeric(1))
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
      vars <- vapply(pred_val, function(x) x$var, numeric(1))
      se <- sqrt(vars)
      # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
      tstar <- qnorm(1 - (1 - level) / 2)
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
    # tstar <- qt(1 - (1 - level) / 2, df = object$n - object$p)
    tstar <- qnorm(1 - (1 - level) / 2)
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

#' Title
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
get_pred_bigdata_ssn_lm <- function(newdata_list, se.fit, interval, formula,
                                    obdata, params_object,
                                    random, reform_bar2_list,
                                    Z_index_obdata_list, reform_bar1_list, Z_val_obdata_list, partition_factor,
                                    reform_bar2, partition_index_obdata, cov_lowchol,
                                    Xmat, y, offset, betahat, cov_betahat,
                                    contrasts, local, xlevels, network_index_obs, pid_obs, ssn.object, newdata_name, additive, anisotropy) {




  # storing partition vector
  partition_vector <- partition_vector(partition_factor,
                                       data = obdata,
                                       newdata = newdata_list$row, reform_bar2 = reform_bar2,
                                       partition_index_data = partition_index_obdata
  )

  # subsetting partition vector (efficient but causes problems later with
  # random effect subsetting)
  if (!is.null(partition_vector) && local$method %in% c("distance", "covariance") &&
      (is.null(random) || !labels(terms(partition_factor)) %in% labels(terms(random)))) {
    partition_index <- as.vector(partition_vector) == 1
    Z_index_obdata_list <- lapply(Z_index_obdata_list, function(x) {
      x$reform_bar2_vals <- x$reform_bar2_vals[partition_index]
      x
    })
    obdata <- obdata[partition_index, , drop = FALSE]
    partition_vector <- Matrix(1, nrow = 1, ncol = NROW(obdata))
  }

  dist_pred_object <- get_dist_pred_object_bigdata(newdata_list$row, network_index_obs, pid_obs, ssn.object, newdata_name, params_object, additive, anisotropy)

  # making random vector if necessary
  randcov_params_val <- params_object$randcov
  if (!is.null(randcov_params_val)) {
    randcov_vector_val <- randcov_vector(randcov_params_val, obdata, newdata_list$row, reform_bar2_list, Z_index_obdata_list)
  } else {
    randcov_vector_val <- NULL
  }

  # making the covariance vector
  cov_vector_val <- get_cov_vector(params_object, dist_pred_object, obdata, newdata_list$row, partition_factor = partition_factor, anisotropy)

  if (local$method == "covariance") {
    n <- length(cov_vector_val)
    # want the largest covariance here and order goes from smallest first to largest last (keep last values which are largest covariance)
    cov_index <- order(as.numeric(cov_vector_val))[seq(from = n, to = max(1, n - local$size + 1))] # use abs() here?
    obdata <- obdata[cov_index, , drop = FALSE]
    obdata <- obdata[order(cov_index), , drop = FALSE]
    cov_vector_val <- cov_vector_val[cov_index]
  }

  if (local$method %in% c("covariance")) {
    if (!is.null(random)) {
      randcov_names <- get_randcov_names(random)
      xlev_list <- lapply(Z_index_obdata_list, function(x) x$reform_bar2_xlev)
      randcov_Zs <- get_randcov_Zs(obdata, randcov_names, xlev_list = xlev_list)
    } else {
      randcov_Zs <- NULL
    }
    partition_matrix_val <- partition_matrix(partition_factor, obdata)

    new_ssn.object <- list(obs = obdata, path = ssn.object$path) # original order
    new_initial_object <- get_initial_cheap(params_object)
    # new_pid_obs <- pid_obs[cov_index]
    dist_object_oblist <- get_dist_object_bigdata(
      new_ssn.object,
      new_initial_object,
      additive,
      anisotropy,
      local_index = rep(1, length(cov_index)), # just to make it work with big data observed function
      observed_index = rep(TRUE, length(cov_index)) # just to make it work with big data observed function
    )$dist_matlist[[1]]
    dist_object_oblist$network_index <- rep(1, length(cov_index)) # just so nugget cov works
    cov_matrix_val <- get_cov_matrix(params_object, dist_object_oblist, randcov_list = randcov_Zs, partition_list = partition_matrix_val, anisotropy = anisotropy, de_scale = 1, diagtol = 0)
    cov_matrix_val <- cov_matrix_val[order(order(cov_index)), order(order(cov_index))]
    cov_lowchol <- t(Matrix::chol(Matrix::forceSymmetric(cov_matrix_val)))
    model_frame <- model.frame(formula, obdata[order(order(cov_index)), , drop = FALSE], drop.unused.levels = TRUE, na.action = na.pass, xlev = xlevels)
    Xmat <- model.matrix(formula, model_frame, contrasts = contrasts)
    y <- model.response(model_frame)
    offset <- model.offset(model_frame)
  }

  # handle offset
  if (!is.null(offset)) {
    y <- y - offset
  }




  total_var <- cov_lowchol[1, 1]^2

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

# extra lists for use with loocv (do nothing but clean with predict)
get_extra_randcov_list <- function(object, obdata, newdata) {
  # random stuff
  if (!is.null(object$random)) {
    randcov_names <- get_randcov_names(object$random)
    # this causes a memory leak and was not even needed
    # randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
    # comment out here for simple
    reform_bar_list <- lapply(randcov_names, function(randcov_name) {
      bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
      reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
      if (bar_split[[1]] != "1") {
        reform_bar1 <- reformulate(bar_split[[1]], intercept = FALSE)
      } else {
        reform_bar1 <- NULL
      }
      list(reform_bar2 = reform_bar2, reform_bar1 = reform_bar1)
    })
    reform_bar2_list <- lapply(reform_bar_list, function(x) x$reform_bar2)
    names(reform_bar2_list) <- randcov_names
    reform_bar1_list <- lapply(reform_bar_list, function(x) x$reform_bar1)
    names(reform_bar1_list) <- randcov_names
    Z_index_obdata_list <- lapply(reform_bar2_list, function(reform_bar2) {

      reform_bar2_mf <- model.frame(reform_bar2, obdata)
      reform_bar2_terms <- terms(reform_bar2_mf)
      reform_bar2_xlev <- .getXlevels(reform_bar2_terms, reform_bar2_mf)
      reform_bar2_mx <- model.matrix(reform_bar2, obdata)
      reform_bar2_names <- colnames(reform_bar2_mx)
      reform_bar2_split <- split(reform_bar2_mx, seq_len(NROW(reform_bar2_mx)))
      reform_bar2_vals <- reform_bar2_names[vapply(reform_bar2_split, function(y) which(as.logical(y)), numeric(1))]

      # adding dummy levels if newdata observations of random effects are not in original data
      # terms object is unchanged if levels change
      # reform_bar2_mf_new <- model.frame(reform_bar2, newdata)
      # reform_bar2_mf_full <- model.frame(reform_bar2, merge(obdata, newdata, all = TRUE))
      # reform_bar2_terms_full <- terms(rbind(reform_bar2_mf, reform_bar2_mf_new))
      reform_bar2_xlev_full <- .getXlevels(reform_bar2_terms, rbind(reform_bar2_mf, model.frame(reform_bar2, newdata)))
      if (!identical(reform_bar2_xlev, reform_bar2_xlev_full)) {
        reform_bar2_xlev <- reform_bar2_xlev_full
      }

      list(reform_bar2_vals = reform_bar2_vals, reform_bar2_xlev = reform_bar2_xlev)
    })
    # Z_index_obdata_list <- lapply(reform_bar2_list, function(reform_bar2) as.vector(model.matrix(reform_bar2, obdata)))
    names(Z_index_obdata_list) <- randcov_names
    Z_val_obdata_list <- lapply(reform_bar1_list, function(reform_bar1) {
      if (is.null(reform_bar1)) {
        return(NULL)
      } else {
        return(as.vector(model.matrix(reform_bar1, obdata)))
      }
    })
    names(Z_val_obdata_list) <- randcov_names
  } else {
    reform_bar2_list <- NULL
    Z_index_obdata_list <- NULL
    reform_bar1_list <- NULL
    Z_val_obdata_list <- NULL
  }
  list(reform_bar1_list = reform_bar1_list, reform_bar2_list = reform_bar2_list,
       Z_val_obdata_list = Z_val_obdata_list, Z_index_obdata_list = Z_index_obdata_list)
}

get_extra_partition_list <- function(object, obdata, newdata) {
  # partition factor stuff
  if (!is.null(object$partition_factor)) {
    partition_factor_val <- get_partition_name(labels(terms(object$partition_factor)))
    bar_split <- unlist(strsplit(partition_factor_val, " | ", fixed = TRUE))
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
    p_reform_bar2_mf <- model.frame(reform_bar2, obdata)
    p_reform_bar2_terms <- terms(p_reform_bar2_mf)
    p_reform_bar2_xlev <- .getXlevels(p_reform_bar2_terms, p_reform_bar2_mf)
    p_reform_bar2_mx <- model.matrix(reform_bar2, obdata)
    p_reform_bar2_names <- colnames(p_reform_bar2_mx)
    p_reform_bar2_split <- split(p_reform_bar2_mx, seq_len(NROW(p_reform_bar2_mx)))
    p_reform_bar2_vals <- p_reform_bar2_names[vapply(p_reform_bar2_split, function(y) which(as.logical(y)), numeric(1))]


    # adding dummy levels if newdata observations of random effects are not in original data
    # terms object is unchanged if levels change
    # p_reform_bar2_mf_new <- model.frame(reform_bar2, newdata)
    # reform_bar2_mf_full <- model.frame(reform_bar2, merge(obdata, newdata, all = TRUE))
    # p_reform_bar2_terms_full <- terms(rbind(p_reform_bar2_mf, p_reform_bar2_mf_new))
    p_reform_bar2_xlev_full <- .getXlevels(p_reform_bar2_terms, rbind(p_reform_bar2_mf, model.frame(reform_bar2, newdata)))
    if (!identical(p_reform_bar2_xlev, p_reform_bar2_xlev_full)) {
      p_reform_bar2_xlev <- p_reform_bar2_xlev_full
    }

    partition_index_obdata <- list(reform_bar2_vals = p_reform_bar2_vals, reform_bar2_xlev = p_reform_bar2_xlev)
    # partition_index_obdata <- as.vector(model.matrix(reform_bar2, obdata))
  } else {
    reform_bar2 <- NULL
    partition_index_obdata <- NULL
  }
  list(reform_bar2 = reform_bar2, partition_index_obdata = partition_index_obdata)
}
