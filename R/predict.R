#' Model predictions (Kriging)
#'
#' @description Predicted values and intervals based on a fitted model object.
#'
#' @param object A fitted model object from [ssn_lm()], [ssn_glm()],
#'   [ssn_lmRF()], or [ssn_decorrelate()].
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
#'   \code{"prediction"} (for prediction intervals). When \code{interval}
#'   is \code{"none"} or \code{"prediction"}, predictions are returned (and when
#'   requested, their corresponding uncertainties). When \code{interval}
#'   is \code{"confidence"}, mean estimates are returned (and when
#'   requested, their corresponding uncertainties). This \code{"none"} behavior
#'   differs from that of \code{lm()}, as \code{lm()} returns confidence
#'   uncertainties (in \code{.$se.fit}).
#' @param level Tolerance/confidence level. The default is \code{0.95}.
#' @param terms If \code{type} is \code{"terms"}, the type of terms to be returned,
#'   specified via either numeric position or name. The default is all terms are included.
#' @param block A logical indicating whether a block prediction over the entire
#'  region in \code{newdata} should be returned. When \code{block} is \code{TRUE},
#'  \code{newdata} should be a dense grid of prediction locations that span
#'  the entire region. The default is \code{FALSE}, which returns point
#'  predictions for each location in \code{newdata}.
#' @param local A optional logical or list controlling the big data approximation. If omitted, \code{local}
#'   is set to \code{TRUE} or \code{FALSE} based on the observed data sample size (i.e., sample size of the fitted
#'   model object) -- if the sample size exceeds 10,000, \code{local} is
#'   set to \code{TRUE}, otherwise it is set to \code{FALSE}. This default behavior
#'   occurs because for point prediction the main computational
#'   burden of the big data approximation depends almost exclusively on the
#'   observed data sample size, not the number of predictions desired
#'   (which we feel is not intuitive at first glance). For block prediction
#'   (\code{block = TRUE}) the density of the prediction grid also matters, so
#'   \code{local} is additionally set to approximate the theoretical solution
#'   when \code{nrow(newdata)} exceeds 10,000 (see below).
#'   If \code{local} is \code{FALSE}, no big data approximation
#'   is implemented. If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{method}: The big data approximation method. If \code{method = "all"},
#'       all observations are used and \code{size} is ignored.
#'       If \code{method = "covariance"}, the \code{size} data observations
#'       having the largest absolute covariance with each prediction location
#'       are used.
#'       The default is \code{"covariance"}.
#'     \item \code{size}: The number of data observations to use when \code{method}
#'       is \code{"covariance"}. The default is 200.
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. This can significantly speed
#'       up computations even when \code{method = "all"} (i.e., no big data
#'       approximation is used), as predictions
#'       are spread out over multiple cores. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'     \item \code{chunk_size}: Maximum number of prediction rows per
#'       covariance-construction chunk (default 1000). Relevant primarily
#'       with prediction for large \code{newdata}.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 200, method = "covariance", chunk_size = 1000, parallel = FALSE)}.
#'
#'   If \code{block} is \code{TRUE}, \code{local} controls two separate big
#'   data approximations, one for the observed data and one for the prediction
#'   grid (\code{newdata}):
#'   \itemize{
#'     \item \code{method} and \code{size} act on the observed data exactly as
#'       when \code{block} is \code{FALSE} (\code{method} takes \code{"all"} or
#'       \code{"covariance"}; the default \code{method}
#'       is \code{"covariance"} with \code{size} \code{4000}. This default
#'       \code{size} is much larger than when \code{block} is \code{FALSE}
#'       because block prediction averages covariances and explanatory
#'       variables before prediction.
#'     \item \code{method_new} controls the big data prediction grid density.
#'       \code{method_new = "basis"} (the default)
#'       approximates the block variance using a basis of
#'       \code{size_new} well-spread grid nodes. \code{method_new = "subset"}
#'       approximates the block mean and variance by subsampling \code{newdata}
#'       so that its size is only \code{size_new}. The default \code{size_new}
#'       is \code{4000}.
#'     \item \code{ordering} chooses the \code{size_new} nodes from
#'       \code{newdata} and takes the same values as the \code{ordering}
#'       argument of \code{\link{ssn_simulate}()}/\code{\link{conditional}()} under
#'       \code{approximation = "vecchia"} (\code{"pid"}, \code{"maxmin"}, \code{"grts"},
#'       \code{"random"}, \code{"none"}, \code{"middleout"},
#'       \code{"outsidein"}, \code{"coordinate"}). The default is
#'       \code{"pid"}.
#'     \item \code{parallel} and \code{ncores} parallelize the
#'       covariance-chunk calculations. Their defaults are the same as for
#'       point prediction.
#'     \item \code{chunk_size}: Maximum number of prediction rows per
#'       covariance-construction chunk (default 1000). Relevant primarily
#'       with prediction for large \code{newdata}.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(method = "covariance", size = 4000, method_new = "basis",
#'   size_new = 4000, ordering = "pid", chunk_size = 1000)}.
#' @param na.action Missing (\code{NA}) values in \code{newdata} will return an error and should
#'   be removed before proceeding.
#' @param ... Other arguments. Only used for models fit using \code{ssn_lmRF()}
#'   where \code{...} indicates other
#'   arguments to \code{ranger::predict.ranger()}.
#'
#' @details For \code{ssn_lm} and \code{ssn_glm} objects, the (empirical)
#'   best linear unbiased predictions (i.e., Kriging
#'   predictions) at each site are returned when \code{interval} is \code{"none"}
#'   or \code{"prediction"} alongside standard errors. Prediction intervals
#'   are also returned if \code{interval} is \code{"prediction"}. When
#'   \code{interval} is \code{"confidence"}, the estimated mean is returned
#'   alongside standard errors and confidence intervals for the mean.
#' 
#'   For `ssn_lmRF` objects, random forest spatial residual model
#'   predictions combine the random forest prediction with the (empirical)
#'   best linear unbiased prediction for the residual. This approach is called
#'   random forest regression Kriging.
#' 
#'   For \code{decorrelate} objects, the spatial decorrelation transformation
#'   predictions recorrelated to the original scale. For \code{decorrelate_list}
#'   objects, predictions are returned for each list element.
#'
#' @return For \code{ssn_lm} and \code{ssn_glm} objects, if \code{se.fit}
#'   is \code{FALSE}, \code{predict()} returns
#'   a vector of predictions or a matrix of predictions with column names
#'   \code{fit}, \code{lwr}, and \code{upr} if \code{interval} is \code{"confidence"}
#'   or \code{"prediction"}. If \code{se.fit} is \code{TRUE}, a list with the following components is returned:
#'   \itemize{
#'     \item \code{fit}: vector or matrix as above
#'     \item \code{se.fit:} standard error of each fit
#'   }
#'
#'   For \code{ssn_lmRF} and \code{ssn_decorrelate} objects, a vector of
#'   predictions.
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
                           level = 0.95, type = c("response", "terms", "weight"), block = FALSE, local, terms = NULL, na.action = na.fail, ...) {

  # match interval argument so the three display
  interval <- match.arg(interval)
  type <- match.arg(type)

  if (!is.null(scale) && !is.numeric(scale)) {
    stop("scale must be numeric.", call. = FALSE)
  }

  if (type == "weight") {
    if (block) {
      stop("type = \"weight\" is not supported for block prediction (block = TRUE).", call. = FALSE)
    }
    se.fit <- FALSE
    interval <- "none"
  }

  if (missing(newdata)) newdata <- "all"
  if (missing(local)) local <- NULL
  if (block) {
    return(predict_block(
      object, newdata, se.fit = se.fit, scale = scale, df = df,
      interval = interval, level = level, type = type, local = local,
      terms = terms, na.action = na.action, ...
    ))
  }

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

  # if (missing(local)) local <- NULL
  # if (is.null(local)) local <- FALSE

  # new data name
  if (missing(newdata)) newdata <- NULL
  newdata_name <- newdata
  # iterate through prediction names
  newdata_name <- resolve_newdata_name(object, newdata_name)
  if (length(newdata_name) > 1) {
    pred_list <- lapply(newdata_name, function(x) {
      predict(object, x,
        se.fit = se.fit, scale = scale, df = df, interval = interval, level = level,
        type = type, terms = terms, na.action = na.action, local = local, ...
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

  nm <- get_newdata_model_matrix(object, newdata)
  newdata <- nm$newdata
  newdata_model <- nm$newdata_model
  newdata_offset <- nm$offset

  # call terms if needed
  if (type == "terms") {
    return(predict_terms(object, newdata_model, se.fit, scale, df, interval, level, add_newdata_rows, terms, ...))
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
    return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows, object$missing_index))
  }

  # make covariance object, in bounded row-chunks over newdata (mirroring
  # block prediction's chunking) so the full n_obs x n_pred covariance is
  # never materialized at once; ctx also bundles the marginal-variance and
  # random-effect setup below so it is built once per call, not once per row
  ctx <- get_point_pred_context(object, newdata_name, newdata, newdata_model, local_list)
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

  pred_val <- run_pred_dispatch(get_pred, newdata_list, local_list = local_list,
    se.fit = se.fit, interval = interval, formula = object$formula,
    obdata = obdata, cov_matrix_val = cov_matrix_val,
    spatial_nugget_var = spatial_nugget_var, randcov_params = randcov_params, cov_lowchol = cov_lowchol,
    randcov_context = randcov_context,
    Xmat = Xmat, y = y, offset = offset,
    betahat = coefficients(object), cov_betahat = vcov(object),
    contrasts = object$contrasts, local = local_list,
    xlevels = object$xlevels, type = type
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
      se <- sqrt(vars)
      if (!is.null(scale)) {
        se <- se * scale
      }
    }
    return(finalize_interval_none(fit, se, add_newdata_rows, object$missing_index))
  }

  if (interval == "prediction") {
    fit <- vapply(pred_val, function(x) x$fit, numeric(1))
    if (!is.null(newdata_offset)) {
      fit <- fit + newdata_offset
    }
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
    return(finalize_interval_bounds(fit, lwr, upr, se, se.fit, add_newdata_rows, object$missing_index))
  }
}

get_pred <- function(newdata_list, se.fit, interval, formula, obdata, cov_matrix_val, spatial_nugget_var, randcov_params, cov_lowchol,
                     Xmat, y, offset, betahat, cov_betahat, contrasts, local, xlevels, type = "response", randcov_context = NULL) {

  cov_vector_val <- newdata_list$c0
  n <- length(cov_vector_val)
  if (local$method == "covariance") {
    keep <- order(abs(as.numeric(cov_vector_val)))[seq.int(n, max(1L, n - local$size + 1L))]
    cov_vector_val <- cov_vector_val[keep]
    cov_lowchol <- t(chol(cov_matrix_val[keep, keep, drop = FALSE]))
    Xmat <- Xmat[keep, , drop = FALSE]
    y <- y[keep]
    if (!is.null(offset)) offset <- offset[keep]
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

  # handle offset
  if (!is.null(offset)) {
    y <- y - offset
  }

  SqrtSigInv_y <- forwardsolve(cov_lowchol, y)
  residuals_pearson <- SqrtSigInv_y - SqrtSigInv_X %*% betahat

  fit <- as.numeric(x0 %*% betahat + Matrix::crossprod(SqrtSigInv_c0, residuals_pearson))
  H <- x0 - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_X)

  if (se.fit || interval == "prediction") {
    total_var <- spatial_nugget_var + randcov_newvar(randcov_params, newdata_list$row, context = randcov_context)
    var <- as.numeric(total_var - Matrix::crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% Matrix::tcrossprod(cov_betahat, H))
    pred_list <- list(fit = fit, var = var)
  } else {
    pred_list <- list(fit = fit)
  }
  pred_list
}

predict_block <- function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
                          interval = c("none", "confidence", "prediction"), level = 0.95,
                          type = c("response", "terms", "weight"), local = NULL,
                          terms = NULL, na.action = na.fail, ...) {
  interval <- match.arg(interval)
  type <- match.arg(type)
  if (identical(type, "weight")) {
    stop("type = \"weight\" is not supported for block prediction (block = TRUE).", call. = FALSE)
  }
  newdata_name <- resolve_newdata_name(object, newdata)
  if (length(newdata_name) > 1L) {
    result <- lapply(newdata_name, function(name) {
      predict_block(
        object, name, se.fit = se.fit, scale = scale, df = df, interval = interval,
        level = level, type = type, local = local, terms = terms,
        na.action = na.action, ...
      )
    })
    names(result) <- newdata_name
    return(result)
  }

  grid <- object$ssn.object$preds[[newdata_name]]
  if (!NROW(grid)) return(NULL)
  local_unset <- is.null(local)
  if (local_unset) {
    if (object$n > 10000 || NROW(grid) > 10000) {
      local <- list(
        method = if (object$n > 10000) "covariance" else "all",
        size = 4000L, method_new = "basis",
        size_new = if (NROW(grid) > 10000) 4000L else Inf,
        ordering = "pid", chunk_size = 1000L, parallel = FALSE
      )
      message("Because the fitted sample size or block grid exceeds 10,000, we are using a computationally efficient block-prediction approximation. To compute the exact solution instead, rerun predict() with local = FALSE.")
    } else {
      local <- FALSE
    }
  }
  local <- get_local_list_prediction_block(local)

  nm <- get_newdata_model_matrix(object, grid)
  grid_size <- NROW(grid)
  nodes <- if (local$size_new >= grid_size) {
    seq_len(grid_size)
  } else {
    get_block_nodes(grid, local$size_new, local$ordering)
  }
  x0 <- matrix(colMeans(nm$newdata_model), nrow = 1)
  colnames(x0) <- colnames(nm$newdata_model)
  offset_block <- if (is.null(nm$offset)) NULL else mean(nm$offset)

  if (identical(type, "terms")) {
    return(predict_terms(
      object, x0, se.fit, scale, df, interval, level, FALSE, terms, ...
    ))
  }
  if (identical(interval, "confidence")) {
    fit <- as.numeric(x0 %*% coefficients(object))
    if (!is.null(offset_block)) fit <- fit + offset_block
    variance <- as.numeric(x0 %*% vcov(object) %*% t(x0))
    se <- sqrt(pmax(variance, 0))
    if (!is.null(scale)) {
      se <- se * scale
    } else {
      df <- Inf
    }
    tstar <- stats::qt(1 - (1 - level) / 2, df = df)
    output <- cbind(fit = fit, lwr = fit - tstar * se, upr = fit + tstar * se)
    rownames(output) <- "1"
    if (se.fit) return(list(fit = output, se.fit = stats::setNames(se, "1")))
    return(output)
  }

  subset_newdata <- identical(local$method_new, "subset") && length(nodes) < grid_size
  covariance_rows <- if (subset_newdata) nodes else seq_len(grid_size)
  need_s0 <- se.fit || identical(interval, "prediction")
  quantities <- get_block_quantities(
    object, newdata_name, c0_rows = covariance_rows, s0_rows = covariance_rows,
    nodes = nodes, chunk_size = local$chunk_size, compute_s0 = need_s0,
    parallel = local$parallel, ncores = local$ncores
  )
  if (subset_newdata && need_s0) {
    # built once and reused for both diagonal sums below, instead of
    # re-deriving the training-side random-effect grouping per row
    randcov_context <- get_randcov_context(
      object$coefficients$params_object$randcov, object$ssn.object$obs, grid
    )
    quantities$s0 <- quantities$s0 -
      get_block_diagonal_sum(object, newdata_name, nodes, randcov_context = randcov_context) / length(nodes)^2 +
      get_block_diagonal_sum(object, newdata_name, seq_len(grid_size), randcov_context = randcov_context) / grid_size^2
  }
  c0 <- quantities$c0
  covariance <- covmatrix(object)
  Xmat <- model.matrix(object)
  y <- model.response(model.frame(object))
  offset <- model.offset(model.frame(object))
  if (!is.null(offset)) y <- y - as.vector(offset)
  if (identical(local$method, "covariance")) {
    cov_index <- order(as.numeric(c0))[seq(from = object$n, to = max(1, object$n - local$size + 1L))]
    c0 <- c0[cov_index]
    covariance <- covariance[cov_index, cov_index, drop = FALSE]
    Xmat <- Xmat[cov_index, , drop = FALSE]
    y <- y[cov_index]
  }

  cov_lowchol <- t(chol(covariance))
  sqrt_siginv_x <- forwardsolve(cov_lowchol, Xmat)
  sqrt_siginv_y <- forwardsolve(cov_lowchol, y)
  residuals_pearson <- sqrt_siginv_y - sqrt_siginv_x %*% coefficients(object)
  sqrt_siginv_c0 <- forwardsolve(cov_lowchol, c0)
  fit <- as.numeric(x0 %*% coefficients(object) + Matrix::crossprod(sqrt_siginv_c0, residuals_pearson))
  if (!is.null(offset_block)) fit <- fit + offset_block
  names(fit) <- "1"
  if (!se.fit && identical(interval, "none")) return(fit)

  H <- x0 - Matrix::crossprod(sqrt_siginv_c0, sqrt_siginv_x)
  variance <- as.numeric(
    quantities$s0 - Matrix::crossprod(sqrt_siginv_c0, sqrt_siginv_c0) +
      H %*% Matrix::tcrossprod(vcov(object), H)
  )
  se <- sqrt(pmax(variance, 0))
  if (!is.null(scale)) {
    se <- se * scale
  } else {
    df <- Inf
  }
  if (identical(interval, "prediction")) {
    tstar <- stats::qt(1 - (1 - level) / 2, df = df)
    output <- cbind(fit = fit, lwr = fit - tstar * se, upr = fit + tstar * se)
    rownames(output) <- "1"
  } else {
    output <- fit
  }
  if (se.fit) return(list(fit = output, se.fit = stats::setNames(se, "1")))
  output
}
