#' Perform k-fold cross validation
#'
#' @description Perform k-fold cross validation with options for computationally
#'   efficient approximations for big data. Generalizes [loocv()] (leave-one-out
#'   cross validation) to leaving out \code{k} folds of (approximately) equal size
#'   instead of single observations.
#'
#' @param object A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param k The number of folds. Must be a whole number at least 2 and no more
#'   than the sample size (the number of non-missing observations in \code{object}).
#'   The default is \code{5}. If \code{k} equals the sample size, \code{kcv()}
#'   is equivalent to (and calls) [loocv()]. Ignored when \code{folds_index}
#'   is supplied.
#' @param folds_index An optional vector, the same length as the number of
#'   non-missing observations in \code{object} and in the same order,
#'   assigning each observation to a fold (at least two distinct values
#'   required). When supplied, this fold assignment is used as-is, implying \code{k}
#'   is ignored and no random fold assignment is performed. The default is
#'   \code{NULL}, which randomly assigns \code{k} (approximately) equally
#'   sized folds.
#' @param cv_predict A logical indicating whether the k-fold cross validation fitted values
#'   should be returned. Defaults to \code{FALSE}. If \code{object} is from [ssn_glm()],
#'   the fitted values returned are on the link scale.
#' @param se.fit A logical indicating whether the k-fold cross validation
#'   prediction standard errors should be returned. Defaults to \code{FALSE}.
#'   If \code{object} is from [ssn_glm()],
#'   the standard errors correspond to the fitted values returned on the link scale.
#' @param local A list or logical. If a list, specific list elements described
#'   in [predict.SSN2()] control the big data approximation behavior.
#'   If a logical, \code{TRUE} chooses default list elements for the list version
#'   of \code{local} as specified in [predict.SSN2()]. Defaults to \code{FALSE},
#'   which performs exact computations.
#' @param interval For [ssn_lm()] objects, whether to append empirical coverage
#'   at `level`. The default, `"none"`, retains the legacy 80, 90, and 95
#'   percent coverage columns.
#' @param level Prediction interval level used when `interval = "prediction"`.
#'   The default is `0.95`.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details Observations are randomly partitioned into \code{k} folds of
#'   (approximately) equal size. Each fold is held out from the data set in turn
#'   and the remaining data are used to make predictions for the held-out fold. This is
#'   compared to the true values of the held-out observations and several fit
#'   statistics are (sometimes optionally) computed: bias, mean-squared-prediction error (MSPE),
#'   root-mean-squared-prediction error (RMSPE), and the squared correlation
#'   (cor2) between the observed data and k-fold cross validation predictions
#'   (regarded as a prediction version of r-squared appropriate for comparing
#'   across spatial and nonspatial models), , and 
#'   prediction interval coverage (cover.XX). Generally,
#'   bias should be near zero and prediction interval coverage at the
#'   intended level for well-fitting models. The lower the MSPE and RMSPE, the better the model
#'   fit (according to the k-fold cross validation criterion). The higher the
#'   cor2, the better the model fit (according to the k-fold cross validation
#'   criterion). cor2 and cover.XX are not returned when \code{object} was fit using
#'   \code{ssn_glm()} because we do not observe the underlying latent mean.
#'
#'   When \code{object} is from \code{ssn_lm()}, setting
#'   \code{interval = "prediction"} additionally reports the empirical coverage
#'   of the \code{level} (e.g. 95\%) k-fold cross validation prediction interval --
#'   the proportion of held-out observations whose true value falls within
#'   \code{fit +/- qnorm(1 - (1 - level) / 2) * se.fit} (the same normal-quantile
#'   interval [predict.spmodel()] uses by default). This is only available for
#'   \code{ssn_lm()} objects, since \code{ssn_glm()}
#'   have no observed-scale latent mean to compare against.
#'
#' @return If \code{cv_predict = FALSE} and \code{se.fit = FALSE},
#'   a fit statistics tibble (with bias, MSPE, RMSPE, and cor2; see Details).
#'   If \code{cv_predict = TRUE} or \code{se.fit = TRUE},
#'   a list with elements: \code{stats}, a fit statistics tibble
#'   (with bias, MSPE, RMSPE, and cor2; see Details); \code{cv_predict}, a numeric vector
#'   with k-fold cross validation predictions for each observation (if \code{cv_predict = TRUE});
#'   and \code{se.fit}, a numeric vector with k-fold cross validation prediction standard
#'   errors for each observation (if \code{se.fit = TRUE}). When \code{object} is from
#'   \code{ssn_lm()} and \code{interval = "prediction"}, the fit
#'   statistics tibble also has \code{cover.XX} columns (e.g. \code{cover.95}
#'   for \code{level = 0.95}; see Details).
#'
#' @name kcv.SSN2
#' @method kcv ssn_lm
#' @export
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf <- ssn_import(temp_path, overwrite = TRUE)
#' fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf,
#'   tailup_type = "exponential", additive = "afvArea"
#' )
#' set.seed(1)
#' kcv(fit, k = 3)
kcv.ssn_lm <- function(object, k = 5, cv_predict = FALSE, se.fit = FALSE,
                       local, interval = c("none", "prediction"),
                       level = 0.95, folds_index, ...) {
  if (missing(folds_index)) folds_index <- NULL
  if (missing(local)) local <- NULL
  interval <- match.arg(interval)
  local <- resolve_cv_auto_local(local, object$n, "kcv")
  local <- resolve_kcv_local(local)
  validate_loocv_level(level)

  X <- model.matrix(object)
  folds <- resolve_kcv_folds(object, k, folds_index, X)
  if (is_loocv_folds(folds$fold_list)) {
    return(loocv(
      object, cv_predict = cv_predict, se.fit = se.fit, local = local,
      interval = interval, level = level, ...
    ))
  }

  # Standardized and coverage statistics always require fold standard errors.
  se_needed <- TRUE
  if (local$method == "all") {
    cv <- get_kcv.ssn_lm(object, folds$fold_list, se.fit = se_needed, local = local)
  } else {
    cv <- get_kcv_local.ssn_lm(object, folds$fold_list, se.fit = se_needed, local = local)
  }
  response <- model.response(model.frame(object))
  stats <- get_kcv_lm_stats(cv$cv_predict, response, cv$se.fit, interval, level)

  if (!cv_predict && !se.fit) {
    return(stats)
  }
  output <- list(stats = stats)
  if (cv_predict) output$cv_predict <- cv$cv_predict
  if (se.fit) output$se.fit <- cv$se.fit
  output
}

#' Auto-escalate an unspecified \code{local} argument to \code{TRUE} for large samples
#'
#' Matches spmodel's \code{kcv()}/\code{loocv()} auto-escalation: when the
#' caller omitted \code{local} (already coerced to \code{NULL} via
#' \code{missing()}), a sample size over 5,000 switches to the big-data
#' approximation with a message; otherwise \code{local} resolves to
#' \code{FALSE}. An explicitly supplied \code{local} (including
#' \code{local = FALSE}) is returned unchanged.
#'
#' @param local The caller's \code{local} argument, \code{NULL} if omitted.
#' @param n The fitted model's sample size.
#' @param fun_name The name of the calling exported function (\code{"kcv"} or
#'   \code{"loocv"}), used in the message.
#'
#' @return \code{local}, resolved to \code{TRUE}/\code{FALSE} if it was
#'   \code{NULL}, otherwise unchanged.
#'
#' @noRd
resolve_cv_auto_local <- function(local, n, fun_name) {
  if (is.null(local)) {
    if (n > 5000) {
      local <- TRUE
      message(
        "Because the sample size exceeds 5000, we are setting local = TRUE to perform computationally ",
        "efficient approximations. To override this behavior and compute the exact solution, rerun ",
        fun_name, "() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly ",
        "long computational times."
      )
    } else {
      local <- FALSE
    }
  }
  local
}

#' Validate and normalize the \code{local} argument for \code{kcv()}
#'
#' @param local A logical or list; see the \code{local} argument to
#'   \code{\link{kcv.SSN2}()}.
#'
#' @return The resolved \code{local} list from
#'   \code{\link{get_local_list_prediction}()}.
#'
#' @noRd
resolve_kcv_local <- function(local) {
  if (is.logical(local) && (length(local) != 1 || is.na(local))) {
    stop("local must be TRUE, FALSE, or a local-control list.", call. = FALSE)
  }
  if (!is.logical(local) && !is.list(local)) {
    stop("local must be TRUE, FALSE, or a local-control list.", call. = FALSE)
  }
  local_list <- get_local_list_prediction(local)
  if (local_list$method == "covariance" &&
      (!is.numeric(local_list$size) || length(local_list$size) != 1 ||
       !is.finite(local_list$size) || local_list$size < 1 ||
       local_list$size != as.integer(local_list$size))) {
    stop("local$size must be a positive integer.", call. = FALSE)
  }
  local_list
}

#' Resolve, validate, and split fold assignments for \code{kcv()}
#'
#' Builds random \code{k} folds when \code{folds_index} is \code{NULL},
#' otherwise validates and accepts a supplied fold assignment (in either
#' fitted-row or original-row order), then splits it into per-fold row-index
#' lists and checks the resulting training design.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param k The number of folds (ignored when \code{folds_index} is
#'   supplied).
#' @param folds_index A supplied fold assignment vector, or \code{NULL}.
#' @param X The fitted model's design matrix, used to validate each fold's
#'   training design via \code{\link{validate_kcv_training_design}()}.
#'
#' @return A list with \code{fold_id} (one value per fitted observation) and
#'   \code{fold_list} (a named list of row-index vectors, one per fold).
#'
#' @noRd
resolve_kcv_folds <- function(object, k, folds_index, X) {
  n <- object$n
  original_n <- max(c(object$observed_index, object$missing_index, 0L))
  if (is.null(folds_index)) {
    check_kcv_k(k, n)
    fold_id <- get_kcv_folds(k, n)
  } else {
    if (!is.atomic(folds_index) || !is.null(dim(folds_index))) {
      stop("folds_index must be a vector.", call. = FALSE)
    }
    if (length(folds_index) == n) {
      fold_id <- folds_index
    } else if (length(folds_index) == original_n) {
      fold_id <- folds_index[object$observed_index]
    } else {
      stop("folds_index must have one value per fitted observation or per original obs row.", call. = FALSE)
    }
    check_kcv_folds_index(fold_id, n)
  }
  fold_list <- split(seq_len(n), fold_id, drop = TRUE)
  validate_kcv_training_design(X, fold_list)
  list(fold_id = fold_id, fold_list = fold_list)
}

#' Validate \code{kcv()}'s \code{k} argument
#'
#' @param k The number of folds to validate.
#' @param n The number of fitted observations.
#'
#' @return \code{NULL}, invisibly, if \code{k} is a whole number from 2 to
#'   \code{n}; otherwise an error.
#'
#' @noRd
check_kcv_k <- function(k, n) {
  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) ||
      k != as.integer(k) || k < 2) {
    stop("k must be a single whole number greater than or equal to 2.", call. = FALSE)
  }
  if (k > n) {
    stop("k cannot exceed the number of fitted observations.", call. = FALSE)
  }
  invisible(NULL)
}

#' Randomly assign observations to \code{k} (approximately) equal folds
#'
#' @param k The number of folds.
#' @param n The number of observations to assign.
#'
#' @return An integer vector of length \code{n}, each entry in
#'   \code{seq_len(k)}.
#'
#' @noRd
get_kcv_folds <- function(k, n) {
  sample(rep(seq_len(k), length.out = n))
}

#' Validate a resolved, per-fitted-observation fold assignment
#'
#' @param folds_index A fold assignment vector (already resolved to one
#'   value per fitted observation).
#' @param n The number of fitted observations.
#'
#' @return \code{NULL}, invisibly, if \code{folds_index} has length \code{n},
#'   no missing values, and at least two distinct fold values; otherwise an
#'   error.
#'
#' @noRd
check_kcv_folds_index <- function(folds_index, n) {
  if (length(folds_index) != n) {
    stop("folds_index must have one assignment for every fitted observation.", call. = FALSE)
  }
  if (anyNA(folds_index)) {
    stop("folds_index must assign every fitted observation to a nonmissing fold.", call. = FALSE)
  }
  if (length(unique(folds_index)) < 2) {
    stop("folds_index must assign observations to at least two nonempty folds.", call. = FALSE)
  }
  invisible(NULL)
}

#' Check that every fold leaves a valid, full-rank fixed-effect training design
#'
#' @param X The fitted model's design matrix.
#' @param fold_list A named list of held-out row-index vectors, one per
#'   fold.
#'
#' @return \code{NULL}, invisibly, if every fold leaves more training rows
#'   than fixed effects and a full-rank training design; otherwise an error
#'   naming the offending fold.
#'
#' @noRd
validate_kcv_training_design <- function(X, fold_list) {
  p <- NCOL(X)
  n <- NROW(X)
  for (i in seq_along(fold_list)) {
    train <- setdiff(seq_len(n), fold_list[[i]])
    label <- names(fold_list)[i]
    if (length(train) <= p) {
      stop(
        paste0("Holding out fold '", label, "' leaves ", length(train),
               " training observations for ", p, " fixed-effect coefficients."),
        call. = FALSE
      )
    }
    if (qr(X[train, , drop = FALSE])$rank < p) {
      stop(
        paste0("Holding out fold '", label,
               "' leaves a rank-deficient fixed-effect training design."),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

#' Check whether a fold assignment is equivalent to leave-one-out
#'
#' @param fold_list A named list of per-fold row-index vectors.
#'
#' @return \code{TRUE} if every fold has exactly one observation (so
#'   \code{kcv()} should dispatch to \code{\link{loocv}()}), \code{FALSE}
#'   otherwise.
#'
#' @noRd
is_loocv_folds <- function(fold_list) {
  length(fold_list) == sum(lengths(fold_list)) && all(lengths(fold_list) == 1)
}

#' Compute Gaussian k-fold cross-validation error/coverage statistics
#'
#' @param cv_predict A vector of k-fold cross-validation predictions.
#' @param response The observed response vector.
#' @param se.fit A vector of k-fold cross-validation prediction standard
#'   errors.
#' @param interval Whether to also append coverage at \code{level}; one of
#'   \code{"none"} or \code{"prediction"}.
#' @param level The prediction interval level used when
#'   \code{interval = "prediction"}.
#'
#' @return A tibble with \code{bias}, \code{std.bias}, \code{MSPE},
#'   \code{RMSPE}, \code{std.MSPE}, \code{RAV}, \code{cor2},
#'   \code{cover.80}, \code{cover.90}, \code{cover.95}, and (if requested) the
#'   \code{level}-specific coverage column.
#'
#' @noRd
get_kcv_lm_stats <- function(cv_predict, response, se.fit, interval, level) {
  validate_loocv_se(se.fit)
  error <- response - cv_predict
  zstat <- abs(error / se.fit)
  stats <- tibble(
    bias = mean(error),
    std.bias = mean(error / se.fit),
    MSPE = mean(error^2),
    RMSPE = sqrt(mean(error^2)),
    std.MSPE = mean(error^2 / se.fit^2),
    RAV = sqrt(mean(se.fit^2)),
    cor2 = cor(cv_predict, response)^2,
    cover.80 = mean(zstat < qnorm(0.90)),
    cover.90 = mean(zstat < qnorm(0.95)),
    cover.95 = mean(zstat < qnorm(0.975))
  )
  if (interval == "prediction") {
    cover_name <- loocv_coverage_name(level)
    if (!cover_name %in% names(stats)) {
      stats[[cover_name]] <- mean(zstat < qnorm(1 - (1 - level) / 2))
    }
  }
  stats
}

#' Compute exact (block-update) k-fold cross-validation predictions for a Gaussian model
#'
#' Computes the full covariance matrix and its inverse once, then applies
#' \code{\link{get_kcv}()}'s partitioned-inverse block update per fold
#' (optionally in parallel).
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param fold_list A named list of per-fold row-index vectors.
#' @param se.fit Whether to also return prediction standard errors.
#' @param local A resolved \code{local} list; \code{local$parallel}/
#'   \code{local$ncores} control parallel execution across folds.
#'
#' @return A list with \code{cv_predict} and \code{se.fit} (\code{NULL} if
#'   not requested).
#'
#' @noRd
get_kcv.ssn_lm <- function(object, fold_list, se.fit, local) {
  Sigma <- covmatrix(object)
  Sigma_inv <- chol2inv(chol(forceSymmetric(Sigma)))
  X <- model.matrix(object)
  frame <- model.frame(object)
  offset <- model.offset(frame)
  y <- model.response(frame)
  y_krige <- if (is.null(offset)) y else y - as.vector(offset)
  yX <- cbind(y_krige, X)
  Sigma_inv_yX <- Sigma_inv %*% yX

  values <- run_pred_dispatch(
    get_kcv, fold_list, local,
    Sig = Sigma, SigInv = Sigma_inv, Xmat = X, y = y_krige,
    yX = yX, SigInv_yX = Sigma_inv_yX, se.fit = se.fit
  )

  pred <- numeric(object$n)
  if (se.fit) se <- numeric(object$n)
  for (i in seq_along(fold_list)) {
    held <- fold_list[[i]]
    pred[held] <- values[[i]]$pred
    if (se.fit) se[held] <- values[[i]]$se.fit
  }
  if (!is.null(offset)) pred <- pred + as.vector(offset)
  list(cv_predict = pred, se.fit = if (se.fit) se else NULL)
}

#' Get kcv fold residual
#'
#' Block generalization of \code{\link{get_loocv}()}: rather than holding out a
#' single observation \code{obs}, holds out a whole fold (a vector of row
#' indices) at once via the same partitioned-matrix (Sherman-Morrison-type)
#' update, replacing division by the scalar \code{SigInv[obs, obs]} with
#' \code{solve()} against the \code{m x m} block \code{SigInv[fold, fold]}
#' (\code{m} = fold size). At \code{m = 1} this reduces to exactly
#' \code{get_loocv()}'s formula.
#'
#' @param fold A vector of row indices to leave out together
#' @param Sig The full covariance matrix
#' @param SigInv The full inverse covariance matrix
#' @param Xmat Model matrix
#' @param y response vector
#'
#' @return A kcv residual
#'
#' @noRd
get_kcv <- function(fold, Sig, SigInv, Xmat, y, yX, SigInv_yX, se.fit) {
  train <- -fold
  SigInv_mm <- SigInv[fold, fold, drop = FALSE] # an m x m block (a scalar when m = 1)
  SigInv_om <- SigInv[train, fold, drop = FALSE]
  new_X <- Xmat[train, , drop = FALSE]
  new_yX <- yX[train, , drop = FALSE]

  # SigInv %*% yX for the data with "fold" removed, obtained via the partitioned
  # inverse update rather than recomputing SigInv from scratch
  new_SigInv_oo_yX <- SigInv_yX[train, , drop = FALSE] -
    SigInv_om %*% yX[fold, , drop = FALSE]
  new_SigInv_yX <- new_SigInv_oo_yX -
    SigInv_om %*% solve(SigInv_mm, crossprod(SigInv_om, new_yX))
  new_SigInv_X <- new_SigInv_yX[, -1, drop = FALSE]
  new_SigInv_y <- new_SigInv_yX[, 1, drop = FALSE]
  # refit betahat using only the remaining observations
  new_covbetahat <- chol2inv(chol(forceSymmetric(crossprod(new_X, new_SigInv_X))))
  new_betahat <- new_covbetahat %*% crossprod(new_X, new_SigInv_y)
  held_c <- Sig[fold, train, drop = FALSE]
  # kriging predictor at "fold" using the refit betahat and the covariance
  # between the fold and the remaining observations (universal kriging equation)
  pred <- Xmat[fold, , drop = FALSE] %*% new_betahat +
    held_c %*% (new_SigInv_y - new_SigInv_X %*% new_betahat)

  if (se.fit) {
    Q <- Xmat[fold, , drop = FALSE] - held_c %*% new_SigInv_X
    new_SigInv_train_c <- tcrossprod(SigInv[train, train, drop = FALSE], held_c) -
      SigInv_om %*% solve(SigInv_mm, crossprod(SigInv_om, t(held_c)))
    # kriging prediction covariance among the fold's own points: marginal
    # covariance minus variance explained by the observed data, plus extra
    # uncertainty from estimating betahat -- only the diagonal (marginal SE per
    # point) is returned, matching get_loocv()'s per-observation output shape
    variance <- Sig[fold, fold, drop = FALSE] - held_c %*% new_SigInv_train_c +
      Q %*% tcrossprod(new_covbetahat, Q)
    se <- sqrt(diag(variance))
  } else {
    se <- NULL
  }
  list(pred = as.numeric(pred), se.fit = se)
}

#' Build one fold's held-out-as-missing SSN object for a local refit
#'
#' Restores the fitted model's original missing-response rows (if any) and
#' marks the given fold's response values missing, so a refit on the result
#' treats the fold as its prediction set.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param fold A vector of row indices (into the fitted rows) to hold out.
#'
#' @return A list with \code{ssn} (the modified SSN object) and
#'   \code{held_original} (the fold's row positions in the original,
#'   unsubset data).
#'
#' @noRd
get_kcv_fold_ssn <- function(object, fold) {
  ssn <- object$ssn.object
  if (length(object$missing_index) > 0) {
    missing <- ssn$preds[[".missing"]]
    if (is.null(missing)) {
      stop("The fitted object does not retain its original missing-response rows.", call. = FALSE)
    }
    original_index <- c(object$observed_index, object$missing_index)
    ssn$obs <- rbind(ssn$obs, missing)[order(original_index), , drop = FALSE]
  }
  held_original <- object$observed_index[fold]
  response_vars <- all.vars(object$formula[[2]])
  for (variable in response_vars) {
    if (!variable %in% names(ssn$obs)) {
      stop("Could not identify every response column while creating a k-fold training set.", call. = FALSE)
    }
    ssn$obs[[variable]][held_original] <- NA
  }
  list(ssn = ssn, held_original = held_original)
}

#' Build known \code{*_initial()} objects from a fitted model's own covariance parameters
#'
#' Used by local k-fold refitting to fix covariance (and, for GLMs,
#' dispersion) parameters at the full fit's estimates while refitting only
#' the fixed effects per fold.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param glm Whether to also build a known \code{dispersion_initial}
#'   object.
#'
#' @return A list with \code{tailup_type}, \code{taildown_type},
#'   \code{euclid_type}, \code{nugget_type} and their corresponding
#'   \code{*_initial} objects (all known), \code{randcov_initial}
#'   (\code{NULL} if there are no random effects), and, if \code{glm},
#'   \code{dispersion_initial}.
#'
#' @noRd
get_kcv_known_initials <- function(object, glm = FALSE) {
  params <- object$coefficients$params_object
  tailup_type <- remove_covtype(class(params$tailup)[1])
  taildown_type <- remove_covtype(class(params$taildown)[1])
  euclid_type <- remove_covtype(class(params$euclid)[1])
  nugget_type <- remove_covtype(class(params$nugget)[1])
  output <- list(
    tailup_type = tailup_type,
    taildown_type = taildown_type,
    euclid_type = euclid_type,
    nugget_type = nugget_type,
    tailup_initial = do.call(
      tailup_initial,
      c(list(tailup_type = tailup_type), as.list(params$tailup), list(known = "given"))
    ),
    taildown_initial = do.call(
      taildown_initial,
      c(list(taildown_type = taildown_type), as.list(params$taildown), list(known = "given"))
    ),
    euclid_initial = do.call(
      euclid_initial,
      c(list(euclid_type = euclid_type), as.list(params$euclid), list(known = "given"))
    ),
    nugget_initial = do.call(
      nugget_initial,
      c(list(nugget_type = nugget_type), as.list(params$nugget), list(known = "given"))
    )
  )
  output$randcov_initial <- if (is.null(object$random)) {
    NULL
  } else {
    do.call(randcov_initial, c(as.list(params$randcov), list(known = "given")))
  }
  if (glm) {
    output$dispersion_initial <- dispersion_initial(
      object$family, dispersion = as.vector(params$dispersion), known = "given"
    )
  }
  output
}

#' Build the \code{local} control list for one fold's local refit
#'
#' Reuses the fitted model's own local fitting groups (via
#' \code{object$local_index}, restricted to the retained rows) when present,
#' otherwise builds fresh k-means groups sized to \code{local$size} (capped
#' at the number of retained rows). Always uses theoretical variance
#' adjustment and forwards parallel controls to fitting.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param fold A vector of row indices (into the fitted rows) held out for
#'   this fold.
#' @param local A resolved \code{local} list (uses \code{local$size}).
#'
#' @return A \code{local} control list for the retained-fold refit.
#'
#' @noRd
get_kcv_estimation_local <- function(object, fold, local) {
  fit_local <- list(var_adjust = "theoretical", parallel = local$parallel)
  if (isTRUE(local$parallel)) fit_local$ncores <- local$ncores
  if (!is.null(object$local_index)) {
    fit_local$index <- object$local_index[-fold]
  } else {
    fit_local$method <- "kmeans"
    fit_local$size <- min(local$size, object$n - length(fold))
  }
  fit_local
}

#' Extract and align one fold's held-out predictions from a refit's \code{predict()} output
#'
#' @param prediction The output of \code{predict()} on a fold's local refit
#'   (a named vector, or a list with \code{fit}/\code{se.fit} when
#'   \code{se.fit} is requested).
#' @param held_original The fold's row positions in the original data, used
#'   to match against \code{prediction}'s \code{"NetworkID/pid"} names.
#' @param se.fit Whether \code{prediction} includes standard errors.
#'
#' @return A list with \code{pred} and, if \code{se.fit}, \code{se.fit},
#'   each aligned to \code{held_original}'s order.
#'
#' @noRd
extract_kcv_fold_prediction <- function(prediction, held_original, se.fit) {
  fit <- if (se.fit) prediction$fit else prediction
  position <- match(as.character(held_original), names(fit))
  if (anyNA(position)) {
    stop("The fold predictions could not be aligned to the original observation rows.", call. = FALSE)
  }
  output <- list(pred = as.numeric(fit[position]))
  if (se.fit) {
    se_position <- match(as.character(held_original), names(prediction$se.fit))
    if (anyNA(se_position)) {
      stop("The fold prediction standard errors could not be aligned to the original observation rows.", call. = FALSE)
    }
    output$se.fit <- as.numeric(prediction$se.fit[se_position])
  }
  output
}

#' Refit an \code{ssn_lm} model with one fold held out, for local k-fold cross-validation
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param fold A vector of row indices (into the fitted rows) to hold out.
#' @param local A resolved \code{local} list.
#'
#' @return The refit \code{ssn_lm} model object, with the fold's rows as its
#'   \code{".missing"} prediction set and covariance parameters fixed at the
#'   original fit's estimates.
#'
#' @noRd
get_kcv_local_lm_refit <- function(object, fold, local) {
  fold_data <- get_kcv_fold_ssn(object, fold)
  initial <- get_kcv_known_initials(object)
  do.call(
    ssn_lm,
    c(
      list(
        formula = object$formula, ssn.object = fold_data$ssn,
        tailup_type = initial$tailup_type, taildown_type = initial$taildown_type,
        euclid_type = initial$euclid_type, nugget_type = initial$nugget_type,
        tailup_initial = initial$tailup_initial, taildown_initial = initial$taildown_initial,
        euclid_initial = initial$euclid_initial, nugget_initial = initial$nugget_initial,
        additive = object$additive, estmethod = object$estmethod,
        anisotropy = object$anisotropy, random = object$random,
        randcov_initial = initial$randcov_initial,
        partition_factor = object$partition_factor,
        local = get_kcv_estimation_local(object, fold, local), contrasts = object$contrasts
      )
    )
  )
}

#' Compute local (big-data) k-fold cross-validation predictions for a Gaussian model
#'
#' For each fold, refits via \code{\link{get_kcv_local_lm_refit}()} and
#' predicts the held-out rows, extracting/aligning the result via
#' \code{\link{extract_kcv_fold_prediction}()}.
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param fold_list A named list of per-fold row-index vectors.
#' @param se.fit Whether to also return prediction standard errors.
#' @param local A resolved \code{local} list.
#'
#' @return A list with \code{cv_predict} and \code{se.fit} (\code{NULL} if
#'   not requested).
#'
#' @noRd
get_kcv_local.ssn_lm <- function(object, fold_list, se.fit, local) {
  pred <- numeric(object$n)
  if (se.fit) se <- numeric(object$n)
  for (i in seq_along(fold_list)) {
    held <- fold_list[[i]]
    refit <- get_kcv_local_lm_refit(object, held, local)
    prediction <- predict(
      refit, newdata = ".missing", se.fit = se.fit, interval = "none", local = local
    )
    fold_prediction <- extract_kcv_fold_prediction(
      prediction, object$observed_index[held], se.fit
    )
    pred[held] <- fold_prediction$pred
    if (se.fit) se[held] <- fold_prediction$se.fit
  }
  list(cv_predict = pred, se.fit = if (se.fit) se else NULL)
}
