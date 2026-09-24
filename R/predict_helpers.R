
#' Pair each newdata row with its model-matrix row and observed-covariance vector
#'
#' Shared by \code{predict.ssn_lm()} and \code{predict.ssn_glm()}: builds the
#' per-row list \code{get_pred()}/\code{get_pred_glm()} iterate over (in serial
#' or parallel dispatch alike).
#'
#' @param newdata A prediction data frame.
#' @param newdata_model The corresponding prediction model matrix.
#' @param cov_vector_list A list of observed-covariance vectors, one per row
#'   of \code{newdata} (from \code{\link{get_point_pred_cov_vector_list}()}).
#'
#' @return A list of length \code{NROW(newdata)}, each element a list with
#'   \code{row} (the newdata row), \code{x0} (the model matrix row), and
#'   \code{c0} (the covariance vector).
#'
#' @noRd
get_newdata_pred_list <- function(newdata, newdata_model, cov_vector_list) {
  newdata_rows_list <- split(newdata, seq_len(NROW(newdata)))
  newdata_model_list <- split(newdata_model, seq_len(NROW(newdata)))
  mapply(
    x = newdata_rows_list, y = newdata_model_list, c = cov_vector_list,
    FUN = function(x, y, c) list(row = x, x0 = y, c0 = c), SIMPLIFY = FALSE
  )
}

#' Resolve a requested prediction set name to every stored prediction set when omitted
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The requested prediction set name, \code{"all"}, or
#'   \code{NULL}.
#'
#' @return \code{newdata_name} unchanged, or (if \code{NULL}/\code{"all"})
#'   every name in \code{object$ssn.object$preds}.
#'
#' @noRd
resolve_newdata_name <- function(object, newdata_name) {
  if (is.null(newdata_name) || identical(newdata_name, "all")) {
    newdata_name <- names(object$ssn.object$preds)
  }
  newdata_name
}

#' Fetch a fitted model's observed data and one named prediction set
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of one prediction set in
#'   \code{object$ssn.object$preds}.
#'
#' @return A list with \code{obdata}, \code{newdata}, and
#'   \code{add_newdata_rows} (whether \code{newdata_name} is
#'   \code{".missing"}, the fitted model's own missing-response rows).
#'
#' @noRd
get_prediction_newdata <- function(object, newdata_name) {
  list(
    obdata = object$ssn.object$obs,
    newdata = object$ssn.object$preds[[newdata_name]],
    add_newdata_rows = identical(newdata_name, ".missing")
  )
}

#' Build a prediction-set model matrix, aligned to the fitted model's own columns
#'
#' Builds \code{newdata}'s model matrix using the fitted model's terms,
#' factor levels, and contrasts, then restricts it to the columns present in
#' the fitted model's own design matrix. Works around a \code{model.frame()}
#' bug with a degree-2 orthogonal polynomial term
#' (e.g. \code{poly(x, y, degree = 2)}) and exactly one prediction row by
#' duplicating that row before building the frame/matrix, then taking only
#' the first row of the result.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata A prediction data frame.
#'
#' @return A list with \code{newdata} (unchanged, except restored to one row
#'   if the workaround applied), \code{newdata_model} (the prediction design
#'   matrix), and \code{offset}.
#'
#' @noRd
get_newdata_model_matrix <- function(object, newdata) {
  formula_newdata <- delete.response(terms(object))
  # fix model frame bug with degree 2 basic polynomial and one prediction row
  # e.g. poly(x, y, degree = 2) and newdata has one row
  if (any(grepl("nmatrix.", attributes(formula_newdata)$dataClasses, fixed = TRUE)) &&
    NROW(newdata) == 1) {
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

  if (any(!complete.cases(newdata_model))) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }

  list(newdata = newdata, newdata_model = newdata_model, offset = offset)
}

#' Build the observed-covariance list for point-level prediction
#'
#' Prefers \code{.bmat}, matching block prediction, and passes the selected
#' backend to each chunk reader. Indexed reads limit each chunk's distance
#' input and covariance construction, but all covariance vectors are retained:
#' total covariance storage remains \code{O(n_obs * n_pred)}.
#'
#' Dense \code{.RData} files require whole-matrix deserialization, so dense
#' fallback uses one \code{covmatrix()} call to avoid repeated full-file reads.
#' Models without stream covariance also use that unchunked path.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set.
#' @param chunk_size The maximum number of newdata rows processed per chunk.
#'   Ignored when the resolved backend is not \code{"bigdata"}.
#' @param randcov_context A precomputed random-effect training-side context
#'   from \code{\link{get_randcov_context}()}, or \code{NULL} (built here).
#'
#' @return A list of length \code{NROW(object$ssn.object$preds[[newdata_name]])},
#'   each element a length-\code{object$n} numeric covariance vector, in
#'   \code{newdata_name}'s original row order.
#'
#' @details Also builds a partition-factor context (when the model has one)
#'   once here and passes it to every chunk's \code{\link{get_block_obs_covariance}()}
#'   call, the same way \code{randcov_context} is reused -- neither is built
#'   by (or forwarded to) the dense/unchunked early-return path above, since
#'   that path makes only one covariance call and has nothing to reuse across.
#'
#' @noRd
get_point_pred_cov_vector_list <- function(object, newdata_name, chunk_size, randcov_context = NULL) {
  backend <- get_block_pred_backend(object, newdata_name, prefer = "bigdata")
  if (!identical(backend, "bigdata")) {
    cov_vector <- covmatrix(object, newdata_name)
    cov_vector_list <- split(cov_vector, seq_len(NROW(cov_vector)))
    names(cov_vector_list) <- NULL
    return(cov_vector_list)
  }

  n_pred <- NROW(object$ssn.object$preds[[newdata_name]])
  chunks <- get_prediction_chunks(seq_len(n_pred), chunk_size)
  if (is.null(randcov_context)) {
    # built once and reused for every chunk, instead of re-deriving the
    # training-side random-effect grouping on each call
    randcov_context <- get_randcov_context(
      object$coefficients$params_object$randcov, object$ssn.object$obs,
      object$ssn.object$preds[[newdata_name]]
    )
  }
  # same idea as randcov_context above, for the partition factor: built once
  # (using every chunk's rows, so its extended factor levels already cover
  # what any individual chunk needs) instead of every chunk's
  # get_block_obs_covariance() -> get_cov_vector() -> partition_vector() call
  # re-deriving the training-side grouping from scratch
  partition_context <- if (is.null(object$partition_factor)) {
    NULL
  } else {
    get_partition_context(
      object$partition_factor, object$ssn.object$obs,
      object$ssn.object$preds[[newdata_name]]
    )
  }
  cov_vector_list <- unlist(lapply(chunks, function(rows) {
    chunk_cov <- get_block_obs_covariance(
      object, newdata_name, rows,
      backend = backend, randcov_context = randcov_context, partition_context = partition_context
    )
    split(chunk_cov, seq_len(NROW(chunk_cov)))
  }), recursive = FALSE)
  names(cov_vector_list) <- NULL
  cov_vector_list
}

#' Assemble point prediction's per-operation shared context, once
#'
#' Bundles the setup common to \code{predict.ssn_lm()}/\code{predict.ssn_glm()}
#' -- the observed-covariance list, per-row prediction list, full observed
#' covariance matrix (when \code{local$method == "all"}), spatial/nugget
#' marginal variance, random-effect variance parameters and their reusable
#' training-side context, and the observed design matrix/response/offset --
#' so it is built once per prediction call by a shared helper instead of
#' independently re-derived by the two predict methods.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set.
#' @param newdata The (already model-matrix-validated) prediction data frame.
#' @param local_list A resolved local-prediction settings list (see
#'   \code{\link{get_local_list_prediction}()}).
#'
#' @return A list with \code{cov_vector_list}, \code{newdata_list},
#'   \code{cov_matrix_val}, \code{spatial_nugget_var}, \code{randcov_params},
#'   \code{cov_lowchol}, \code{randcov_context}, \code{Xmat}, \code{y}, and
#'   \code{offset}.
#'
#' @noRd
get_point_pred_context <- function(object, newdata_name, newdata, newdata_model, local_list) {
  params_object <- object$coefficients$params_object
  randcov_params <- params_object$randcov
  # built once here and reused both by get_point_pred_cov_vector_list()'s
  # chunk loop and by every row's randcov_newvar() call below, instead of
  # each independently re-deriving its own copy
  randcov_context <- get_randcov_context(randcov_params, object$ssn.object$obs, newdata)

  cov_vector_list <- get_point_pred_cov_vector_list(object, newdata_name, local_list$chunk_size, randcov_context = randcov_context)
  newdata_list <- get_newdata_pred_list(newdata, newdata_model, cov_vector_list)
  cov_matrix_val <- covmatrix(object)
  spatial_nugget_var <- get_spatial_nugget_var(params_object, object$diagtol)
  cov_lowchol <- if (local_list$method == "all") t(chol(cov_matrix_val)) else NULL

  list(
    cov_vector_list = cov_vector_list, newdata_list = newdata_list,
    cov_matrix_val = cov_matrix_val, spatial_nugget_var = spatial_nugget_var,
    randcov_params = randcov_params, cov_lowchol = cov_lowchol,
    randcov_context = randcov_context,
    Xmat = model.matrix(object), y = model.response(model.frame(object)),
    offset = model.offset(model.frame(object))
  )
}

#' Dispatch point prediction or cross-validation folds, in parallel or serially
#'
#' Uses serial or parallel evaluation with cleanup on success and error.
#' Workers load the same development source or installation as the parent
#' so namespace-resolved helpers remain consistent across processes.
#'
#' @param fn The prediction or cross-validation function.
#' @param data_list A list or vector of prediction rows or held-out indices.
#' @param local_list A resolved local-prediction settings list;
#'   \code{local_list$parallel} selects \code{parallel::parLapply()} over
#'   \code{local_list$ncores} workers, otherwise plain \code{lapply()}.
#' @param ... Every other named argument \code{fn} needs, forwarded unchanged.
#'
#' @return A list with one element per input, as returned by
#'   \code{fn}.
#'
#' @noRd
run_pred_dispatch <- function(fn, data_list, local_list, ...) {
  if (!local_list$parallel) {
    return(lapply(data_list, fn, ...))
  }
  cl <- make_ssn_cluster(local_list$ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::parLapply(cl, data_list, fn, ...)
}

#' Create workers using the parent's SSN2 source or installation
#'
#' @param ncores Number of workers.
#' @return A cluster that the caller must stop after use.
#' @noRd
make_ssn_cluster <- function(ncores) {
  cl <- parallel::makeCluster(ncores)
  ready <- FALSE
  on.exit(if (!ready) parallel::stopCluster(cl), add = TRUE)
  dev_path <- NULL
  if (requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("SSN2")) {
    source_file <- utils::getSrcFilename(utils::getSrcref(make_ssn_cluster), full.names = TRUE)
    candidate <- normalizePath(file.path(dirname(source_file), ".."), mustWork = TRUE)
    if (file.exists(file.path(candidate, "DESCRIPTION"))) dev_path <- candidate
  }
  parallel::clusterCall(cl, function(dev_path) {
    if (is.null(dev_path)) {
      loadNamespace("SSN2")
    } else {
      pkgload::load_all(dev_path, quiet = TRUE)
    }
    NULL
  }, dev_path)
  ready <- TRUE
  cl
}

#' Apply row naming and pick the return shape for an \code{interval = "none"} prediction
#'
#' The last step of \code{predict.ssn_lm()}'s/\code{predict.ssn_glm()}'s
#' \code{interval == "none"} branch is identical once \code{fit}/\code{se} are
#' on their final scale (post-offset, post-\code{invlink()}/delta-method
#' adjustment where applicable) -- name the rows if requested and return
#' either \code{fit} alone or \code{list(fit, se.fit)}, depending on whether
#' standard errors were requested at all.
#'
#' @param fit The point predictions, already on their final scale.
#' @param se The standard errors, already on their final scale, or
#'   \code{NULL} when standard errors were not requested -- this
#'   \code{NULL}-ness is itself the signal for which return shape to use.
#' @param add_newdata_rows Whether to name the returned values using
#'   \code{missing_index}.
#' @param missing_index Row labels to apply (\code{object$missing_index}).
#'
#' @return \code{list(fit = fit, se.fit = se)} if \code{se} was supplied,
#'   otherwise \code{fit} alone.
#'
#' @noRd
finalize_interval_none <- function(fit, se, add_newdata_rows, missing_index) {
  if (!is.null(se)) {
    if (add_newdata_rows) {
      names(fit) <- missing_index
      names(se) <- missing_index
    }
    list(fit = fit, se.fit = se)
  } else {
    if (add_newdata_rows) {
      names(fit) <- missing_index
    }
    fit
  }
}

#' Assemble the fit/lwr/upr matrix and pick the return shape for an interval prediction
#'
#' The last step of \code{predict.ssn_lm()}'s/\code{predict.ssn_glm()}'s
#' \code{interval == "confidence"}/\code{"prediction"} branches is identical
#' once \code{fit}/\code{lwr}/\code{upr}/\code{se} are on their final scale
#' (post-offset, post-\code{invlink()}/delta-method adjustment where
#' applicable) -- bind them into the returned matrix, name the rows if
#' requested, and return either the matrix alone or \code{list(fit, se.fit)}.
#' Unlike \code{finalize_interval_none()}, \code{se} is always a real vector
#' here (needed upstream to build \code{lwr}/\code{upr} regardless of whether
#' standard errors were actually requested), so \code{se.fit} must be passed
#' explicitly rather than inferred from \code{se}'s \code{NULL}-ness.
#'
#' @param fit,lwr,upr The point predictions and interval bounds, already on
#'   their final scale.
#' @param se The standard errors, already on their final scale.
#' @param se.fit Whether standard errors were requested.
#' @param add_newdata_rows Whether to name the returned values using
#'   \code{missing_index}.
#' @param missing_index Row labels to apply (\code{object$missing_index}).
#'
#' @return \code{list(fit = <n x 3 matrix>, se.fit = se)} if \code{se.fit},
#'   otherwise the \code{<n x 3 matrix>} alone.
#'
#' @noRd
finalize_interval_bounds <- function(fit, lwr, upr, se, se.fit, add_newdata_rows, missing_index) {
  fit <- cbind(fit, lwr, upr)
  row.names(fit) <- seq_len(NROW(fit))
  if (se.fit) {
    if (add_newdata_rows) {
      row.names(fit) <- missing_index
      names(se) <- missing_index
    }
    list(fit = fit, se.fit = se)
  } else {
    if (add_newdata_rows) {
      row.names(fit) <- missing_index
    }
    fit
  }
}
