#' Package fitted coefficients into a fitted model's \code{coefficients} element
#'
#' @param betahat The fixed-effect coefficient estimates.
#' @param params_object A joint covariance parameter object.
#'
#' @return A list with \code{fixed} and \code{params_object}.
#'
#' @noRd
get_coefficients <- function(betahat, params_object) {
  # store fixed effect and random coefficients as list
  list(fixed = betahat, params_object = params_object)
}

#' Build one covariance component's per-group covariance matrix list
#'
#' @param params A single covariance component's parameter object (e.g.
#'   \code{params_object$tailup}).
#' @param data_object A model data object with \code{dist_object_oblist},
#'   \code{anisotropy}, and (optionally) \code{partition_list}.
#'
#' @return A list of covariance matrices, one per group in
#'   \code{data_object$dist_object_oblist}, with the partition factor applied
#'   if present.
#'
#' @noRd
get_component_covariance_list <- function(params, data_object) {
  covariance <- lapply(data_object$dist_object_oblist, function(x) {
    cov_matrix(params, x, anisotropy = data_object$anisotropy)
  })
  if (!is.null(data_object$partition_list)) {
    covariance <- Map(`*`, covariance, data_object$partition_list)
  }
  covariance
}

#' Compute fitted values, decomposed by covariance component
#'
#' Computes the overall mean-structure fitted response (plus offset), and
#' each active covariance component's (tailup, taildown, euclid, nugget,
#' random effects) contribution to the fitted values via its covariance
#' matrix applied to the whitened residual.
#'
#' @param betahat The fixed-effect coefficient estimates.
#' @param params_object A joint covariance parameter object.
#' @param data_object A model data object.
#' @param eigenprods_list A list of per-group eigen-decomposition products
#'   (with \code{SigInv_y}, \code{SigInv_X}), one per group.
#'
#' @return A list with \code{response}, \code{tailup}, \code{taildown},
#'   \code{euclid}, \code{nugget} (each \code{NULL} if that component is
#'   inactive), and \code{randcov} (a named list, one element per random
#'   effect, or \code{NULL}).
#'
#' @noRd
get_fitted <- function(betahat, params_object, data_object, eigenprods_list) {
  # find mean fitted values
  fitted_response <- as.numeric(do.call("rbind", lapply(data_object$X_list, function(x) x %*% betahat)))

  # incorporate offset if necessary
  if (!is.null(data_object$offset)) {
    fitted_response <- fitted_response + data_object$offset
  }

  # find SigInv times residuals product used throughout
  SigInv_r_list <- lapply(eigenprods_list, function(x) x$SigInv_y - x$SigInv_X %*% betahat)

  # find tailup fitted values (NULL if not used)
  tailup_none <- inherits(params_object$tailup, "tailup_none")
  if (tailup_none) {
    fitted_tailup <- NULL
  } else {
    tailup_list <- get_component_covariance_list(params_object$tailup, data_object)
    fitted_tailup <- as.numeric(do.call("rbind", mapply(
      s = tailup_list, r = SigInv_r_list,
      function(s, r) s %*% r, SIMPLIFY = FALSE
    )))
  }

  # find taildown fitted values (NULL if not used)
  taildown_none <- inherits(params_object$taildown, "taildown_none")
  if (taildown_none) {
    fitted_taildown <- NULL
  } else {
    taildown_list <- get_component_covariance_list(params_object$taildown, data_object)
    fitted_taildown <- as.numeric(do.call("rbind", mapply(
      s = taildown_list, r = SigInv_r_list,
      function(s, r) s %*% r, SIMPLIFY = FALSE
    )))
  }

  # find euclid fitted values (NULL if not used)
  euclid_none <- inherits(params_object$euclid, "euclid_none")
  if (euclid_none) {
    fitted_euclid <- NULL
  } else {
    euclid_list <- get_component_covariance_list(params_object$euclid, data_object)
    fitted_euclid <- as.numeric(do.call("rbind", mapply(
      s = euclid_list, r = SigInv_r_list,
      function(s, r) s %*% r, SIMPLIFY = FALSE
    )))
  }

  # find nugget fitted values (NULL if not used)
  nugget_none <- inherits(params_object$nugget, "nugget_none")
  if (nugget_none) {
    fitted_nugget <- NULL
  } else {
    fitted_nugget <- as.numeric(params_object$nugget[["nugget"]] * do.call("rbind", SigInv_r_list))
  }

  # find random effect fitted values (NULL if not used)
  if (is.null(names(params_object$randcov))) {
    fitted_randcov <- NULL
  } else {
    fitted_randcov <- lapply(names(params_object$randcov), function(x) {
      fitted_val <- params_object$randcov[[x]] * do.call("rbind", mapply(
        z = data_object$randcov_list,
        r = SigInv_r_list,
        function(z, r) {
          crossprod(z[[x]][["Z"]], r)
        }
      ))
      fitted_val <- tapply(fitted_val, rownames(fitted_val), function(x) {
        val <- mean(x[x != 0])
        if (length(val) == 0) { # replace if all zeros somehow
          val <- rep(0, length(x))
          names(val) <- names(x)
        }
        val
      })
      # all combinations yields values with many zeros -- don't want to include these in the mean
      names_fitted_val <- rownames(fitted_val)
      fitted_val <- as.numeric(fitted_val)
      names(fitted_val) <- names_fitted_val
      fitted_val
    })
    names(fitted_randcov) <- names(params_object$randcov)
  }

  # return all as list
  fitted_values <- list(
    response = as.numeric(fitted_response),
    tailup = as.numeric(fitted_tailup),
    taildown = as.numeric(fitted_taildown),
    euclid = as.numeric(fitted_euclid),
    nugget = as.numeric(fitted_nugget),
    randcov = fitted_randcov
  )
}

#' Compute leverage (hat) values on the whitened scale
#'
#' @param cov_betahat The fixed-effect coefficient covariance matrix.
#' @param SqrtSigInv_X The whitened design matrix
#'   (\eqn{\Sigma^{-1/2}}\code{X}).
#'
#' @return A numeric vector of hat values (diagonal of the whitened hat
#'   matrix).
#'
#' @noRd
get_hatvalues <- function(cov_betahat, SqrtSigInv_X) {
  # only the diagonal of the whitened hat matrix is needed, so skip forming
  # the full n x n product
  get_diag_XVXt(SqrtSigInv_X, cov_betahat)
}

#' Compute response, Pearson, and standardized residuals
#'
#' @param betahat The fixed-effect coefficient estimates.
#' @param data_object A model data object with \code{X_list}/\code{y_list}.
#' @param eigenprods_list A list of per-group eigen-decomposition products
#'   (with \code{SqrtSigInv_y}, \code{SqrtSigInv_X}), one per group.
#' @param hatvalues Leverage values from \code{\link{get_hatvalues}()}.
#'
#' @return A list with \code{response}, \code{pearson}, and
#'   \code{standardized} residual vectors.
#'
#' @noRd
get_residuals <- function(betahat, data_object, eigenprods_list, hatvalues) {
  # first find response residuals
  residuals_response <- as.numeric(do.call("rbind", mapply(
    y = data_object$y_list, x = data_object$X_list,
    function(y, x) y - x %*% betahat, SIMPLIFY = FALSE
  )))

  # then find pearson residuals (pre multiplied by inverse square root)
  residuals_pearson <- as.numeric(do.call(
    "rbind",
    lapply(eigenprods_list, function(x) x$SqrtSigInv_y - x$SqrtSigInv_X %*% betahat)
  ))

  # then find standardized
  residuals_standardized <- residuals_pearson / sqrt(1 - hatvalues) # (I - H on bottom)

  # return as list
  list(response = as.numeric(residuals_response), pearson = as.numeric(residuals_pearson), standardized = as.numeric(residuals_standardized))
}

#' Compute Cook's distance
#'
#' @param residuals A residuals list from \code{\link{get_residuals}()} (uses
#'   \code{residuals$standardized}).
#' @param hatvalues Leverage values from \code{\link{get_hatvalues}()}.
#' @param p The number of fixed effects.
#'
#' @return A numeric vector of Cook's distances.
#'
#' @noRd
get_cooks_distance <- function(residuals, hatvalues, p) {
  # find cook's distance
  residuals$standardized^2 * hatvalues / (p * (1 - hatvalues))
}

#' Assemble fitted-model statistics for an IID (no spatial dependence, no
#' random effects) fit, exact or grouped/big-data
#'
#' Shared by \code{get_model_stats_iid()} and
#' \code{get_model_stats_bigdata_iid()}, which are identical apart from which
#' of \code{data_object$order}/\code{data_object$order_bigdata} restores
#' original row order; \code{order} takes the place of whichever one applies.
#'
#' @param cov_est_object A fitted covariance estimation object, as returned by
#'   \code{get_gloglik_iid()} (or \code{use_gloglik_known()} for a fully
#'   fixed nugget)
#' @param data_object The data object
#' @param order The row-order vector to restore original data order with
#'   (\code{data_object$order} for the exact case,
#'   \code{data_object$order_bigdata} for the grouped/big-data case).
#'
#' @return The same statistics as \code{get_model_stats()}/
#'   \code{get_model_stats_bigdata()}, computed directly from a QR
#'   decomposition of \code{X} rather than a full covariance matrix, since
#'   with no spatial dependence or random effects the covariance matrix is
#'   simply \code{nugget * I}
#'
#' @noRd
get_model_stats_iid_core <- function(cov_est_object, data_object, order) {
  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)

  # with Sigma = nugget * I, GLS reduces to OLS, so betahat and its
  # covariance come directly from the QR decomposition of X rather than a
  # full n x n covariance matrix
  qr_val <- qr(X)
  R_val <- qr.R(qr_val)
  nugget <- cov_est_object$params_object$nugget[["nugget"]]
  cor_betahat <- chol2inv(chol(crossprod(R_val, R_val)))
  cov_betahat <- nugget * cor_betahat
  betahat <- as.numeric(backsolve(R_val, qr.qty(qr_val, y)))
  names(betahat) <- colnames(data_object$X_list[[1]])
  cov_betahat <- as.matrix(cov_betahat)
  rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat) <- colnames(data_object$X_list[[1]])

  # return coefficients
  coefficients <- get_coefficients(betahat, cov_est_object$params_object)

  # return fitted
  fitted_response <- as.numeric(X %*% betahat)
  resid <- as.numeric(y - X %*% betahat)
  if (!is.null(data_object$offset)) {
    fitted_response <- fitted_response + data_object$offset
  }
  # tailup/taildown/euclid are inactive, matching the general path's own NA
  # placeholder for a "none" component (see get_fitted())
  fitted <- list(
    response = fitted_response,
    tailup = rep(NA_real_, data_object$n),
    taildown = rep(NA_real_, data_object$n),
    euclid = rep(NA_real_, data_object$n),
    nugget = resid,
    randcov = NULL
  )

  # return hat values (only the diagonal is needed)
  hatvalues <- get_diag_XVXt(X, cor_betahat)

  # return residuals
  residuals <- list(
    response = resid,
    pearson = resid / sqrt(nugget)
  )
  residuals$standardized <- residuals$pearson / sqrt(1 - hatvalues)

  # return cooks distance
  cooks_distance <- get_cooks_distance(residuals, hatvalues, data_object$p)

  # local estimation partitions the data, so all these vectors are currently
  # ordered by partition rather than original row order; order records how
  # to map partition order back to the original data order
  # reorder relevant quantities
  model_stats_names <- data_object$pid[data_object$observed_index]
  fitted$response <- fitted$response[order(order)]
  names(fitted$response) <- model_stats_names
  fitted$tailup <- fitted$tailup[order(order)]
  names(fitted$tailup) <- model_stats_names
  fitted$taildown <- fitted$taildown[order(order)]
  names(fitted$taildown) <- model_stats_names
  fitted$euclid <- fitted$euclid[order(order)]
  names(fitted$euclid) <- model_stats_names
  fitted$nugget <- fitted$nugget[order(order)]
  names(fitted$nugget) <- model_stats_names
  hatvalues <- hatvalues[order(order)]
  names(hatvalues) <- model_stats_names
  residuals$response <- residuals$response[order(order)]
  names(residuals$response) <- model_stats_names
  residuals$pearson <- residuals$pearson[order(order)]
  names(residuals$pearson) <- model_stats_names
  residuals$standardized <- residuals$standardized[order(order)]
  names(residuals$standardized) <- model_stats_names
  cooks_distance <- cooks_distance[order(order)]
  names(cooks_distance) <- model_stats_names

  # return variance covariance matrices
  vcov <- get_vcov(cov_betahat)

  # return deviance
  deviance <- as.numeric(crossprod(residuals$pearson, residuals$pearson))

  # generalized r squared
  muhat <- mean(y)
  pearson_null <- as.numeric(y - muhat) / sqrt(nugget)
  deviance_null <- as.numeric(crossprod(pearson_null, pearson_null))
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)

  # set null model R2 equal to zero (no covariates)
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }

  # return npar (number of estimated covariance parameters)
  npar <- sum(unlist(lapply(cov_est_object$is_known, function(x) length(x) - sum(x))))

  # return list
  list(
    coefficients = coefficients,
    fitted = fitted,
    hatvalues = hatvalues,
    residuals = residuals,
    cooks_distance = cooks_distance,
    vcov = vcov,
    deviance = deviance,
    pseudoR2 = pseudoR2,
    npar = npar
  )
}

#' Assemble fitted-model statistics for a general (non-IID) Gaussian fit,
#' exact or grouped/big-data
#'
#' Shared by \code{get_model_stats()} and \code{get_model_stats_bigdata()},
#' which are identical apart from which of
#' \code{data_object$order}/\code{data_object$order_bigdata} restores original
#' row order; \code{order} takes the place of whichever one applies.
#'
#' @param cov_est_object A fitted covariance estimation object.
#' @param data_object The data object.
#' @param order The row-order vector to restore original data order with
#'   (\code{data_object$order} for the exact case,
#'   \code{data_object$order_bigdata} for the grouped/big-data case).
#'
#' @return The same statistics list \code{get_model_stats()}/
#'   \code{get_model_stats_bigdata()} return.
#'
#' @noRd
get_model_stats_core <- function(cov_est_object, data_object, order) {
  # store the covariance matrix list
  cov_matrix_list <- get_cov_matrix_list(cov_est_object$params_object, data_object)

  # compute eigenproducts for later use
  if (data_object$parallel) {
    cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
      cluster_list_element <- list(
        c = cov_matrix_list[[l]],
        x = data_object$X_list[[l]],
        y = data_object$y_list[[l]],
        o = data_object$ones_list[[l]]
      )
    })
    eigenprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_eigenprods_parallel)
    names(eigenprods_list) <- names(cov_matrix_list)
  } else {
    eigenprods_list <- mapply(
      c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, o = data_object$ones_list,
      function(c, x, y, o) get_eigenprods(c, x, y, o),
      SIMPLIFY = FALSE
    )
  }

  # get inverse cov beta hat list and add together
  invcov_betahat_list <- lapply(eigenprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_X))
  invcov_betahat_sum <- Reduce("+", invcov_betahat_list)
  # find unadjusted cov beta hat matrix
  cov_betahat_noadjust <- chol2inv(chol(forceSymmetric(invcov_betahat_sum)))
  # put it in a list the number of times there are unique local indices
  cov_betahat_noadjust_list <- rep(list(cov_betahat_noadjust), times = length(invcov_betahat_list))

  # get relevant necessary product list-wise
  Xt_SigInv_y_list <- lapply(eigenprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_y))

  # get betahat list-wise
  betahat_list <- mapply(
    l = cov_betahat_noadjust_list, r = Xt_SigInv_y_list,
    function(l, r) l %*% r,
    SIMPLIFY = FALSE
  )

  # get global beta hat and name
  betahat <- as.numeric(cov_betahat_noadjust %*%
    Reduce("+", Xt_SigInv_y_list))
  names(betahat) <- colnames(data_object$X_list[[1]])

  # Account for dependence between local fitting groups.
  cov_betahat <- cov_betahat_adjust(
    invcov_betahat_list,
    betahat_list, betahat,
    eigenprods_list, data_object,
    cov_est_object$params_object,
    cov_betahat_noadjust, data_object$var_adjust
  )

  # and name it
  cov_betahat <- as.matrix(cov_betahat)
  rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat) <- colnames(data_object$X_list[[1]])

  # return fixed and random coefficients
  coefficients <- get_coefficients(betahat, cov_est_object$params_object)

  # return fixed effect fitted values and random fitted values
  fitted <- get_fitted(betahat, cov_est_object$params_object, data_object, eigenprods_list)

  # return hat values (leverage)
  hatvalues <- as.numeric(unlist(lapply(eigenprods_list, function(x) get_hatvalues(cov_betahat_noadjust, x$SqrtSigInv_X))))

  # return residuals (response, pearson, standardized)
  residuals <- get_residuals(betahat, data_object, eigenprods_list, hatvalues)

  # return cooks distance (influence)
  cooks_distance <- get_cooks_distance(residuals, hatvalues, data_object$p)

  # reorder relevant quantities to match data order
  ## fitted values
  model_stats_names <- data_object$pid[data_object$observed_index]
  fitted$response <- fitted$response[order(order)]
  names(fitted$response) <- model_stats_names
  fitted$tailup <- fitted$tailup[order(order)]
  names(fitted$tailup) <- model_stats_names
  fitted$taildown <- fitted$taildown[order(order)]
  names(fitted$taildown) <- model_stats_names
  fitted$euclid <- fitted$euclid[order(order)]
  names(fitted$euclid) <- model_stats_names
  fitted$nugget <- fitted$nugget[order(order)]
  names(fitted$nugget) <- model_stats_names
  ## hat values
  hatvalues <- hatvalues[order(order)]
  names(hatvalues) <- model_stats_names
  ## residuals
  residuals$response <- residuals$response[order(order)]
  names(residuals$response) <- model_stats_names
  residuals$pearson <- residuals$pearson[order(order)]
  names(residuals$pearson) <- model_stats_names
  residuals$standardized <- residuals$standardized[order(order)]
  names(residuals$standardized) <- model_stats_names
  ## cook's distance
  cooks_distance <- cooks_distance[order(order)]
  names(cooks_distance) <- model_stats_names

  # get variance covariance matrices
  vcov <- get_vcov(cov_betahat)

  # get deviance
  deviance <- as.numeric(crossprod(residuals$pearson, residuals$pearson))

  # generalized r squared
  ## find covariance matrix of fixed effect in null (intercept-only) model
  SqrtSigInv_ones <- as.numeric(do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_ones)))
  cov_muhat <- 1 / crossprod(SqrtSigInv_ones, SqrtSigInv_ones)
  ## find intercept estimate
  SqrtSigInv_y <- do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_y))
  muhat <- as.vector(cov_muhat * crossprod(SqrtSigInv_ones, SqrtSigInv_y))
  ## compute relevant residuals
  SqrtSigInv_rmuhat <- as.numeric(do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_y - x$SqrtSigInv_ones * muhat)))
  ## compute deviance for null model
  deviance_null <- as.numeric(crossprod(SqrtSigInv_rmuhat, SqrtSigInv_rmuhat))
  ## get pseudoR2 as 1 - deviance ratio
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)
  ## if no covariances set pseudoR2 to 0
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }

  # return npar (number of estimated covariance parameters)
  npar <- sum(unlist(lapply(cov_est_object$is_known, function(x) length(x) - sum(x))))

  # return list
  list(
    coefficients = coefficients,
    fitted = fitted,
    hatvalues = hatvalues,
    residuals = residuals,
    cooks_distance = cooks_distance,
    vcov = vcov,
    deviance = deviance,
    pseudoR2 = pseudoR2,
    npar = npar
  )
}

#' Assemble fitted-model statistics for a general (non-IID) GLM fit, exact or
#' grouped/big-data
#'
#' Shared by \code{get_model_stats_glm()} and
#' \code{get_model_stats_bigdata_glm()}, which are identical apart from (1)
#' which of \code{data_object$order}/\code{data_object$order_bigdata} restores
#' original row order, and (2) how \code{eigenprods_list} is computed --
#' \code{get_model_stats_glm()} always builds it serially, while
#' \code{get_model_stats_bigdata_glm()} optionally parallelizes that step.
#' Both differences are resolved before this core is called: each wrapper
#' builds its own \code{eigenprods_list} with its own serial/parallel policy,
#' then passes it in along with the applicable \code{order}.
#'
#' @param cov_est_object A fitted covariance estimation object.
#' @param data_object The data object.
#' @param estmethod The estimation method.
#' @param order The row-order vector to restore original data order with
#'   (\code{data_object$order} for the exact case,
#'   \code{data_object$order_bigdata} for the grouped/big-data case).
#' @param eigenprods_list Each group's eigen-decomposition products, already
#'   computed by the caller (serially or in parallel).
#'
#' @return The same statistics list \code{get_model_stats_glm()}/
#'   \code{get_model_stats_bigdata_glm()} return.
#'
#' @noRd
get_model_stats_glm_core <- function(cov_est_object, data_object, estmethod, order, eigenprods_list) {
  # find model components
  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)

  SigInv_list <- lapply(eigenprods_list, function(x) x$SigInv)
  SigInv <- Matrix::bdiag(SigInv_list)
  SigInv_X <- do.call("rbind", lapply(eigenprods_list, function(x) x$SigInv_X))

  # get inverse cov beta hat list and add together
  invcov_betahat_list <- lapply(eigenprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_X))
  invcov_betahat_sum <- Reduce("+", invcov_betahat_list)
  # find unadjusted cov beta hat matrix
  cov_betahat_noadjust <- chol2inv(chol(forceSymmetric(invcov_betahat_sum)))
  # put it in a list the number of times there are unique local indices
  cov_betahat_noadjust_list <- rep(list(cov_betahat_noadjust), times = length(invcov_betahat_list))

  # dispersion
  dispersion <- as.vector(cov_est_object$params_object$dispersion)

  # Newton-Raphson
  w_and_H <- get_w_and_H(data_object, dispersion,
    SigInv_list, SigInv_X, cov_betahat_noadjust,
    invcov_betahat_sum, estmethod,
    ret_mHInv = TRUE
  )

  w <- as.vector(w_and_H$w)
  # H <- w_and_H$H

  # put w in eigenprods
  w_list <- split(w, sort(data_object$local_index))

  Xt_SigInv_w_list <- mapply(
    x = eigenprods_list, w = w_list,
    function(x, w) crossprod(x$SigInv_X, w),
    SIMPLIFY = FALSE
  )

  betahat_list <- mapply(
    l = cov_betahat_noadjust_list, r = Xt_SigInv_w_list,
    function(l, r) l %*% r,
    SIMPLIFY = FALSE
  )

  betahat <- as.numeric(cov_betahat_noadjust %*%
    Reduce("+", Xt_SigInv_w_list))
  names(betahat) <- colnames(data_object$X_list[[1]])

  # Account for dependence between local fitting groups.
  cov_betahat <- cov_betahat_adjust(
    invcov_betahat_list,
    betahat_list, betahat,
    eigenprods_list, data_object,
    cov_est_object$params_object,
    cov_betahat_noadjust, data_object$var_adjust
  )

  cov_betahat <- as.matrix(cov_betahat)
  wts_beta <- tcrossprod(cov_betahat, SigInv_X)
  betawtsvarw <- wts_beta %*% w_and_H$mHInv %*% t(wts_beta)

  cov_betahat_uncorrected <- cov_betahat # save uncorrected cov beta hat
  rownames(cov_betahat_uncorrected) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat_uncorrected) <- colnames(data_object$X_list[[1]])

  cov_betahat <- as.matrix(cov_betahat + betawtsvarw)
  rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
  colnames(cov_betahat) <- colnames(data_object$X_list[[1]])

  # return fixed and random coefficients
  coefficients <- get_coefficients_glm(betahat, cov_est_object$params_object)

  # return fitted
  fitted <- get_fitted_glm(w_list, betahat, cov_est_object$params_object, data_object, eigenprods_list)

  # return hat values
  hatvalues <- get_hatvalues_glm(w, X, data_object, dispersion)

  # return deviance i
  deviance_i <- get_deviance_glm(data_object$family, y, fitted$response, data_object$size, dispersion)
  deviance_i <- pmax(deviance_i, 0) # sometimes numerical instability can cause these to be slightly negative

  # storing relevant products
  SigInv_X_null <- do.call("rbind", lapply(eigenprods_list, function(x) x$SigInv_ones))
  ## lower chol %*% X
  SqrtSigInv_X_null <- do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_ones))
  # covariance of beta hat
  ## t(X) %*% sigma_inverse %*% X
  Xt_SigInv_X_null <- crossprod(SqrtSigInv_X_null, SqrtSigInv_X_null)
  ## t(X) %*% sigma_inverse %*% X)^(-1)
  Xt_SigInv_X_upchol_null <- chol(Xt_SigInv_X_null)
  cov_betahat_null <- chol2inv(Xt_SigInv_X_upchol_null)

  # Newton-Raphson
  w_and_H_null <- get_w_and_H(
    data_object, dispersion,
    SigInv_list, SigInv_X_null, cov_betahat_null, Xt_SigInv_X_null, estmethod
  )

  w_null <- as.vector(w_and_H_null$w)

  fitted_null <- get_fitted_null(w_null, data_object)

  # return deviance i
  deviance_i_null <- get_deviance_glm(data_object$family, y, fitted_null, data_object$size, dispersion)
  deviance_i_null <- pmax(deviance_i_null, 0) # sometimes numerical instability can cause these to be slightly non-negative

  deviance <- sum(deviance_i)
  deviance_null <- sum(deviance_i_null)
  pseudoR2 <- as.numeric(1 - deviance / deviance_null)

  # should always be non-negative
  pseudoR2 <- pmax(0, pseudoR2)
  # set null model R2 equal to zero (no covariates)
  if (length(labels(terms(data_object$formula))) == 0) {
    pseudoR2 <- 0
  }

  # return residuals
  residuals <- get_residuals_glm(w, y, data_object, deviance_i, hatvalues, dispersion)

  # return cooks distance
  cooks_distance <- get_cooks_distance_glm(residuals, hatvalues, data_object$p)

  # return variance covariance matrices
  vcov <- get_vcov_glm(cov_betahat, cov_betahat_uncorrected) # note first argument is the adjusted one

  # reorder relevant quantities to match data order
  model_stats_glm_names <- data_object$pid[data_object$observed_index]
  ## fitted values
  fitted$response <- fitted$response[order(order)]
  names(fitted$response) <- model_stats_glm_names
  fitted$link <- fitted$link[order(order)]
  names(fitted$link) <- model_stats_glm_names
  fitted$tailup <- fitted$tailup[order(order)]
  names(fitted$tailup) <- model_stats_glm_names
  fitted$taildown <- fitted$taildown[order(order)]
  names(fitted$taildown) <- model_stats_glm_names
  fitted$euclid <- fitted$euclid[order(order)]
  names(fitted$euclid) <- model_stats_glm_names
  fitted$nugget <- fitted$nugget[order(order)]
  names(fitted$nugget) <- model_stats_glm_names
  ## hat values
  hatvalues <- hatvalues[order(order)]
  names(hatvalues) <- model_stats_glm_names
  ## residuals
  residuals$response <- residuals$response[order(order)]
  names(residuals$response) <- model_stats_glm_names
  residuals$deviance <- residuals$deviance[order(order)]
  names(residuals$deviance) <- model_stats_glm_names
  residuals$pearson <- residuals$pearson[order(order)]
  names(residuals$pearson) <- model_stats_glm_names
  residuals$standardized <- residuals$standardized[order(order)]
  names(residuals$standardized) <- model_stats_glm_names
  ## cook's distance
  cooks_distance <- cooks_distance[order(order)]
  names(cooks_distance) <- model_stats_glm_names
  y <- y[order(order)]
  if (is.null(data_object$size)) {
    size <- NULL
  } else {
    size <- data_object$size[order(order)]
    names(size) <- model_stats_glm_names
  }

  # return npar (number of estimated covariance parameters)
  npar <- sum(unlist(lapply(cov_est_object$is_known, function(x) length(x) - sum(x))))

  # return list
  list(
    coefficients = coefficients,
    fitted = fitted,
    hatvalues = hatvalues,
    residuals = residuals,
    cooks_distance = cooks_distance,
    vcov = vcov,
    deviance = deviance,
    pseudoR2 = pseudoR2,
    npar = npar,
    w = w,
    y = y, # problems with model.response later if not here
    size = size
  )
}

#' Package the fixed-effect covariance matrix into a fitted model's \code{vcov} element
#'
#' Warns if any diagonal entry is negative (numerically unstable fit).
#'
#' @param cov_betahat The fixed-effect coefficient covariance matrix.
#'
#' @return A list with \code{fixed}.
#'
#' @noRd
get_vcov <- function(cov_betahat) {
  vcov_fixed <- cov_betahat
  if (any(diag(vcov_fixed) < 0)) {
    warning("Model fit potentially unstable. Consider fixing nugget (via nugget_initial) at some non-zero value greater than 1e-4 and refitting the model.", call. = FALSE)
  }
  vcov <- list(fixed = vcov_fixed)
}
