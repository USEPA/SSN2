#' Get relevant model fit statistics and diagnostics
#'
#' @param cov_est_object Covariance parameter estimation object.
#' @param data_object Data object.
#' @param estmethod Estimation method.
#'
#' @noRd
get_model_stats <- function(cov_est_object, data_object, estmethod) {
  # store the covariance matrix list
  cov_matrix_list <- get_cov_matrix_list(cov_est_object$params_object, data_object)

  # compute eigenproducts for later use
  if (data_object$parallel) { # put back when local implemented
    # cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
    #   cluster_list_element <- list(
    #     c = cov_matrix_list[[l]],
    #     x = data_object$X_list[[l]],
    #     y = data_object$y_list[[l]],
    #     o = data_object$ones_list[[l]]
    #   )
    # })
    # eigenprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_eigenprods_parallel)
    # names(eigenprods_list) <- names(cov_matrix_list)
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

  # adjust the global covariance of beta hat (revisit when local implemented)
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
  fitted$response <- fitted$response[order(data_object$order)]
  names(fitted$response) <- model_stats_names
  fitted$tailup <- fitted$tailup[order(data_object$order)]
  names(fitted$tailup) <- model_stats_names
  fitted$taildown <- fitted$taildown[order(data_object$order)]
  names(fitted$taildown) <- model_stats_names
  fitted$euclid <- fitted$euclid[order(data_object$order)]
  names(fitted$euclid) <- model_stats_names
  fitted$nugget <- fitted$nugget[order(data_object$order)]
  names(fitted$nugget) <- model_stats_names
  ## hat values
  hatvalues <- hatvalues[order(data_object$order)]
  names(hatvalues) <- model_stats_names
  ## residuals
  residuals$response <- residuals$response[order(data_object$order)]
  names(residuals$response) <- model_stats_names
  residuals$pearson <- residuals$pearson[order(data_object$order)]
  names(residuals$pearson) <- model_stats_names
  residuals$standardized <- residuals$standardized[order(data_object$order)]
  names(residuals$standardized) <- model_stats_names
  ## cook's distance
  cooks_distance <- cooks_distance[order(data_object$order)]
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

###############################################################################
############### helpers to store relevant statistics as list objects
###############################################################################

get_coefficients <- function(betahat, params_object) {
  # store fixed effect and random coefficients as list
  list(fixed = betahat, params_object = params_object)
}

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
    tailup_list <- lapply(data_object$dist_object_oblist, function(x) cov_matrix(params_object$tailup, x))
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
    taildown_list <- lapply(data_object$dist_object_oblist, function(x) cov_matrix(params_object$taildown, x))
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
    euclid_list <- lapply(data_object$dist_object_oblist, function(x) cov_matrix(params_object$euclid, x, data_object$anisotropy))
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

get_hatvalues <- function(cov_betahat, SqrtSigInv_X) {
  # the hat matrix of the whitened residuals
  hatvalues <- diag(SqrtSigInv_X %*% tcrossprod(cov_betahat, SqrtSigInv_X))
  as.numeric(hatvalues)
}

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

get_cooks_distance <- function(residuals, hatvalues, p) {
  # find cook's distance
  residuals$pearson^2 * hatvalues / (p * (1 - hatvalues))
}

get_vcov <- function(cov_betahat) {
  # only operational for the fixed effects right now
  vcov_fixed <- cov_betahat
  if (any(diag(vcov_fixed) < 0)) {
    warning("Model fit potentially unstable. Consider fixing nugget (via nugget_initial) at some non-zero value greater than 1e-4 and refitting the model.", call. = FALSE)
  }
  vcov <- list(fixed = vcov_fixed)
}
