#' Get relevant model fit statistics and diagnostics for glms
#'
#' @param cov_est_object Covariance parameter estimation object.
#' @param data_object Data object.
#' @param estmethod Estimation method.
#'
#' @noRd
get_model_stats_glm <- function(cov_est_object, data_object, estmethod) {
  cov_matrix_list <- get_cov_matrix_list(cov_est_object$params_object, data_object)

  # find model components
  X <- do.call("rbind", data_object$X_list)
  y <- do.call("rbind", data_object$y_list)

  # eigen products (put back when local impelmented)
  # if (data_object$parallel) {
  #   cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
  #     cluster_list_element <- list(
  #       c = cov_matrix_list[[l]],
  #       x = data_object$X_list[[l]],
  #       y = data_object$y_list[[l]],
  #       o = data_object$ones_list[[l]]
  #     )
  #   })
  #   eigenprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_eigenprods_glm_parallel)
  #   names(eigenprods_list) <- names(cov_matrix_list)
  # } else {
  #   eigenprods_list <- mapply(
  #     c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, o = data_object$ones_list,
  #     function(c, x, y, o) get_eigenprods_glm(c, x, y, o),
  #     SIMPLIFY = FALSE
  #   )
  # }

  eigenprods_list <- mapply(
    c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, o = data_object$ones_list,
    function(c, x, y, o) get_eigenprods_glm(c, x, y, o),
    SIMPLIFY = FALSE
  )

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

  # newton rhapson
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

  # adjust the global covariance of beta hat (revisit when local implemented)
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

  # newton rhapson
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
  fitted$response <- fitted$response[order(data_object$order)]
  names(fitted$response) <- model_stats_glm_names
  fitted$link <- fitted$link[order(data_object$order)]
  names(fitted$link) <- model_stats_glm_names
  fitted$tailup <- fitted$tailup[order(data_object$order)]
  names(fitted$tailup) <- model_stats_glm_names
  fitted$taildown <- fitted$taildown[order(data_object$order)]
  names(fitted$taildown) <- model_stats_glm_names
  fitted$euclid <- fitted$euclid[order(data_object$order)]
  names(fitted$euclid) <- model_stats_glm_names
  fitted$nugget <- fitted$nugget[order(data_object$order)]
  names(fitted$nugget) <- model_stats_glm_names
  ## hat values
  hatvalues <- hatvalues[order(data_object$order)]
  names(hatvalues) <- model_stats_glm_names
  ## residuals
  residuals$response <- residuals$response[order(data_object$order)]
  names(residuals$response) <- model_stats_glm_names
  residuals$deviance <- residuals$deviance[order(data_object$order)]
  names(residuals$deviance) <- model_stats_glm_names
  residuals$pearson <- residuals$pearson[order(data_object$order)]
  names(residuals$pearson) <- model_stats_glm_names
  residuals$standardized <- residuals$standardized[order(data_object$order)]
  names(residuals$standardized) <- model_stats_glm_names
  ## cook's distance
  cooks_distance <- cooks_distance[order(data_object$order)]
  names(cooks_distance) <- model_stats_glm_names
  y <- y[order(data_object$order)]
  if (is.null(data_object$size)) {
    size <- NULL
  } else {
    size <- data_object$size[order(data_object$order)]
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

###############################################################################
############### helpers to store relevant statistics as list objects
###############################################################################

get_coefficients_glm <- function(betahat, params_object) {
  list(fixed = betahat, params_object = params_object)
}

get_fitted_glm <- function(w_list, betahat, params_object, data_object, eigenprods_list) {
  fitted_link <- unname(do.call("c", w_list)) # unlist(w_list, use.names = FALSE)
  # add offset
  if (!is.null(data_object$offset)) {
    fitted_link <- fitted_link + as.vector(data_object$offset)
  }
  fitted_response <- invlink(fitted_link, data_object$family, data_object$size)

  SigInv_r_list <- mapply(
    x = eigenprods_list, w = w_list, function(x, w) x$SigInv %*% as.matrix(w, ncol = 1) - x$SigInv_X %*% betahat,
    SIMPLIFY = FALSE
  )

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
    link = as.numeric(fitted_link),
    tailup = as.numeric(fitted_tailup),
    taildown = as.numeric(fitted_taildown),
    euclid = as.numeric(fitted_euclid),
    nugget = as.numeric(fitted_nugget),
    randcov = fitted_randcov
  )
}

get_fitted_null <- function(w, data_object) {
  fitted_link <- as.numeric(w)
  # add offset
  if (!is.null(data_object$offset)) {
    fitted_link <- fitted_link + data_object$offset
  }
  # fitted_link
  fitted_response <- invlink(fitted_link, data_object$family, data_object$size)
}

invlink <- function(fitted_link, family, size) {
  if (family == "poisson") {
    fitted <- exp(fitted_link)
  } else if (family == "binomial") {
    if (is.null(size)) size <- 1
    fitted <- size * expit(fitted_link)
  } else if (family == "nbinomial") {
    fitted <- exp(fitted_link)
  } else if (family == "Gamma") {
    # fitted <- 1 / fitted_link
    fitted <- exp(fitted_link)
  } else if (family == "inverse.gaussian") {
    fitted <- exp(fitted_link)
  } else if (family == "beta") {
    fitted <- expit(fitted_link)
  }
  fitted
}

get_hatvalues_glm <- function(w, X, data_object, dispersion) {
  # the hat matrix of the whitened residuals
  V <- get_V(w, data_object$family, data_object$size, dispersion)
  SqrtVInv_X <- sqrt(V) * X # same as diag(sqrt(V)) %*% X
  cov_vhat <- chol2inv(chol(Matrix::forceSymmetric(crossprod(SqrtVInv_X, SqrtVInv_X))))
  hatvalues <- diag(SqrtVInv_X %*% tcrossprod(cov_vhat, SqrtVInv_X))
  if (any(hatvalues > 0.999)) {
    hatvalues_sum <- sum(hatvalues)
    hatvalues[hatvalues > 0.999] <- 0.999
    hatvalues <- hatvalues * (hatvalues_sum / sum(hatvalues))
  }
  as.numeric(hatvalues)
}

get_V <- function(w, family, size, dispersion) {
  # V = (1 / dispersion) * (dmu / deta)^2 * (V(mu))
  # when canonical link used, dmu / deta = dmu / dtheta = V(mu)
  # so V = (1 / dispersion) * V(mu)
  # V is such that cov(betahat) = (XtVX/dispersion)^{-1}
  # and hence V^{-1}dispersion = Sigma used in fitting
  # but V equals var(y) when dispersion is one
  if (family == "poisson") {
    mu <- exp(w)
    V <- mu
  } else if (family == "binomial") {
    mu <- expit(w)
    V <- size * mu * (1 - mu)
  } else if (family == "nbinomial") {
    mu <- exp(w)
    V <- mu / (1 + (mu / dispersion)) # from Ver Hoef and Boveng 2007
  } else if (family == "Gamma") {
    mu <- exp(w)
    V <- mu^2
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    V <- mu^3
  } else if (family == "beta") {
    mu <- expit(w)
    V <- mu * (1 - mu)
  }
  V
}

get_var_y <- function(w, family, size, dispersion) {
  # var(y) = dispersion * var(mu)
  # when dispersion = 1, var(y) = var(mu)
  if (family == "poisson") {
    var_y <- get_V(w, family, size, dispersion)
  } else if (family == "binomial") {
    var_y <- get_V(w, family, size, dispersion)
  } else if (family == "nbinomial") {
    mu <- exp(w)
    var_y <- mu + mu^2 / dispersion
  } else if (family == "Gamma") {
    dispersion_true <- 1 / dispersion
    var_y <- get_V(w, family, size, dispersion) * dispersion_true
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    dispersion_true <- 1 / (mu * dispersion)
    var_y <- get_V(w, family, size, dispersion) * dispersion_true
  } else if (family == "beta") {
    dispersion_true <- 1 / (1 + dispersion)
    var_y <- get_V(w, family, size, dispersion) * dispersion_true
  }
  var_y
}

get_deviance_glm <- function(family, y, fitted_response, size, dispersion) {
  # if (!is.null(offset)) {
  #   fitted_link <- fitted_link + offset # undo w = w - offset for deviance to match glm
  # }

  # fitted_response <- invlink(fitted_link, family, size)

  # faraway p 157
  # y <- pmax(y, 1e-8) # so deviance  Inf is not calculated
  if (family == "poisson") {
    half_deviance_i <- ifelse(y == 0, 0, y * log(y / fitted_response)) - (y - fitted_response)
    # half_deviance_i <- y * pmax(-1e10, log(y / fitted_response)) - (y - fitted_response)
  } else if (family == "binomial") {
    half_deviance_i <- ifelse(y == 0, 0, y * log(y / fitted_response)) +
      ifelse(size - y == 0, 0, (size - y) * log((size - y) / (size - fitted_response)))
  } else if (family == "nbinomial") {
    # hand derived
    half_deviance_i <- ifelse(y == 0, 0, y * (log(y / (y + dispersion)) - log(fitted_response / (fitted_response + dispersion)))) +
      dispersion * (log(fitted_response + dispersion) - log(y + dispersion))
  } else if (family == "Gamma") {
    half_deviance_i <- -log(y / fitted_response) + (y - fitted_response) / fitted_response
  } else if (family == "inverse.gaussian") {
    half_deviance_i <- 0.5 * (y - fitted_response)^2 / (y * fitted_response^2)
  } else if (family == "beta") {
    constant <- log(gamma(fitted_response * dispersion)) + log(gamma((1 - fitted_response) * dispersion)) - log(gamma(y * dispersion)) - log(gamma((1 - y) * dispersion))
    half_deviance_i <- constant + (y - fitted_response) * dispersion * log(y) + ((1 - y) - (1 - fitted_response)) * dispersion * log(1 - y)
  }
  deviance_i <- 2 * half_deviance_i
  as.numeric(deviance_i)
}

get_residuals_glm <- function(w, y, data_object, deviance_i, hatvalues, dispersion) {
  # add offset
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }

  residuals_response <- y - invlink(w, data_object$family, data_object$size)

  residuals_deviance <- sign(residuals_response) * sqrt(deviance_i)

  residuals_pearson <- residuals_response / sqrt(get_var_y(w, data_object$family, data_object$size, dispersion))

  residuals_standardized <- residuals_deviance / sqrt(1 - hatvalues) # (I - H on bottom)
  list(
    response = as.numeric(residuals_response), deviance = as.numeric(residuals_deviance),
    pearson = as.numeric(residuals_pearson), standardized = as.numeric(residuals_standardized)
  )
}

get_cooks_distance_glm <- function(residuals, hatvalues, p) {
  residuals$standardized^2 * hatvalues / (p * (1 - hatvalues))
}

get_vcov_glm <- function(cov_betahat_corrected, cov_betahat_uncorrected) {
  if (any(diag(cov_betahat_corrected) < 0)) {
    warning("Model fit potentially unstable. Consider fixing ie (via spcov_initial) at some non-zero value greater than 1e-4 and refitting the model.", call. = FALSE)
  }
  vcov_fixed <- list(corrected = cov_betahat_corrected, uncorrected = cov_betahat_uncorrected)
  vcov <- list(fixed = vcov_fixed)
}
