get_model_stats_glm <- function(cov_est_object, data_object, estmethod) {
  cov_matrix_list <- get_cov_matrix_list(cov_est_object$params_object, data_object)

  eigenprods_list <- mapply(
    c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, o = data_object$ones_list,
    function(c, x, y, o) get_eigenprods_glm(c, x, y, o),
    SIMPLIFY = FALSE
  )

  get_model_stats_glm_core(cov_est_object, data_object, estmethod, data_object$order, eigenprods_list)
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

get_delta_se <- function(fit, se.fit, family, newdata_size = 1) {
  if (family %in% c("poisson", "nbinomial", "Gamma", "inverse.gaussian")) {
    g <- exp(fit)
  } else if (family %in% c("binomial", "beta")) {
    fit <- expit(fit)
    g <- fit * (1 - fit)
  }
  val <- se.fit * g # the square root of the delta method variance

  # binomial response is a proportion on the link scale; rescale back up to
  # the count scale by multiplying by the number of trials, mirroring
  # invlink()'s own size * expit(fitted_link)
  if (family == "binomial") {
    if (is.null(newdata_size)) newdata_size <- 1
    val <- val * newdata_size
  }

  val
}

get_hatvalues_glm <- function(w, X, data_object, dispersion) {
  if (!is.null(data_object$offset)) {
    w <- as.vector(w + data_object$offset)
  }

  # the hat matrix of the whitened residuals
  V <- get_V(w, data_object$family, data_object$size, dispersion)
  SqrtVInv_X <- sqrt(V) * X # same as diag(sqrt(V)) %*% X
  cov_vhat <- chol2inv(chol(Matrix::forceSymmetric(crossprod(SqrtVInv_X, SqrtVInv_X))))
  # only the diagonal of the whitened hat matrix is needed, so skip forming
  # the full n x n product
  hatvalues <- get_diag_XVXt(SqrtVInv_X, cov_vhat)
  if (any(hatvalues > 0.999)) {
    hatvalues_sum <- sum(hatvalues)
    hatvalues[hatvalues > 0.999] <- 0.999
    hatvalues <- hatvalues * (hatvalues_sum / sum(hatvalues))
  }
  as.numeric(hatvalues)
}

get_V <- function(w, family, size, dispersion) {
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
    # E(-D_ii) = dispersion, free of mu (this is glm(family = Gamma(link = "log"))'s hat matrix)
    V <- rep(dispersion, length(w))
  } else if (family == "inverse.gaussian") {
    # E(-D_ii) = dispersion + 1/2, free of mu
    V <- rep(dispersion + 0.5, length(w))
  } else if (family == "beta") {
    mu <- expit(w)
    V <- dispersion^2 * (mu * (1 - mu))^2 *
      (trigamma(mu * dispersion) + trigamma((1 - mu) * dispersion))
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
    # mu^2 / dispersion; computed directly (not via get_V()) since get_V()'s
    # Gamma branch is now the leverage weight (constant in mu), not V(mu)
    mu <- exp(w)
    var_y <- mu^2 / dispersion
  } else if (family == "inverse.gaussian") {
    # mu^2 / dispersion; computed directly for the same reason as Gamma above
    mu <- exp(w)
    var_y <- mu^2 / dispersion
  } else if (family == "beta") {
    mu <- expit(w)
    var_y <- mu * (1 - mu) / (1 + dispersion)
  }
  var_y
}

get_dispersion_factor <- function(w, family, size, dispersion) {
  n <- length(w)
  if (family %in% c("poisson", "binomial", "nbinomial", "beta")) {
    a_phi <- rep(1, n)
  } else if (family == "Gamma") {
    a_phi <- rep(1 / dispersion, n)
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    a_phi <- 1 / (mu * dispersion)
  }
  a_phi
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
    # has NA problem for large dispersion
    # constant <- log(gamma(fitted_response * dispersion)) + log(gamma((1 - fitted_response) * dispersion)) - log(gamma(y * dispersion)) - log(gamma((1 - y) * dispersion))
    constant <- lgamma(fitted_response * dispersion) + lgamma((1 - fitted_response) * dispersion) - lgamma(y * dispersion) - lgamma((1 - y) * dispersion)
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

  a_phi <- get_dispersion_factor(w, data_object$family, data_object$size, dispersion)
  residuals_standardized <- residuals_deviance / sqrt(a_phi * (1 - hatvalues)) # (I - H on bottom)
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
    warning("Model fit potentially unstable. Consider fixing nugget (via nugget_initial) at some non-zero value greater than 1e-4 and refitting the model.", call. = FALSE)
  }
  vcov_fixed <- list(corrected = cov_betahat_corrected, uncorrected = cov_betahat_uncorrected)
  vcov <- list(fixed = vcov_fixed)
}
