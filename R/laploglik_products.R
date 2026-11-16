#' Get products involved in Laplace log likelihood (glms)
#'
#' @param params_object Covariance parameter object
#' @param data_object Data object
#' @param estmethod Estimation Method
#'
#' @noRd
laploglik_products <- function(params_object, data_object, estmethod) {

  cov_matrix_list <- get_cov_matrix_list(params_object, data_object)

  # cholesky products (no local)
  if (data_object$parallel) {
    cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
      cluster_list_element <- list(
        c = cov_matrix_list[[l]],
        x = data_object$X_list[[l]],
        y = data_object$y_list[[l]]
      )
    })
    cholprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_cholprods_glm_parallel)
    names(cholprods_list) <- names(cov_matrix_list)
  } else {
    cholprods_list <- mapply(
      c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
      function(c, x, y) get_cholprods_glm(c, x, y),
      SIMPLIFY = FALSE
    )
  }

  SigInv_list <- lapply(cholprods_list, function(x) x$SigInv)
  SigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SigInv_X))

  # storing relevant products
  ## lower chol %*% X
  SqrtSigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_X))
  ## lower chol %*% y
  SqrtSigInv_y <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_y))
  # covariance of beta hat
  ## t(X) %*% sigma_inverse %*% X
  Xt_SigInv_X <- crossprod(SqrtSigInv_X, SqrtSigInv_X)
  ## t(X) %*% sigma_inverse %*% X)^(-1)
  Xt_SigInv_X_upchol <- chol(Xt_SigInv_X)
  cov_betahat <- chol2inv(Xt_SigInv_X_upchol)

  # find dispersion
  dispersion <- as.vector(params_object$dispersion) # take class away

  # Newton-Raphson
  w_and_H <- get_w_and_H(
    data_object, dispersion,
    SigInv_list, SigInv_X, cov_betahat, Xt_SigInv_X, estmethod
  )

  w <- w_and_H$w
  # H <- w_and_H$H
  mHldet <- w_and_H$mHldet

  betahat <- tcrossprod(cov_betahat, SigInv_X) %*% w
  X <- do.call("rbind", data_object$X_list)
  r <- w - X %*% betahat
  sizes <- vapply(SigInv_list, nrow, 1L)
  idx <- split(seq_len(sum(sizes)), rep(seq_along(sizes), sizes))
  rt_SigInv_r <- sum(mapply(function(S, i) as.numeric(crossprod(r[i], S %*% r[i])), SigInv_list, idx))

  # get wolfinger objects
  y <- as.vector(do.call("rbind", data_object$y_list))
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }
  l00 <- get_l00(data_object$family, w, y, data_object$size, dispersion)
  l01 <- mHldet
  l1 <- sum(unlist(lapply(cholprods_list, function(x) 2 * sum(log(diag(x$Sig_lowchol))))))
  l2 <- as.numeric(rt_SigInv_r)

  # returning relevant quantities
  if (estmethod == "reml") {
    l3 <- 2 * sum(log(diag(Xt_SigInv_X_upchol)))
    return(list(l00 = l00, l01 = l01, l1 = l1, l2 = l2, l3 = l3))
  }

  if (estmethod == "ml") {
    return(list(l00 = l00, l01 = l01, l1 = l1, l2 = l2))
  }
}

#' Get minus twice the Laplace log likelihood
#'
#' @param laploglik_products Laplace log likelihood products
#' @param data_object Data object
#' @param estmethod Estimation method
#'
#' @noRd
get_minustwolaploglik <- function(laploglik_products, data_object, estmethod) {
# On the -2log-likelihood scale, the Laplace constant cancels the Gaussian n log(2 pi)
  if (estmethod == "reml") {
    minustwolaploglik <- as.numeric(laploglik_products$l00 + laploglik_products$l01 + laploglik_products$l1 +
      laploglik_products$l2 + laploglik_products$l3 +
      -data_object$p * log(2 * pi))
  } else if (estmethod == "ml") {
    minustwolaploglik <- as.numeric(laploglik_products$l00 + laploglik_products$l01 + laploglik_products$l1 +
      laploglik_products$l2)
  }

  minustwolaploglik
}


#' Get prediction and hessian of latent effects
#'
#' @param data_object Data object
#' @param dispersion Dispersion parameter
#' @param SigInv_list List of inverse covariance matrices for each block
#' @param SigInv_X Product between inverse covariance matrix and model matrix
#' @param cov_betahat Covariance of fixed effects
#' @param cov_betahat_Inv Inverse covariance of fixed effects (precision)
#' @param estmethod Estimation method
#' @param ret_mHInv Should -(H)^-1 be retained for later use?
#'
#' @noRd
get_w_and_H <- function(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
  family <- data_object$family
  y <- as.vector(do.call("rbind", data_object$y_list))
  size <- data_object$size
  w <- get_w_init(family, y, dispersion)
  wdiffmax <- Inf
  iter <- 0

  # The offset is a known, non-estimated shift on the link scale, so w and the
  # linear predictor are not the same vector: w is the offset-free latent
  # process the optimizer solves for (and the scale of the
  # spatial covariance), while the family log-likelihood is always evaluated at
  # the linear predictor w + offset. Every get_d()/get_D() call below therefore
  # takes w + offset, and every Ptheta product takes w alone. Defaulting the
  # offset to 0 keeps that distinction visible in one place instead of
  # duplicating it across if/else branches, matching glm().
  offset <- if (is.null(data_object$offset)) 0 else as.vector(data_object$offset)

  # single-partition case: the Hessian is small enough to solve directly
  if (length(SigInv_list) == 1) {
    SigInv <- as.matrix(SigInv_list[[1]])
    SigInv_X <- as.matrix(SigInv_X)
    cov_betahat <- as.matrix(cov_betahat)
    Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
    while (iter < 50 && wdiffmax > 1e-4) {
      iter <- iter + 1
      # compute the d vector
      d <- get_d(family, w + offset, y, size, dispersion)
      # and then the gradient vector
      g <- as.vector(d - Ptheta %*% w)
      # Next, compute H
      D <- get_D(family, w + offset, y, size, dispersion)
      mH <- Ptheta
      diag(mH) <- diag(mH) - as.vector(diag(D))
      R <- tryCatch(chol(mH), error = function(e) NULL)
      solveHg <- if (is.null(R)) solve(-mH, g) else -backsolve(R, forwardsolve(t(R), g))
      solveHg <- as.vector(solveHg)
      wnew <- w - solveHg
      # check overshoot on loglik surface
      dnew <- get_d(family, wnew + offset, y, size, dispersion)
      gnew <- dnew - Ptheta %*% wnew
      if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing nugget at a known, non-zero value (via nugget_initial), or fixing dispersion at one (via dispersion_initial).", call. = FALSE)
      if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
      wdiffmax <- max(abs(wnew - w))
      w <- wnew
    }

    # ML holds beta fixed in the integral; Newton uses the profiled Hessian.
    if (estmethod == "ml") {
      m <- SigInv
      diag(m) <- diag(m) - as.vector(diag(D))
      Rm <- tryCatch(chol(m), error = function(e) NULL)
      mHldet <- if (is.null(Rm)) as.numeric(determinant(m, logarithm = TRUE)$modulus) else 2 * sum(log(diag(Rm)))
    } else {
      mHldet <- if (is.null(R)) as.numeric(determinant(mH, logarithm = TRUE)$modulus) else 2 * sum(log(diag(R)))
    }
    w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
    if (ret_mHInv) {
      w_and_H_list$mHInv <- if (is.null(R)) solve(mH) else chol2inv(R)
    }
  } else {
    SigInv_list <- lapply(SigInv_list, as.matrix)
    SigInv_X <- as.matrix(SigInv_X)
    cov_betahat <- as.matrix(cov_betahat)
    cov_betahat_Inv <- as.matrix(cov_betahat_Inv)
    sizes <- vapply(SigInv_list, nrow, 1L)
    idx <- split(seq_len(sum(sizes)), rep(seq_along(sizes), sizes))
    stopifnot(identical(unname(lengths(idx)),
      unname(lengths(split(seq_along(y), sort(data_object$local_index))))))
    Pw <- function(v) {
      Sv <- unlist(Map(function(S, i) S %*% v[i], SigInv_list, idx), use.names = FALSE)
      Sv - as.vector(SigInv_X %*% (cov_betahat %*% crossprod(SigInv_X, v)))
    }
    p <- ncol(SigInv_X)
    # add cov_betahat_Inv stability by same diagonal tolerance as this can have problems too
    diag(cov_betahat_Inv) <- diag(cov_betahat_Inv) + data_object$diagtol

    while (iter < 50 && wdiffmax > 1e-4) {
      iter <- iter + 1
      # compute the d vector
      d <- get_d(family, w + offset, y, size, dispersion)
      # and then the gradient vector
      g <- d - Pw(w)
      # Next, compute H
      D <- get_D(family, w + offset, y, size, dispersion)
      F <- factor_laplace_blocks(SigInv_list, as.vector(diag(D)), idx,
        data_object$parallel, data_object$cl)
      AUg <- apply_Ainv(F, cbind(SigInv_X, g), idx)
      AU <- AUg[, seq_len(p), drop = FALSE]
      Ag <- AUg[, p + 1L]
      mid <- cov_betahat_Inv + crossprod(SigInv_X, AU)
      # Woodbury product; see spmodel's technical.Rmd #sec:bigdata-glm.
      solveHg <- as.vector(Ag - AU %*% solve(mid, crossprod(SigInv_X, Ag)))
      wnew <- w - solveHg
      # check overshoot on loglik surface
      dnew <- get_d(family, wnew + offset, y, size, dispersion)
      gnew <- dnew - Pw(wnew)
      if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing nugget at a known, non-zero value (via nugget_initial), or fixing dispersion at one (via dispersion_initial).", call. = FALSE)
      if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
      wdiffmax <- max(abs(wnew - w))
      # update w
      w <- wnew
    }

    mHldet <- logdet_A(F)
    if (estmethod == "reml") {
      mHldet <- mHldet + 2 * sum(log(diag(chol(cov_betahat)))) +
        as.numeric(determinant(mid, logarithm = TRUE)$modulus)
    }
    w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
    if (ret_mHInv) {
      w_and_H_list$mHInv <- structure(list(F = F, idx = idx, AU = AU, mid_inv = solve(mid)),
        class = "laplace_mhinv")
    }
  }

  w_and_H_list
}

factor_laplace_block <- function(block) {
  m <- block$S
  diag(m) <- diag(m) - block$D
  R <- tryCatch(chol(m), error = function(e) NULL)
  if (is.null(R)) list(type = "lu", A = -m) else list(type = "chol", R = R)
}

factor_laplace_blocks <- function(S_list, D, idx, parallel, cl) {
  blocks <- Map(function(S, i) list(S = S, D = D[i]), S_list, idx)
  if (parallel) parallel::parLapply(cl, blocks, factor_laplace_block) else lapply(blocks, factor_laplace_block)
}

apply_Ainv <- function(F, V, idx) {
  V <- as.matrix(V)
  do.call(rbind, Map(function(f, i) {
    if (f$type == "chol") {
      -backsolve(f$R, forwardsolve(t(f$R), V[i, , drop = FALSE]))
    } else {
      solve(f$A, V[i, , drop = FALSE])
    }
  }, F, idx))
}

logdet_A <- function(F) {
  sum(vapply(F, function(f) {
    if (f$type == "chol") 2 * sum(log(diag(f$R))) else as.numeric(determinant(f$A, logarithm = TRUE)$modulus)
  }, numeric(1)))
}

mhinv_apply <- function(M, V) {
  if (!inherits(M, "laplace_mhinv")) return(as.matrix(M %*% V))
  V <- as.matrix(V)
  -(apply_Ainv(M$F, V, M$idx) - M$AU %*% (M$mid_inv %*% crossprod(M$AU, V)))
}

mhinv_diag <- function(M) {
  if (!inherits(M, "laplace_mhinv")) return(Matrix::diag(M))
  diag_Ainv <- unlist(lapply(M$F, function(f) {
    if (f$type == "chol") -diag(chol2inv(f$R)) else diag(solve(f$A))
  }), use.names = FALSE)
  -(diag_Ainv - rowSums((M$AU %*% M$mid_inv) * M$AU))
}

#' Get the gradient of w (called d)
#'
#' @param family The glm family
#' @param w The latent effect (link mean)
#' @param y The response
#' @param size Number of trails (for binomial response)
#' @param dispersion Dispersion parameter
#'
#' @noRd
get_d <- function(family, w, y, size, dispersion) {
  if (family == "poisson") {
    d <- -exp(w) + y
  } else if (family == "nbinomial") {
    d <- dispersion * (y - exp(w)) / (dispersion + exp(w))
  } else if (family == "binomial") {
    d <- y - size * expit(w)
  } else if (family == "Gamma") {
    d <- -dispersion + dispersion * y * exp(-w)
  } else if (family == "inverse.gaussian") {
    # d <- 1 / dispersion * (y - exp(w)) / exp(2 * w)
    d <- dispersion * (y / (2 * exp(w)) - exp(w) / (2 * y)) + 1 / 2
  } else if (family == "beta") {
    one_expw <- 1 + exp(w)
    k0 <- digamma(dispersion * exp(w) / one_expw) - digamma(dispersion / one_expw) + log(1 / y - 1)
    d <- -dispersion * exp(w) * k0 / one_expw^2
  }
  d
}

#' Get the Hessian of w (called D)
#'
#' @param family The glm family
#' @param w The latent effect (link mean)
#' @param y The response
#' @param size Number of trails (for binomial response)
#' @param dispersion Dispersion parameter
#'
#' @noRd
get_D <- function(family, w, y, size, dispersion) {
  w <- as.vector(w)

  if (family == "poisson") {
    D_vec <- -exp(w)
  } else if (family == "nbinomial") {
    D_vec <- -(dispersion * exp(w) * (dispersion + y)) / ((dispersion + exp(w))^2)
  } else if (family == "binomial") {
    D_vec <- -size * expit(w) / (1 + exp(w))
  } else if (family == "Gamma") {
    D_vec <- -dispersion * y * exp(-w)
  } else if (family == "inverse.gaussian") {
    # D_vec <- 1 / dispersion * (exp(w) - 2 * y) / exp(2 * w)
    D_vec <- -dispersion * (exp(2 * w) + y^2) / (2 * y * exp(w))
  } else if (family == "beta") {
    one_expw <- 1 + exp(w)
    k0 <- digamma(dispersion * exp(w) / one_expw) - digamma(dispersion / one_expw) + log(1 / y - 1)
    k1 <- dispersion * (trigamma(dispersion * exp(w) / one_expw) + trigamma(dispersion / one_expw)) - 2 * sinh(w) * k0
    D_vec <- -dispersion * exp(2 * w) * k1 / one_expw^4
  }
  D <- Diagonal(x = D_vec)
}

get_D3 <- function(family, w, y, size, dispersion) {
  w <- as.numeric(w)
  mu <- exp(w)
  if (family == "poisson") {
    value <- -mu
  } else if (family == "nbinomial") {
    p <- stats::plogis(w - log(dispersion))
    value <- -(dispersion + y) * p * (1 - p) * (1 - 2 * p)
  } else if (family == "binomial") {
    p <- expit(w)
    value <- -size * p * (1 - p) * (1 - 2 * p)
  } else if (family == "Gamma") {
    value <- dispersion * y * exp(-w)
  } else if (family == "inverse.gaussian") {
    value <- -dispersion * (mu - y^2 * exp(-w)) / (2 * y)
  } else if (family == "beta") {
    p <- expit(w)
    a <- p * (1 - p)
    b <- 1 - 2 * p
    k0 <- digamma(dispersion * p) - digamma(dispersion * (1 - p)) + log((1 - y) / y)
    k1 <- trigamma(dispersion * p) + trigamma(dispersion * (1 - p))
    k2 <- psigamma(dispersion * p, 2) - psigamma(dispersion * (1 - p), 2)
    value <- -dispersion * a * (b^2 - 2 * a) * k0 -
      3 * dispersion^2 * a^2 * b * k1 - dispersion^3 * a^3 * k2
  }
  as.numeric(value)
}

get_laplace_correction <- function(mHInv, family, eta, y, size, dispersion) {
  v <- mhinv_diag(mHInv)
  t3 <- get_D3(family, eta, y, size, dispersion)
  delta <- 0.5 * as.numeric(mhinv_apply(mHInv, t3 * v))
  list(delta = delta)
}

get_latent_center <- function(object, fit_correct = TRUE) {
  if (!is.logical(fit_correct) || length(fit_correct) != 1L || is.na(fit_correct)) {
    stop("fit_correct must be TRUE or FALSE.", call. = FALSE)
  }
  eta <- fitted(object, type = "link")
  active <- fit_correct && !is.null(object$laplace)
  if (!active) return(list(eta = eta, betahat = coef(object), active = FALSE))
  list(eta = eta + object$laplace$delta[names(eta)], betahat = object$laplace$betahat, active = TRUE)
}

#' Get initial values for w
#'
#' @param w The latent effect (link mean)
#' @param y The response
#' @param dispersion Dispersion parameter
#'
#' @noRd
get_w_init <- function(family, y, dispersion) {
  if (family == "poisson") {
    w_init <- 0.5 * log(y + 1)
  } else if (family == "nbinomial") {
    w_init <- 0.5 * log(y + 1)
  } else if (family == "binomial") {
    w_init <- rep(0, times = length(y))
  } else if (family == "Gamma") {
    w_init <- 0.5 * log(y + 1)
  } else if (family == "inverse.gaussian") {
    w_init <- 0.5 * log(y + 1)
  } else if (family == "beta") {
    w_init <- rep(0, times = length(y))
  }
  w_init
}

#' Get glm distribution piece of Laplace log likelihood
#'
#' @param family The glm family
#' @param w The latent effect (link mean)
#' @param y The response
#' @param size Number of trails (for binomial response)
#' @param dispersion Dispersion parameter
#'
#' @noRd
get_l00 <- function(family, w, y, size, dispersion) {
  w <- as.vector(w)
  y <- as.vector(y)
  # -2 is for -2ll constant
  if (family == "poisson") {
    mu <- exp(w)
    l00 <- -2 * sum(dpois(y, lambda = mu, log = TRUE))
  } else if (family == "nbinomial") {
    mu <- exp(w)
    l00 <- -2 * sum(dnbinom(x = y, mu = mu, size = dispersion, log = TRUE))
  } else if (family == "binomial") {
    mu <- expit(w)
    l00 <- -2 * sum(dbinom(y, size, mu, log = TRUE))
  } else if (family == "Gamma") {
    mu <- exp(w)
    # disp_recip <- 1 / dispersion
    # l00 <- -2 * sum(dgamma(y, shape = disp_recip, scale = dispersion * mu, log = TRUE))
    l00 <- -2 * sum(dgamma(y, shape = dispersion, scale = mu / dispersion, log = TRUE))
  } else if (family == "inverse.gaussian") {
    mu <- exp(w)
    # disp_recip <- 1 / dispersion
    # l00 <- -2 * sum((log(disp_recip) - log(2 * pi) - 3 * log(y)) / 2 - (disp_recip * (y - mu)^2 / (2 * y * mu^2)))
    l00 <- -2 * sum(1 / 2 * (log(dispersion) + log(exp(w)) - log(2 * pi) - log(y^3)) - dispersion * (y - exp(w))^2 / (2 * exp(w) * y))
  } else if (family == "beta") {
    mu <- expit(w)
    a <- mu * dispersion
    b <- (1 - mu) * dispersion
    l00 <- -2 * sum(dbeta(x = y, shape1 = a, shape2 = b, log = TRUE))
  }
  l00
}
