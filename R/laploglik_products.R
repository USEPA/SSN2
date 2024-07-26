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
  # if (data_object$parallel) {
  #   cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
  #     cluster_list_element <- list(
  #       c = cov_matrix_list[[l]],
  #       x = data_object$X_list[[l]],
  #       y = data_object$y_list[[l]]
  #     )
  #   })
  #   cholprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_cholprods_glm_parallel)
  #   names(cholprods_list) <- names(cov_matrix_list)
  # } else {
  #   cholprods_list <- mapply(
  #     c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
  #     function(c, x, y) get_cholprods_glm(c, x, y),
  #     SIMPLIFY = FALSE
  #   )
  # }

  cholprods_list <- mapply(
    c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list,
    function(c, x, y) get_cholprods_glm(c, x, y),
    SIMPLIFY = FALSE
  )

  SigInv_list <- lapply(cholprods_list, function(x) x$SigInv)
  SigInv <- Matrix::bdiag(SigInv_list)
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

  # newton rhapson
  w_and_H <- get_w_and_H(
    data_object, dispersion,
    SigInv_list, SigInv_X, cov_betahat, Xt_SigInv_X, estmethod
  )

  w <- w_and_H$w
  # H <- w_and_H$H
  mHldet <- w_and_H$mHldet

  betahat <- tcrossprod(cov_betahat, SigInv_X) %*% w
  # reset w after finding betahat
  if (!is.null(data_object$offset)) {
    w <- w + data_object$offset
  }
  X <- do.call("rbind", data_object$X_list)
  r <- w - X %*% betahat
  rt_SigInv_r <- crossprod(r, SigInv) %*% r

  # get wolfinger objects
  y <- as.vector(do.call("rbind", data_object$y_list))
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
  if (estmethod == "reml") {
    minustwolaploglik <- as.numeric(laploglik_products$l00 + laploglik_products$l01 + laploglik_products$l1 +
      laploglik_products$l2 + laploglik_products$l3 +
      (data_object$n - data_object$p) * log(2 * pi))
  } else if (estmethod == "ml") {
    minustwolaploglik <- as.numeric(laploglik_products$l00 + laploglik_products$l01 + laploglik_products$l1 +
      laploglik_products$l2 +
      data_object$n * log(2 * pi))
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
  SigInv <- Matrix::bdiag(SigInv_list)
  Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
  y <- as.vector(do.call("rbind", data_object$y_list))
  size <- data_object$size
  w <- get_w_init(family, y, dispersion)
  wdiffmax <- Inf
  iter <- 0



  if (length(SigInv_list) == 1) {
    while (iter < 50 && wdiffmax > 1e-4) {
      iter <- iter + 1
      # if (family %in% c("binomial", "beta")) {
      #   w <- pmax(pmin(w, 8), -8)
      # }
      # compute the d vector
      d <- get_d(family, w, y, size, dispersion)
      # and then the gradient vector
      g <- d - Ptheta %*% w
      # Next, compute H
      D <- get_D(family, w, y, size, dispersion)
      H <- D - Ptheta # not PD but -H is
      solveHg <- solve(H, g)
      wnew <- w - solveHg
      # mH_upchol <- chol(Matrix::forceSymmetric(-H))
      # solveHg <- backsolve(mH_upchol, forwardsolve(t(mH_upchol), g))
      # wnew <- w + solveHg # + because -H is already applied
      # check overshoot on loglik surface
      dnew <- get_d(family, wnew, y, size, dispersion)
      gnew <- dnew - Ptheta %*% wnew
      if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing ie at a known, non-zero value (via spcov_initial), or fixing dispersion at one (via dispersion_initial).", call. = FALSE)
      if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg # + because -H is already applied
      # if (max(abs(gnew)) > max(abs(g))) wnew <- w + 0.1 * solveHg
      wdiffmax <- max(abs(wnew - w))
      # update w
      w <- wnew
    }

    # if (family %in% c("binomial", "beta")) {
    #   w <- pmax(pmin(w, 8), -8)
    # }

    mHldet <- as.numeric(determinant(-H, logarithm = TRUE)$modulus)
    # mHldet <- 2 * sum(log(diag(mH_upchol)))
    w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
    if (ret_mHInv) {
      # not done above because this is only for model stats and solve(H) slower than solve(H, g)
      HInv <- solve(H)
      w_and_H_list$mHInv <- -HInv
      # mHInv <- chol2inv(mH_upchol)
      # w_and_H_list$mHInv <- mHInv
    }
  } else {
    # # add cov_betahat_Inv stability by same diagonal tolerance as this can have problems too
    # diag(cov_betahat_Inv) <- diag(cov_betahat_Inv) + data_object$diagtol
    #
    # while (iter < 50 && wdiffmax > 1e-4) {
    #   iter <- iter + 1
    #   # compute the d vector
    #   d <- get_d(family, w, y, size, dispersion)
    #   # and then the gradient vector
    #   g <- d - Ptheta %*% w
    #   # Next, compute H
    #   D <- get_D(family, w, y, size, dispersion)
    #   D_diag <- diag(D)
    #   D_list <- lapply(split(D_diag, sort(data_object$local_index)), function(x) Diagonal(x = x))
    #   # cholesky products (while local not implemented)
    #   # if (data_object$parallel) {
    #   #   cluster_list <- lapply(seq_along(D_list), function(l) {
    #   #     cluster_list_element <- list(
    #   #       D = D_list[[l]],
    #   #       S = SigInv_list[[l]]
    #   #     )
    #   #   })
    #   #   DSigInv_list <- parallel::parLapply(data_object$cl, cluster_list, get_DSigInv_parallel)
    #   #   names(DSigInv_list) <- names(D_list)
    #   # } else {
    #   #   DSigInv_list <- mapply(
    #   #     D = D_list, S = SigInv_list,
    #   #     function(D, S) get_DSigInv(D, S),
    #   #     SIMPLIFY = FALSE
    #   #   )
    #   # }
    #
    #   DSigInv_list <- mapply(
    #     D = D_list, S = SigInv_list,
    #     function(D, S) get_DSigInv(D, S),
    #     SIMPLIFY = FALSE
    #   )
    #
    #   # while local not impelmented
    #   # if (data_object$parallel) {
    #   #   cluster_list <- DSigInv_list
    #   #   DSigInv_Inv_list <- parallel::parLapply(data_object$cl, cluster_list, solve)
    #   #   names(DSigInv_Inv_list) <- names(D_list)
    #   # } else {
    #   #   DSigInv_Inv_list <- lapply(DSigInv_list, function(x) solve(x))
    #   # }
    #
    #   DSigInv_Inv_list <- lapply(DSigInv_list, function(x) solve(x))
    #
    #   DSigInv_Inv <- Matrix::bdiag(DSigInv_Inv_list)
    #   HInv <- smw_HInv(AInv = DSigInv_Inv, U = SigInv_X, CInv = cov_betahat_Inv)
    #   solveHg <- HInv %*% g
    #   wnew <- w - solveHg
    #   # check overshoot on loglik surface
    #   dnew <- get_d(family, wnew, y, size, dispersion)
    #   gnew <- dnew - Ptheta %*% wnew
    #   if (any(is.na(gnew) | is.infinite(gnew))) stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing ie at a known, non-zero value (via spcov_initial), or fixing dispersion at one (via dispersion_initial).", call. = FALSE)
    #   if (max(abs(gnew)) > max(abs(g))) wnew <- w - 0.1 * solveHg
    #   wdiffmax <- max(abs(wnew - w))
    #   # update w
    #   w <- wnew
    # }
    #
    # mHldet <- smw_mHldet(A_list = DSigInv_list, AInv = DSigInv_Inv, U = SigInv_X, C = cov_betahat, CInv = cov_betahat_Inv)
    # w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
    # if (ret_mHInv) {
    #   w_and_H_list$mHInv <- -HInv
    # }
  }


  # handle offset
  if (!is.null(data_object$offset)) {
    w_and_H_list$w <- w_and_H_list$w - data_object$offset
  }

  w_and_H_list
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
    k1 <- dispersion * (trigamma(dispersion * exp(w) / one_expw) + trigamma(dispersion / one_expw)) - 2 * sinh(w) * (k0 + 2 * atanh(1 - 2 * y))
    D_vec <- -dispersion * exp(2 * w) * k1 / one_expw^4
  }
  D <- Diagonal(x = D_vec)
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

# smw_HInv <- function(AInv, U, CInv) {
#   mid <- CInv + t(U) %*% AInv %*% U
#   # solve_mid <- tryCatch(solve(mid), error = function(e) {
#   #   diag(mid) <- diag(mid) + 1e-4 # inverse stability
#   #   solve(mid)
#   # })
#   # diag(mid) <- diag(mid) + 1e-4
#   # if (all(mid == 0)) diag(mid) <- diag(mid) + 1e-4
#   AInv - (AInv %*% U) %*% solve(mid) %*% (t(U) %*% AInv)
# }
#
# smw_mHldet <- function(A_list, AInv, U, C, CInv) {
#   Aldet <- sum(unlist(lapply(A_list, function(x) determinant(x, logarithm = TRUE)$modulus))) # must be positive det for -H
#   Cldet <- 2 * sum(log(diag(t(chol(C)))))
#   mid <- CInv + t(U) %*% AInv %*% U
#   # diag(mid) <- diag(mid) + 1e-4
#   midldet <- determinant(mid, logarithm = TRUE)$modulus
#   as.numeric(Aldet + Cldet + midldet)
# }
#
# get_DSigInv <- function(D, SigInv) {
#   D - SigInv
# }
#
# get_DSigInv_parallel <- function(cluster_list) {
#   D <- cluster_list$D
#   S <- cluster_list$S
#   get_DSigInv(D, S)
# }
