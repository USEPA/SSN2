skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"))

# Frozen SSN2 0.4.0 functions keep the comparison independent of the new helpers.
v0.4.0_reference <- new.env(parent = asNamespace("SSN2"))

v0.4.0_reference$laploglik_products <-
function(params_object, data_object, estmethod) {
    cov_matrix_list <- get_cov_matrix_list(params_object, data_object)
    if (data_object$parallel) {
        cluster_list <- lapply(seq_along(cov_matrix_list), function(l) {
            cluster_list_element <- list(c = cov_matrix_list[[l]], x = data_object$X_list[[l]], y = data_object$y_list[[l]])
        })
        cholprods_list <- parallel::parLapply(data_object$cl, cluster_list, get_cholprods_glm_parallel)
        names(cholprods_list) <- names(cov_matrix_list)
    }
    else {
        cholprods_list <- mapply(c = cov_matrix_list, x = data_object$X_list, y = data_object$y_list, 
            function(c, x, y) get_cholprods_glm(c, x, y), SIMPLIFY = FALSE)
    }
    SigInv_list <- lapply(cholprods_list, function(x) x$SigInv)
    SigInv <- Matrix::bdiag(SigInv_list)
    SigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SigInv_X))
    SqrtSigInv_X <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_X))
    SqrtSigInv_y <- do.call("rbind", lapply(cholprods_list, function(x) x$SqrtSigInv_y))
    Xt_SigInv_X <- crossprod(SqrtSigInv_X, SqrtSigInv_X)
    Xt_SigInv_X_upchol <- chol(Xt_SigInv_X)
    cov_betahat <- chol2inv(Xt_SigInv_X_upchol)
    dispersion <- as.vector(params_object$dispersion)
    w_and_H <- get_w_and_H(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, Xt_SigInv_X, 
        estmethod)
    w <- w_and_H$w
    mHldet <- w_and_H$mHldet
    betahat <- tcrossprod(cov_betahat, SigInv_X) %*% w
    X <- do.call("rbind", data_object$X_list)
    r <- w - X %*% betahat
    rt_SigInv_r <- crossprod(r, SigInv) %*% r
    y <- as.vector(do.call("rbind", data_object$y_list))
    if (!is.null(data_object$offset)) {
        w <- w + data_object$offset
    }
    l00 <- get_l00(data_object$family, w, y, data_object$size, dispersion)
    l01 <- mHldet
    l1 <- sum(unlist(lapply(cholprods_list, function(x) 2 * sum(log(diag(x$Sig_lowchol))))))
    l2 <- as.numeric(rt_SigInv_r)
    if (estmethod == "reml") {
        l3 <- 2 * sum(log(diag(Xt_SigInv_X_upchol)))
        return(list(l00 = l00, l01 = l01, l1 = l1, l2 = l2, l3 = l3))
    }
    if (estmethod == "ml") {
        return(list(l00 = l00, l01 = l01, l1 = l1, l2 = l2))
    }
}
environment(v0.4.0_reference$laploglik_products) <- v0.4.0_reference

v0.4.0_reference$get_w_and_H <-
function(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
    family <- data_object$family
    SigInv <- Matrix::bdiag(SigInv_list)
    Ptheta <- SigInv - SigInv_X %*% tcrossprod(cov_betahat, SigInv_X)
    y <- as.vector(do.call("rbind", data_object$y_list))
    size <- data_object$size
    w <- get_w_init(family, y, dispersion)
    wdiffmax <- Inf
    iter <- 0
    offset <- if (is.null(data_object$offset)) 
        0
    else as.vector(data_object$offset)
    if (length(SigInv_list) == 1) {
        while (iter < 50 && wdiffmax > 1e-04) {
            iter <- iter + 1
            d <- get_d(family, w + offset, y, size, dispersion)
            g <- d - Ptheta %*% w
            D <- get_D(family, w + offset, y, size, dispersion)
            H <- D - Ptheta
            solveHg <- solve(H, g)
            wnew <- w - solveHg
            dnew <- get_d(family, wnew + offset, y, size, dispersion)
            gnew <- dnew - Ptheta %*% wnew
            if (any(is.na(gnew) | is.infinite(gnew))) 
                stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing nugget at a known, non-zero value (via nugget_initial), or fixing dispersion at one (via dispersion_initial).", 
                  call. = FALSE)
            if (max(abs(gnew)) > max(abs(g))) 
                wnew <- w - 0.1 * solveHg
            wdiffmax <- max(abs(wnew - w))
            w <- wnew
        }
        mHldet <- as.numeric(determinant(if (estmethod == "ml") SigInv - D else -H, logarithm = TRUE)$modulus)
        w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
        if (ret_mHInv) {
            HInv <- solve(H)
            w_and_H_list$mHInv <- -HInv
        }
    }
    else {
        diag(cov_betahat_Inv) <- diag(cov_betahat_Inv) + data_object$diagtol
        while (iter < 50 && wdiffmax > 1e-04) {
            iter <- iter + 1
            d <- get_d(family, w + offset, y, size, dispersion)
            g <- d - Ptheta %*% w
            D <- get_D(family, w + offset, y, size, dispersion)
            D_diag <- diag(D)
            D_list <- lapply(split(D_diag, sort(data_object$local_index)), function(x) Diagonal(x = x))
            if (data_object$parallel) {
                cluster_list <- lapply(seq_along(D_list), function(l) {
                  cluster_list_element <- list(D = D_list[[l]], S = SigInv_list[[l]])
                })
                DSigInv_list <- parallel::parLapply(data_object$cl, cluster_list, get_DSigInv_parallel)
                names(DSigInv_list) <- names(D_list)
            }
            else {
                DSigInv_list <- mapply(D = D_list, S = SigInv_list, function(D, S) get_DSigInv(D, S), 
                  SIMPLIFY = FALSE)
            }
            if (data_object$parallel) {
                cluster_list <- DSigInv_list
                DSigInv_Inv_list <- parallel::parLapply(data_object$cl, cluster_list, solve)
                names(DSigInv_Inv_list) <- names(D_list)
            }
            else {
                DSigInv_Inv_list <- lapply(DSigInv_list, function(x) solve(x))
            }
            DSigInv_Inv <- Matrix::bdiag(DSigInv_Inv_list)
            HInv <- smw_HInv(AInv = DSigInv_Inv, U = SigInv_X, CInv = cov_betahat_Inv)
            solveHg <- HInv %*% g
            wnew <- w - solveHg
            dnew <- get_d(family, wnew + offset, y, size, dispersion)
            gnew <- dnew - Ptheta %*% wnew
            if (any(is.na(gnew) | is.infinite(gnew))) 
                stop("Convergence problem. Try using a different family, removing extreme observations, rescaling the response variable (if continuous), fixing nugget at a known, non-zero value (via nugget_initial), or fixing dispersion at one (via dispersion_initial).", 
                  call. = FALSE)
            if (max(abs(gnew)) > max(abs(g))) 
                wnew <- w - 0.1 * solveHg
            wdiffmax <- max(abs(wnew - w))
            w <- wnew
        }
        if (estmethod == "ml") {
            mHldet <- sum(vapply(DSigInv_list, function(x) {
                as.numeric(determinant(x, logarithm = TRUE)$modulus)
            }, numeric(1)))
        }
        else {
            mHldet <- smw_mHldet(A_list = DSigInv_list, AInv = DSigInv_Inv, U = SigInv_X, C = cov_betahat, 
                CInv = cov_betahat_Inv)
        }
        w_and_H_list <- list(w = w, H = NULL, mHldet = mHldet)
        if (ret_mHInv) {
            w_and_H_list$mHInv <- -HInv
        }
    }
    w_and_H_list
}
environment(v0.4.0_reference$get_w_and_H) <- v0.4.0_reference

v0.4.0_reference$get_d <-
function(family, w, y, size, dispersion) {
    if (family == "poisson") {
        d <- -exp(w) + y
    }
    else if (family == "nbinomial") {
        d <- dispersion * (y - exp(w))/(dispersion + exp(w))
    }
    else if (family == "binomial") {
        d <- y - size * expit(w)
    }
    else if (family == "Gamma") {
        d <- -dispersion + dispersion * y * exp(-w)
    }
    else if (family == "inverse.gaussian") {
        d <- dispersion * (y/(2 * exp(w)) - exp(w)/(2 * y)) + 1/2
    }
    else if (family == "beta") {
        one_expw <- 1 + exp(w)
        k0 <- digamma(dispersion * exp(w)/one_expw) - digamma(dispersion/one_expw) + log(1/y - 1)
        d <- -dispersion * exp(w) * k0/one_expw^2
    }
    d
}
environment(v0.4.0_reference$get_d) <- v0.4.0_reference

v0.4.0_reference$get_D <-
function(family, w, y, size, dispersion) {
    w <- as.vector(w)
    if (family == "poisson") {
        D_vec <- -exp(w)
    }
    else if (family == "nbinomial") {
        D_vec <- -(dispersion * exp(w) * (dispersion + y))/((dispersion + exp(w))^2)
    }
    else if (family == "binomial") {
        D_vec <- -size * expit(w)/(1 + exp(w))
    }
    else if (family == "Gamma") {
        D_vec <- -dispersion * y * exp(-w)
    }
    else if (family == "inverse.gaussian") {
        D_vec <- -dispersion * (exp(2 * w) + y^2)/(2 * y * exp(w))
    }
    else if (family == "beta") {
        one_expw <- 1 + exp(w)
        k0 <- digamma(dispersion * exp(w)/one_expw) - digamma(dispersion/one_expw) + log(1/y - 1)
        k1 <- dispersion * (trigamma(dispersion * exp(w)/one_expw) + trigamma(dispersion/one_expw)) - 
            2 * sinh(w) * k0
        D_vec <- -dispersion * exp(2 * w) * k1/one_expw^4
    }
    D <- Diagonal(x = D_vec)
}
environment(v0.4.0_reference$get_D) <- v0.4.0_reference

v0.4.0_reference$get_D3 <-
function(family, w, y, size, dispersion) {
    w <- as.numeric(w)
    mu <- exp(w)
    if (family == "poisson") {
        value <- -mu
    }
    else if (family == "nbinomial") {
        p <- stats::plogis(w - log(dispersion))
        value <- -(dispersion + y) * p * (1 - p) * (1 - 2 * p)
    }
    else if (family == "binomial") {
        p <- expit(w)
        value <- -size * p * (1 - p) * (1 - 2 * p)
    }
    else if (family == "Gamma") {
        value <- dispersion * y * exp(-w)
    }
    else if (family == "inverse.gaussian") {
        value <- -dispersion * (mu - y^2 * exp(-w))/(2 * y)
    }
    else if (family == "beta") {
        p <- expit(w)
        a <- p * (1 - p)
        b <- 1 - 2 * p
        k0 <- digamma(dispersion * p) - digamma(dispersion * (1 - p)) + log((1 - y)/y)
        k1 <- trigamma(dispersion * p) + trigamma(dispersion * (1 - p))
        k2 <- psigamma(dispersion * p, 2) - psigamma(dispersion * (1 - p), 2)
        value <- -dispersion * a * (b^2 - 2 * a) * k0 - 3 * dispersion^2 * a^2 * b * k1 - dispersion^3 * 
            a^3 * k2
    }
    as.numeric(value)
}
environment(v0.4.0_reference$get_D3) <- v0.4.0_reference

v0.4.0_reference$get_laplace_correction <-
function(mHInv, family, eta, y, size, dispersion) {
    v <- Matrix::diag(mHInv)
    t3 <- get_D3(family, eta, y, size, dispersion)
    delta <- 0.5 * as.numeric(mHInv %*% (t3 * v))
    list(delta = delta)
}
environment(v0.4.0_reference$get_laplace_correction) <- v0.4.0_reference

v0.4.0_reference$get_w_init <-
function(family, y, dispersion) {
    if (family == "poisson") {
        w_init <- 0.5 * log(y + 1)
    }
    else if (family == "nbinomial") {
        w_init <- 0.5 * log(y + 1)
    }
    else if (family == "binomial") {
        w_init <- rep(0, times = length(y))
    }
    else if (family == "Gamma") {
        w_init <- 0.5 * log(y + 1)
    }
    else if (family == "inverse.gaussian") {
        w_init <- 0.5 * log(y + 1)
    }
    else if (family == "beta") {
        w_init <- rep(0, times = length(y))
    }
    w_init
}
environment(v0.4.0_reference$get_w_init) <- v0.4.0_reference

v0.4.0_reference$get_l00 <-
function(family, w, y, size, dispersion) {
    w <- as.vector(w)
    y <- as.vector(y)
    if (family == "poisson") {
        mu <- exp(w)
        l00 <- -2 * sum(dpois(y, lambda = mu, log = TRUE))
    }
    else if (family == "nbinomial") {
        mu <- exp(w)
        l00 <- -2 * sum(dnbinom(x = y, mu = mu, size = dispersion, log = TRUE))
    }
    else if (family == "binomial") {
        mu <- expit(w)
        l00 <- -2 * sum(dbinom(y, size, mu, log = TRUE))
    }
    else if (family == "Gamma") {
        mu <- exp(w)
        l00 <- -2 * sum(dgamma(y, shape = dispersion, scale = mu/dispersion, log = TRUE))
    }
    else if (family == "inverse.gaussian") {
        mu <- exp(w)
        l00 <- -2 * sum(1/2 * (log(dispersion) + log(exp(w)) - log(2 * pi) - log(y^3)) - dispersion * 
            (y - exp(w))^2/(2 * exp(w) * y))
    }
    else if (family == "beta") {
        mu <- expit(w)
        a <- mu * dispersion
        b <- (1 - mu) * dispersion
        l00 <- -2 * sum(dbeta(x = y, shape1 = a, shape2 = b, log = TRUE))
    }
    l00
}
environment(v0.4.0_reference$get_l00) <- v0.4.0_reference

v0.4.0_reference$smw_HInv <-
function(AInv, U, CInv) {
    mid <- CInv + t(U) %*% AInv %*% U
    AInv - (AInv %*% U) %*% solve(mid) %*% (t(U) %*% AInv)
}
environment(v0.4.0_reference$smw_HInv) <- v0.4.0_reference

v0.4.0_reference$smw_mHldet <-
function(A_list, AInv, U, C, CInv) {
    Aldet <- sum(unlist(lapply(A_list, function(x) determinant(x, logarithm = TRUE)$modulus)))
    Cldet <- 2 * sum(log(diag(t(chol(C)))))
    mid <- CInv + t(U) %*% AInv %*% U
    midldet <- determinant(mid, logarithm = TRUE)$modulus
    as.numeric(Aldet + Cldet + midldet)
}
environment(v0.4.0_reference$smw_mHldet) <- v0.4.0_reference

v0.4.0_reference$get_DSigInv <-
function(D, SigInv) {
    D - SigInv
}
environment(v0.4.0_reference$get_DSigInv) <- v0.4.0_reference

v0.4.0_reference$get_DSigInv_parallel <-
function(cluster_list) {
    D <- cluster_list$D
    S <- cluster_list$S
    get_DSigInv(D, S)
}
environment(v0.4.0_reference$get_DSigInv_parallel) <- v0.4.0_reference

v0.4.0_reference$get_model_stats_glm_core <-
function(cov_est_object, data_object, estmethod, order, eigenprods_list) {
    X <- do.call("rbind", data_object$X_list)
    y <- do.call("rbind", data_object$y_list)
    SigInv_list <- lapply(eigenprods_list, function(x) x$SigInv)
    SigInv <- Matrix::bdiag(SigInv_list)
    SigInv_X <- do.call("rbind", lapply(eigenprods_list, function(x) x$SigInv_X))
    invcov_betahat_list <- lapply(eigenprods_list, function(x) crossprod(x$SqrtSigInv_X, x$SqrtSigInv_X))
    invcov_betahat_sum <- Reduce("+", invcov_betahat_list)
    cov_betahat_noadjust <- chol2inv(chol(forceSymmetric(invcov_betahat_sum)))
    cov_betahat_noadjust_list <- rep(list(cov_betahat_noadjust), times = length(invcov_betahat_list))
    dispersion <- as.vector(cov_est_object$params_object$dispersion)
    w_and_H <- get_w_and_H(data_object, dispersion, SigInv_list, SigInv_X, cov_betahat_noadjust, invcov_betahat_sum, 
        estmethod, ret_mHInv = TRUE)
    w <- as.vector(w_and_H$w)
    w_list <- split(w, sort(data_object$local_index))
    Xt_SigInv_w_list <- mapply(x = eigenprods_list, w = w_list, function(x, w) crossprod(x$SigInv_X, 
        w), SIMPLIFY = FALSE)
    betahat_list <- mapply(l = cov_betahat_noadjust_list, r = Xt_SigInv_w_list, function(l, r) l %*% 
        r, SIMPLIFY = FALSE)
    betahat <- as.numeric(cov_betahat_noadjust %*% Reduce("+", Xt_SigInv_w_list))
    names(betahat) <- colnames(data_object$X_list[[1]])
    eta <- if (is.null(data_object$offset)) 
        w
    else w + data_object$offset
    laplace <- get_laplace_correction(w_and_H$mHInv, data_object$family, eta, y, data_object$size, dispersion)
    center_list <- split(w + laplace$delta, sort(data_object$local_index))
    center_products <- mapply(function(x, w) crossprod(x$SigInv_X, w), eigenprods_list, center_list, 
        SIMPLIFY = FALSE)
    laplace$betahat <- as.numeric(cov_betahat_noadjust %*% Reduce("+", center_products))
    names(laplace$betahat) <- names(betahat)
    cov_betahat <- cov_betahat_adjust(invcov_betahat_list, betahat_list, betahat, eigenprods_list, data_object, 
        cov_est_object$params_object, cov_betahat_noadjust, data_object$var_adjust)
    cov_betahat <- as.matrix(cov_betahat)
    wts_beta <- tcrossprod(cov_betahat_noadjust, SigInv_X)
    betawtsvarw <- wts_beta %*% w_and_H$mHInv %*% t(wts_beta)
    cov_betahat_uncorrected <- cov_betahat
    rownames(cov_betahat_uncorrected) <- colnames(data_object$X_list[[1]])
    colnames(cov_betahat_uncorrected) <- colnames(data_object$X_list[[1]])
    cov_betahat <- as.matrix(cov_betahat + betawtsvarw)
    rownames(cov_betahat) <- colnames(data_object$X_list[[1]])
    colnames(cov_betahat) <- colnames(data_object$X_list[[1]])
    coefficients <- get_coefficients_glm(betahat, cov_est_object$params_object)
    fitted <- get_fitted_glm(w_list, betahat, cov_est_object$params_object, data_object, eigenprods_list)
    hatvalues <- get_hatvalues_glm(w, X, data_object, dispersion)
    deviance_i <- get_deviance_glm(data_object$family, y, fitted$response, data_object$size, dispersion)
    deviance_i <- pmax(deviance_i, 0)
    SigInv_X_null <- do.call("rbind", lapply(eigenprods_list, function(x) x$SigInv_ones))
    SqrtSigInv_X_null <- do.call("rbind", lapply(eigenprods_list, function(x) x$SqrtSigInv_ones))
    Xt_SigInv_X_null <- crossprod(SqrtSigInv_X_null, SqrtSigInv_X_null)
    Xt_SigInv_X_upchol_null <- chol(Xt_SigInv_X_null)
    cov_betahat_null <- chol2inv(Xt_SigInv_X_upchol_null)
    w_and_H_null <- get_w_and_H(data_object, dispersion, SigInv_list, SigInv_X_null, cov_betahat_null, 
        Xt_SigInv_X_null, estmethod)
    w_null <- as.vector(w_and_H_null$w)
    fitted_null <- get_fitted_null(w_null, data_object)
    deviance_i_null <- get_deviance_glm(data_object$family, y, fitted_null, data_object$size, dispersion)
    deviance_i_null <- pmax(deviance_i_null, 0)
    deviance <- sum(deviance_i)
    deviance_null <- sum(deviance_i_null)
    pseudoR2 <- as.numeric(1 - deviance/deviance_null)
    pseudoR2 <- pmax(0, pseudoR2)
    if (length(labels(terms(data_object$formula))) == 0) {
        pseudoR2 <- 0
    }
    residuals <- get_residuals_glm(w, y, data_object, deviance_i, hatvalues, dispersion)
    cooks_distance <- get_cooks_distance_glm(residuals, hatvalues, data_object$p)
    vcov <- get_vcov_glm(cov_betahat, cov_betahat_uncorrected)
    model_stats_glm_names <- data_object$pid[data_object$observed_index]
    laplace$delta <- laplace$delta[order(order)]
    names(laplace$delta) <- model_stats_glm_names
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
    hatvalues <- hatvalues[order(order)]
    names(hatvalues) <- model_stats_glm_names
    residuals$response <- residuals$response[order(order)]
    names(residuals$response) <- model_stats_glm_names
    residuals$deviance <- residuals$deviance[order(order)]
    names(residuals$deviance) <- model_stats_glm_names
    residuals$pearson <- residuals$pearson[order(order)]
    names(residuals$pearson) <- model_stats_glm_names
    residuals$standardized <- residuals$standardized[order(order)]
    names(residuals$standardized) <- model_stats_glm_names
    cooks_distance <- cooks_distance[order(order)]
    names(cooks_distance) <- model_stats_glm_names
    y <- y[order(order)]
    if (is.null(data_object$size)) {
        size <- NULL
    }
    else {
        size <- data_object$size[order(order)]
        names(size) <- model_stats_glm_names
    }
    npar <- sum(unlist(lapply(cov_est_object$is_known, function(x) length(x) - sum(x))))
    list(coefficients = coefficients, fitted = fitted, hatvalues = hatvalues, residuals = residuals, 
        cooks_distance = cooks_distance, vcov = vcov, deviance = deviance, pseudoR2 = pseudoR2, laplace = laplace, 
        npar = npar, w = w, y = y, size = size)
}
environment(v0.4.0_reference$get_model_stats_glm_core) <- v0.4.0_reference

v0.4.0_data <- function(n = 36L) {
  set.seed(2)
  s <- mf04p
  root <- tempfile('m1-newton-')
  dir.create(root)
  stopifnot(file.copy(s$path, root, recursive = TRUE, copy.mode = FALSE))
  s$path <- file.path(root, basename(s$path))
  stopifnot(n <= nrow(s$obs))
  s$obs <- s$obs[seq_len(n), ]
  s$preds <- list(newton = s$preds$pred1km[c(5, 1, 8), ])
  d <- data.frame(x = rnorm(n), x2 = rnorm(n),
    x3 = 1e3 * rnorm(n), off = runif(n, -.15, .15),
    part = factor(rep(c('z', 'b'), length.out = n)), g = factor(rep(1:8, length.out = n)))
  eta <- .2 + .25 * d$x + rnorm(n, sd = .25)
  mu <- exp(eta); pr <- plogis(eta)
  d$poisson <- rpois(n, mu)
  d$nbinomial <- rnbinom(n, mu = mu, size = 4)
  d$binomial <- rbinom(n, 5, pr)
  d$bernoulli <- rbinom(n, 1, pr)
  d$Gamma <- rgamma(n, shape = 4, rate = 4 / mu)
  d$inverse.gaussian <- mu * exp(rnorm(n, sd = .25))
  d$beta <- .2 + .6 * pr
  for (nm in names(d)) s$obs[[nm]] <- d[[nm]]
  s$preds$newton$x <- c(-.5, .2, .8)
  s$preds$newton$off <- c(.1, -.3, .2)
  s$preds$newton$g <- factor(c(1, 2, 3), levels = levels(d$g))
  s$preds$newton$part <- factor(c('z', 'b', 'z'), levels = levels(d$part))
  ssn_create_distmat(s, predpts = 'newton', among_predpts = TRUE, overwrite = TRUE)
  s
}

v0.4.0_cases <- function() {
  d <- v0.4.0_data()
  out <- list()
  for (fam in c('poisson', 'binomial', 'bernoulli', 'nbinomial', 'Gamma', 'inverse.gaussian', 'beta')) {
    for (layout in c('single', 'local', 'partition')) for (method in c('ml', 'reml')) {
      out[[paste(fam, layout, method, sep = '-')]] <- v0.4.0_args(d, fam, layout, method)
    }
  }
  for (method in c('ml', 'reml')) {
    a <- v0.4.0_args(d, 'poisson', 'local', method)
    a$formula <- poisson ~ x + x2 + x3 + offset(off)
    a$local <- list(index = sample(rep(c(7, 2, 9, 4), c(1, 5, 11, 19))))
    out[[paste0('unequal-singleton-scaled-', method)]] <- a
    a <- v0.4.0_args(d, 'poisson', 'local', method)
    a$formula <- poisson ~ x + offset(off)
    a$random <- ~g
    a$randcov_initial <- randcov_initial(g = .1, known = 'given')
    a$partition_factor <- ~part
    out[[paste0('random-partition-offset-', method)]] <- a
    a <- v0.4.0_args(v0.4.0_data(36L), 'bernoulli', 'local', method)
    a$local <- list(index = rep(c(7, 2, 9, 4), each = 9))
    a$nugget_initial <- nugget_initial('nugget', 1e-4, known = 'given')
    out[[paste0('bernoulli-floor-', method)]] <- a
  }
  out
}

v0.4.0_args <- function(d, fam, layout, method, stream = 'mixed') {
  response <- if (fam == 'binomial') 'cbind(binomial, 5 - binomial)' else fam
  f <- as.formula(paste(response, '~ x'))
  a <- list(formula = f, ssn.object = d, family = if (fam == 'bernoulli') 'binomial' else fam,
    additive = 'afvArea', estmethod = method,
    nugget_initial = nugget_initial('nugget', .2, known = 'given'),
    local = if (layout == 'local') list(index = sample(rep(c(7, 2, 9, 4), length.out = nrow(d$obs)))) else FALSE)
  if (stream %in% c('tailup', 'mixed'))
    a$tailup_initial <- tailup_initial('exponential', .2, 10000, known = 'given')
  if (stream %in% c('taildown', 'mixed'))
    a$taildown_initial <- taildown_initial('exponential', .2, 10000, known = 'given')
  if (stream %in% c('euclidean', 'mixed'))
    a$euclid_initial <- euclid_initial('exponential', .4, 12000, known = 'given')
  if (layout == 'partition') a$partition_factor <- ~part
  if (!fam %in% c('poisson', 'binomial', 'bernoulli'))
    a$dispersion_initial <- dispersion_initial(as.character(fam), 4, known = 'given')
  a
}

v0.4.0_capture <- function(args, solver = get_w_and_H) {
  force(solver)
  captured <- NULL
  fit <- testthat::with_mocked_bindings(do.call(ssn_glm, args),
    get_w_and_H = function(data_object, dispersion, SigInv_list, SigInv_X,
      cov_betahat, cov_betahat_Inv, estmethod, ret_mHInv = FALSE) {
      a <- as.list(environment())
      a$captured <- a$solver <- a$args <- a$fit <- NULL
      if (ret_mHInv) captured <<- a
      do.call(solver, a)
    }, .package = 'SSN2')
  list(fit = fit, args = captured)
}

v0.4.0_outputs <- function(fit) {
  list(coefficients = fit$coefficients, logLik = as.numeric(logLik(fit)), AIC = AIC(fit),
    deviance = deviance(fit), vcov = vcov(fit),
    link = fitted(fit, type = 'link'), response = fitted(fit, type = 'response'),
    uncorrected_vcov = vcov(fit, var_correct = FALSE),
    laplace = fit$laplace, prediction = predict(fit, 'newton', local = FALSE, se.fit = TRUE))
}

v0.4.0_scaled_error <- function(new, ref) {
  new <- as.numeric(new)
  ref <- as.numeric(ref)
  if (length(new) != length(ref)) return(Inf)
  finite <- is.finite(new) & is.finite(ref)
  if (!identical(new[!finite], ref[!finite])) return(Inf)
  if (!any(finite)) return(0)
  max(abs(new[finite] - ref[finite])) / max(1, max(abs(ref[finite])))
}

v0.4.0_record <- function(table, row) {
  path <- Sys.getenv('SSN2_v0.4.0_AGREEMENT_DIR')
  if (!nzchar(path)) return(invisible(NULL))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(path, paste0(table, '.csv'))
  utils::write.table(row, file, sep = ',', row.names = FALSE,
    col.names = !file.exists(file), append = file.exists(file))
}

v0.4.0_check_outputs <- function(new, old, case, table) {
  new_out <- v0.4.0_outputs(new)
  old_out <- v0.4.0_outputs(old)
  for (field in names(old_out)) {
    error <- v0.4.0_scaled_error(unlist(new_out[[field]]), unlist(old_out[[field]]))
    v0.4.0_record(table, data.frame(case = case, field = field, scaled_error = error,
      logLik_absolute_error = abs(as.numeric(logLik(new) - logLik(old)))))
    expect_lte(error, 1e-6)
  }
}

v0.4.0_solve <- function(solver, a) {
  force(solver)
  calls <- 0L
  last_D <- NULL
  original_d <- get('get_d', envir = environment(solver))
  original_D <- get('get_D', envir = environment(solver))
  localenv <- new.env(parent = environment(solver))
  localenv$get_d <- function(...) { calls <<- calls + 1L; original_d(...) }
  localenv$get_D <- function(...) { last_D <<- original_D(...); last_D }
  environment(solver) <- localenv
  value <- do.call(solver, a)
  list(value = value, iterations = calls / 2L, D = last_D)
}

v0.4.0_check_solver <- function(a, solver = get_w_and_H, case = 'solver') {
  old <- v0.4.0_solve(v0.4.0_reference$get_w_and_H, a)
  new <- v0.4.0_solve(solver, a)
  v0.4.0_record('solver-agreement', data.frame(case = case, family = a$data_object$family,
    method = a$estmethod, blocks = length(a$SigInv_list), inverse = a$ret_mHInv,
    iterations_new = new$iterations, iterations_ref = old$iterations,
    w_error = v0.4.0_scaled_error(new$value$w, old$value$w),
    determinant_error = v0.4.0_scaled_error(new$value$mHldet, old$value$mHldet)))
  expect_identical(new$iterations, old$iterations)
  expect_lte(v0.4.0_scaled_error(new$value$w, old$value$w), 1e-10)
  expect_lte(v0.4.0_scaled_error(new$value$mHldet, old$value$mHldet), 1e-8)
  if (isTRUE(a$ret_mHInv)) {
    M <- new$value$mHInv
    n <- length(new$value$w)
    expect_lte(v0.4.0_scaled_error(mhinv_apply(M, diag(n)), old$value$mHInv), 1e-10)
    expect_lte(v0.4.0_scaled_error(mhinv_diag(M), Matrix::diag(old$value$mHInv)), 1e-10)
    V <- cbind(seq_along(new$value$w) / length(new$value$w), 1)
    expect_lte(v0.4.0_scaled_error(mhinv_apply(M, V), old$value$mHInv %*% V), 1e-10)
    expect_lte(v0.4.0_scaled_error(mhinv_apply(M, V[, 1]), old$value$mHInv %*% V[, 1]), 1e-10)
    eta <- as.vector(new$value$w) + if (is.null(a$data_object$offset)) 0 else a$data_object$offset
    y <- as.vector(do.call('rbind', a$data_object$y_list))
    delta <- get_laplace_correction(M, a$data_object$family, eta, y, a$data_object$size, a$dispersion)
    ref <- v0.4.0_reference$get_laplace_correction(old$value$mHInv, a$data_object$family, eta, y,
      a$data_object$size, a$dispersion)
    expect_lte(v0.4.0_scaled_error(delta$delta, ref$delta), 1e-10)
    W <- as.matrix(tcrossprod(a$cov_betahat, a$SigInv_X))
    expect_lte(v0.4.0_scaled_error(W %*% mhinv_apply(M, t(W)), W %*% old$value$mHInv %*% t(W)), 1e-10)
    v0.4.0_record('consumer-agreement', data.frame(case = case, family = a$data_object$family,
      method = a$estmethod, blocks = length(a$SigInv_list),
      inverse_error = v0.4.0_scaled_error(mhinv_apply(M, diag(n)), old$value$mHInv),
      diagonal_error = v0.4.0_scaled_error(mhinv_diag(M), Matrix::diag(old$value$mHInv)),
      matrix_error = v0.4.0_scaled_error(mhinv_apply(M, V), old$value$mHInv %*% V),
      vector_error = v0.4.0_scaled_error(mhinv_apply(M, V[, 1]), old$value$mHInv %*% V[, 1]),
      correction_error = v0.4.0_scaled_error(delta$delta, ref$delta),
      coefficient_covariance_error = v0.4.0_scaled_error(W %*% mhinv_apply(M, t(W)),
        W %*% old$value$mHInv %*% t(W))))
    if (length(a$SigInv_list) > 1L) {
      expect_s3_class(M, 'laplace_mhinv')
      expect_named(M, c('F', 'idx', 'AU', 'mid_inv'))
      expect_equal(dim(M$AU), c(n, ncol(a$SigInv_X)))
    }
  }
  invisible(list(old = old, new = new))
}

test_that('block Newton solves preserve all families, methods and model terms', {
  set.seed(2)
  cases <- v0.4.0_cases()
  for (case in names(cases)) {
    args <- cases[[case]]
    captured <- v0.4.0_capture(args)
    v0.4.0_check_solver(captured$args, case = case)
    a <- captured$args
    ends <- cumsum(vapply(a$SigInv_list, nrow, 1L))
    idx <- Map(seq.int, c(1L, head(ends, -1L) + 1L), ends)
    r <- seq_len(nrow(a$SigInv_X)) / 10
    block <- sum(mapply(function(S, i) as.numeric(crossprod(r[i], S %*% r[i])), a$SigInv_list, idx))
    dense <- as.numeric(crossprod(r, Matrix::bdiag(a$SigInv_list)) %*% r)
    expect_lte(v0.4.0_scaled_error(block, dense), 1e-8)
    params <- captured$fit$coefficients$params_object
    old <- v0.4.0_reference$laploglik_products(params, a$data_object, a$estmethod)
    new <- laploglik_products(params, a$data_object, a$estmethod)
    for (field in names(old)) {
      error <- v0.4.0_scaled_error(new[[field]], old[[field]])
      v0.4.0_record('likelihood-agreement', data.frame(case = case, field = field, scaled_error = error))
      expect_lte(error, 1e-8)
    }
    a$ret_mHInv <- FALSE
    v0.4.0_check_solver(a, case = case)
  }
})

test_that('stream covariance paths preserve known fits and every inverse consumer', {
  s <- v0.4.0_data(24L)
  for (stream in c('tailup', 'taildown', 'euclidean', 'mixed')) {
    for (layout in c('single', 'local')) for (method in c('ml', 'reml')) {
      args <- v0.4.0_args(s, 'poisson', layout, method, stream)
      args$formula <- poisson ~ x + offset(off)
      args$random <- ~g
      args$randcov_initial <- randcov_initial(g = .1, known = 'given')
      args$partition_factor <- ~part
      new <- v0.4.0_capture(args)
      old <- with_mocked_bindings(v0.4.0_capture(args, v0.4.0_reference$get_w_and_H),
        laploglik_products = v0.4.0_reference$laploglik_products,
        get_laplace_correction = v0.4.0_reference$get_laplace_correction,
        get_model_stats_glm_core = v0.4.0_reference$get_model_stats_glm_core)
      case <- paste(stream, layout, method, sep = '-')
      v0.4.0_check_solver(new$args, case = case)
      v0.4.0_check_outputs(new$fit, old$fit, case, 'known-fit-agreement')
      for (correct in c(FALSE, TRUE)) {
        expect_lte(v0.4.0_scaled_error(
          unlist(loocv(new$fit, cv_predict = TRUE, se.fit = TRUE, fit_correct = correct)),
          unlist(loocv(old$fit, cv_predict = TRUE, se.fit = TRUE, fit_correct = correct))), 1e-6)
        expect_lte(v0.4.0_scaled_error(
          unlist(kcv(new$fit, k = 3, cv_predict = TRUE, se.fit = TRUE, fit_correct = correct,
            folds_index = rep(1:3, length.out = 24))),
          unlist(kcv(old$fit, k = 3, cv_predict = TRUE, se.fit = TRUE, fit_correct = correct,
            folds_index = rep(1:3, length.out = 24)))), 1e-6)
      }
      set.seed(2)
      draws <- conditional(new$fit, output = c('latent', 'beta'), samples = 5)
      set.seed(2)
      ref <- conditional(old$fit, output = c('latent', 'beta'), samples = 5)
      expect_lte(v0.4.0_scaled_error(unlist(draws), unlist(ref)), 1e-6)
    }
  }
})

test_that('fully estimated stream fits preserve natural outputs and interior parameters', {
  s <- v0.4.0_data(24L)
  families <- c('poisson', 'binomial', 'bernoulli', 'nbinomial', 'Gamma', 'inverse.gaussian', 'beta')
  streams <- c('mixed', 'tailup', 'taildown', 'euclidean', 'tailup', 'taildown', 'euclidean')
  for (i in seq_along(families)) {
    for (layout in c('single', 'local')) for (method in c('ml', 'reml')) {
      args <- v0.4.0_args(s, families[i], layout, method, streams[i])
      args$formula <- update(args$formula, . ~ . + offset(off))
      args$nugget_initial <- nugget_initial('nugget', .2)
      if (!is.null(args$tailup_initial)) args$tailup_initial <- tailup_initial('exponential', .2, 10000)
      if (!is.null(args$taildown_initial)) args$taildown_initial <- taildown_initial('exponential', .2, 10000)
      if (!is.null(args$euclid_initial)) args$euclid_initial <- euclid_initial('exponential', .4, 12000)
      if (!is.null(args$dispersion_initial))
        args$dispersion_initial <- dispersion_initial(as.character(families[i]), 4)
      new <- v0.4.0_capture(args)
      old <- with_mocked_bindings(v0.4.0_capture(args, v0.4.0_reference$get_w_and_H),
        laploglik_products = v0.4.0_reference$laploglik_products,
        get_laplace_correction = v0.4.0_reference$get_laplace_correction,
        get_model_stats_glm_core = v0.4.0_reference$get_model_stats_glm_core)
      case <- paste(families[i], layout, method, sep = '-')
      v0.4.0_check_outputs(new$fit, old$fit, case, 'estimated-fit-agreement')
      p <- unlist(old$fit$coefficients$params_object)
      q <- unlist(new$fit$coefficients$params_object)
      expect_identical(names(q), names(p))
      boundary <- grepl('\\.de$|nugget', names(p)) & pmin(p, q) <= 1e-6
      if (any(boundary)) {
        message('Flat-boundary log parameters (', families[i], '/', layout, '/', method,
          '): old=', paste(log(p[boundary]), collapse = ','),
          '; new=', paste(log(q[boundary]), collapse = ','))
        expect_lte(abs(as.numeric(logLik(new$fit) - logLik(old$fit))), 1e-8)
      }
      interior <- !boundary & p > 0 & q > 0 & is.finite(p) & is.finite(q)
      for (j in seq_along(p)) v0.4.0_record('estimated-parameter-agreement', data.frame(
        case = case, parameter = names(p)[j], new = q[j], reference = p[j],
        log_new = log(q[j]), log_reference = log(p[j]), boundary = boundary[j],
        natural_error = v0.4.0_scaled_error(q[j], p[j]),
        log_error = v0.4.0_scaled_error(log(q[j]), log(p[j]))))
      expect_lte(v0.4.0_scaled_error(log(q[interior]), log(p[interior])), 1e-6)
      v0.4.0_check_solver(new$args, case = paste0('estimated-', case))
    }
  }
})

test_that('block Newton factorization preserves parallel fitting', {
  set.seed(2)
  args <- v0.4.0_args(v0.4.0_data(), 'poisson', 'local', 'reml')
  serial <- v0.4.0_capture(args)
  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
  parallel::clusterEvalQ(cl, library(SSN2))
  expect_true(all(unlist(parallel::clusterEvalQ(cl,
    exists('factor_laplace_block', envir = asNamespace('SSN2'), inherits = FALSE)))))
  a <- serial$args
  a$data_object$parallel <- TRUE
  a$data_object$cl <- cl
  checked <- v0.4.0_check_solver(a)
  args$local$parallel <- TRUE
  args$local$ncores <- 2L
  fit <- do.call(ssn_glm, args)
  expect_lte(v0.4.0_scaled_error(coef(fit), coef(serial$fit)), 1e-10)
  expect_lte(v0.4.0_scaled_error(vcov(fit), vcov(serial$fit)), 1e-10)
})

test_that('block factors preserve singleton and indefinite beta curvature', {
  set.seed(2)
  candidates <- expand.grid(w = seq(-8, 8, length.out = 65), y = c(.01, .1, .5, .9, .99))
  D <- diag(get_D('beta', candidates$w, candidates$y, NULL, 4))
  k <- which.max(D)
  expect_gt(D[k], 0)
  S <- list(matrix(D[k] / 2, 1, 1), matrix(c(2, .2, .2, 2), 2))
  dv <- c(D[k], -1, -2)
  idx <- list(1L, 2:3)
  F <- factor_laplace_blocks(S, dv, idx, FALSE, NULL)
  expect_identical(vapply(F, `[[`, '', 'type'), c('lu', 'chol'))
  A <- Matrix::bdiag(Map(function(s, i) v0.4.0_reference$get_DSigInv(Matrix::Diagonal(x = dv[i]), s), S, idx))
  V <- cbind(c(.5, -1, 2), 1)
  expect_lte(v0.4.0_scaled_error(apply_Ainv(F, V, idx), solve(A, V)), 1e-10)
  expect_lte(v0.4.0_scaled_error(apply_Ainv(F, V[, 1], idx), solve(A, V[, 1])), 1e-10)
  expect_lte(v0.4.0_scaled_error(logdet_A(F), as.numeric(determinant(A)$modulus)), 1e-8)
  U <- matrix(c(.1, .2, -.1), 3, 1)
  AU <- apply_Ainv(F, U, idx)
  mid <- matrix(1, 1, 1) + crossprod(U, AU)
  M <- structure(list(F = F, idx = idx, AU = AU, mid_inv = solve(mid)), class = 'laplace_mhinv')
  ref <- -as.matrix(v0.4.0_reference$smw_HInv(solve(A), U, matrix(1, 1, 1)))
  expect_lte(v0.4.0_scaled_error(mhinv_apply(M, V), ref %*% V), 1e-10)
  expect_lte(v0.4.0_scaled_error(mhinv_apply(M, V[, 1]), ref %*% V[, 1]), 1e-10)
  expect_lte(v0.4.0_scaled_error(mhinv_diag(M), diag(ref)), 1e-10)
  bad <- factor_laplace_blocks(list(matrix(1, 1, 1)), 1, list(1L), FALSE, NULL)
  expect_error(apply_Ainv(bad, matrix(1, 1, 1), list(1L)))
  expect_error(solve(v0.4.0_reference$get_DSigInv(Matrix::Diagonal(x = 1), matrix(1, 1, 1))))
})

test_that('LU block solves preserve the complete Newton iteration', {
  args <- v0.4.0_args(v0.4.0_data(), 'beta', 'local', 'reml')
  a <- v0.4.0_capture(args)$args
  count <- 0L
  with_mocked_bindings(v0.4.0_check_solver(a), factor_laplace_block = function(block) {
    count <<- count + 1L
    m <- block$S
    diag(m) <- diag(m) - block$D
    list(type = 'lu', A = -m)
  })
  expect_gt(count, 0L)
})

test_that('naturally indefinite beta iterations preserve the convergence error', {
  D <- as.numeric(diag(get_D('beta', 1.75, .01, NULL, 4)))
  expect_gt(D, 0)
  S <- list(matrix(D / 2, 1, 1), matrix(c(2, .2, .2, 2), 2))
  X <- matrix(1, 3, 1)
  U <- as.matrix(Matrix::bdiag(S) %*% X)
  CI <- crossprod(X, U)
  a <- list(data_object = list(family = 'beta',
    y_list = list(matrix(.01, 1), matrix(c(.4, .6), 2)), size = NULL,
    local_index = c(1, 2, 2), diagtol = 1e-5, parallel = FALSE),
    dispersion = 4, SigInv_list = S, SigInv_X = U,
    cov_betahat = solve(CI), cov_betahat_Inv = CI, estmethod = 'reml', ret_mHInv = TRUE)
  counts <- integer(2)
  lu <- 0L
  original <- factor_laplace_block
  with_mocked_bindings({
    solvers <- list(v0.4.0_reference$get_w_and_H, get_w_and_H)
    for (i in seq_along(solvers)) {
      solver <- solvers[[i]]
      localenv <- new.env(parent = environment(solver))
      localenv$get_w_init <- function(...) c(1.75, 0, 0)
      original_d <- get('get_d', envir = environment(solver))
      localenv$get_d <- function(...) { counts[i] <<- counts[i] + 1L; original_d(...) }
      environment(solver) <- localenv
      expect_error(suppressWarnings(do.call(solver, a)), 'Convergence problem.*nugget_initial')
    }
  }, factor_laplace_block = function(block) {
    F <- original(block)
    if (F$type == 'lu') lu <<- lu + 1L
    F
  })
  expect_gt(lu, 0L)
  expect_identical(counts[1], counts[2])
  v0.4.0_record('indefinite-agreement', data.frame(case = 'natural-beta-convergence-error',
    derivative_calls_ref = counts[1], derivative_calls_new = counts[2], lu_blocks = lu))
})

test_that('larger Bernoulli blocks preserve the solve at the nugget floor', {
  args <- v0.4.0_args(v0.4.0_data(44L), 'bernoulli', 'local', 'reml')
  args$local <- list(index = rep(c(7, 2, 9, 4), each = 11))
  args$nugget_initial <- nugget_initial('nugget', 1e-4, known = 'given')
  a <- v0.4.0_capture(args, solver = v0.4.0_reference$get_w_and_H)$args
  v0.4.0_check_solver(a)
})

test_that('inverse helpers preserve ordinary and Matrix inputs', {
  V <- cbind(c(.5, -1, 2), 1)
  dense <- crossprod(matrix(seq_len(9) / 10, 3)) + diag(3)
  for (M in list(dense, Matrix::Matrix(dense, sparse = FALSE), Matrix::Matrix(dense, sparse = TRUE))) {
    expect_identical(mhinv_diag(M), Matrix::diag(M))
    expect_identical(mhinv_apply(M, V), as.matrix(M %*% V))
    expect_identical(mhinv_apply(M, V[, 1]), as.matrix(M %*% V[, 1]))
    expect_identical(get_laplace_correction(M, 'poisson', rep(.2, 3), 1:3, NULL, 1),
      v0.4.0_reference$get_laplace_correction(M, 'poisson', rep(.2, 3), 1:3, NULL, 1))
  }
})

test_that('single-block beta solves preserve the general-solve fallback', {
  a <- v0.4.0_capture(v0.4.0_args(v0.4.0_data(), 'beta', 'single', 'reml'))$args
  calls <- 0L
  fallback <- get_w_and_H
  localenv <- new.env(parent = environment(get_w_and_H))
  localenv$chol <- function(...) {
    calls <<- calls + 1L
    stop('Forced Cholesky failure')
  }
  environment(fallback) <- localenv
  for (method in c('ml', 'reml')) {
    a$estmethod <- method
    v0.4.0_check_solver(a, fallback)
  }
  expect_gt(calls, 0L)
})
