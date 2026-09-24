skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("local kcv forwards parallel controls and preserves fitting groups", {
  object <- list(n = 8, local_index = rep(1:2, 4))
  local <- resolve_kcv_local(list(size = 3, parallel = TRUE, ncores = 2))
  controls <- get_kcv_estimation_local(object, c(2, 5), local)
  expect_true(controls$parallel)
  expect_identical(controls$ncores, 2)
  expect_identical(controls$index, object$local_index[-c(2, 5)])
  expect_identical(controls$var_adjust, "theoretical")
  object$local_index <- NULL
  expect_identical(get_kcv_estimation_local(object, 1:6, local)$size, 2)
})

test_that("kcv parallel dispatch uses workers and cleans up after errors", {
  local <- list(parallel = TRUE, ncores = 2)
  pids <- unlist(run_pred_dispatch(function(x) Sys.getpid(), as.list(1:4), local))
  expect_length(unique(pids), 2)
  expect_false(Sys.getpid() %in% pids)
  stopped <- FALSE
  stop_cluster <- parallel::stopCluster
  local_mocked_bindings(stopCluster = function(cl) {
    stopped <<- TRUE
    stop_cluster(cl)
  }, .package = "parallel")
  expect_error(run_pred_dispatch(function(x) stop("worker failure"), list(1), local),
               "worker failure")
  expect_true(stopped)
})

test_that("parallel initialization and fitting errors stop their clusters", {
  stopped <- 0L
  local_mocked_bindings(stopCluster = function(cl) stopped <<- stopped + 1L,
    makeCluster = function(...) list(),
    clusterCall = function(...) stop("initialization failure"), .package = "parallel")
  expect_error(make_ssn_cluster(2), "initialization failure")
  expect_identical(stopped, 1L)
  local_mocked_bindings(
    make_ssn_cluster = function(...) list(),
    cov_estimate_gloglik = function(...) stop("fitting failure"),
    cov_estimate_laploglik = function(...) stop("fitting failure")
  )
  ssn <- mf04p
  ssn_create_bigdist(ssn, overwrite = TRUE)
  ssn$obs$cv_count <- round(ssn$obs$Summer_mn)
  args <- list(formula = cv_count ~ ELEV_DEM, ssn.object = ssn,
    tailup_initial = tailup_initial("exponential", de = 2, range = 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.5, known = "given"),
    additive = "afvArea",
    local = list(index = rep(1:2, length.out = nrow(ssn$obs)), parallel = TRUE, ncores = 2))
  expect_error(do.call(ssn_lm, args), "fitting failure")
  expect_identical(stopped, 2L)
  expect_error(do.call(ssn_glm, c(args, list(family = "poisson"))), "fitting failure")
  expect_identical(stopped, 3L)
})

test_that("singleton kcv delegates to loocv for both model classes", {
  args <- list(formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_initial = tailup_initial("exponential", de = 2, range = 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.5, known = "given"),
    additive = "afvArea")
  gaussian <- do.call(ssn_lm, args)
  args$ssn.object$obs$Summer_mn <- round(args$ssn.object$obs$Summer_mn)
  glm <- do.call(ssn_glm, c(args, list(family = "poisson")))
  for (fit in list(gaussian, glm)) {
    expected <- loocv(fit, cv_predict = TRUE, se.fit = TRUE)
    set.seed(2)
    expect_equal(kcv(fit, k = fit$n, cv_predict = TRUE, se.fit = TRUE), expected)
    expect_equal(kcv(fit, folds_index = rev(seq_len(fit$n)),
                     cv_predict = TRUE, se.fit = TRUE), expected)
    expect_equal(kcv(fit, folds_index = seq_len(fit$n), cv_predict = TRUE, se.fit = TRUE,
      local = list(method = "all", parallel = TRUE, ncores = 2)), expected, tolerance = 1e-8)
    expected_local <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, local = list(size = 12))
    expect_equal(kcv(fit, folds_index = seq_len(fit$n), cv_predict = TRUE, se.fit = TRUE,
      local = list(size = 12, parallel = TRUE, ncores = 2)), expected_local, tolerance = 1e-8)
  }
  received <- NULL
  local_mocked_bindings(loocv = function(object, ...) {
    received <<- list(...)
    "delegated"
  })
  for (fit in list(gaussian, glm)) {
    expect_identical(kcv(fit, folds_index = seq_len(fit$n),
      local = list(parallel = TRUE, ncores = 2)), "delegated")
    expect_true(received$local$parallel)
    expect_identical(received$local$ncores, 2)
  }
})

test_that("local kcv agrees across serial and parallel grouped refits", {
  ssn <- mf04p
  ssn_create_bigdist(ssn, overwrite = TRUE)
  ssn$obs$cv_count <- round(ssn$obs$Summer_mn)
  ssn$obs$cv_offset <- seq_len(nrow(ssn$obs)) / 40
  ssn$obs$cv_count[c(7, 31)] <- NA
  args <- list(formula = cv_count ~ ELEV_DEM + offset(cv_offset), ssn.object = ssn,
    tailup_initial = tailup_initial("exponential", de = 2, range = 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.5, known = "given"),
    additive = "afvArea",
    local = list(index = rep(1:3, length.out = nrow(ssn$obs)), var_adjust = "theoretical"))
  fits <- list(do.call(ssn_lm, args), do.call(ssn_glm, c(args, list(family = "poisson"))))
  folds <- rep(1:2, length.out = nrow(ssn$obs))
  for (fit in fits) {
    set.seed(2)
    serial <- kcv(fit, folds_index = folds, cv_predict = TRUE, se.fit = TRUE,
                  local = list(size = 12))
    set.seed(2)
    parallel <- kcv(fit, folds_index = folds, cv_predict = TRUE, se.fit = TRUE,
                    local = list(size = 12, parallel = TRUE, ncores = 2))
    expect_equal(parallel, serial, tolerance = 1e-8)
    expect_true(all(is.finite(parallel$cv_predict)))
    expect_true(all(is.finite(parallel$se.fit)))
  }
})

# kcv
make_kcv_initials <- function() {
  list(
    tailup = tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range")),
    taildown = taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range")),
    nugget = nugget_initial("nugget", nugget = 0.5, known = "nugget")
  )
}

direct_kcv_lm <- function(Sigma, X, y, offset, folds) {
  y_krige <- y - offset
  lapply(folds, function(held) {
    train <- setdiff(seq_along(y), held)
    Sigma_inv <- chol2inv(chol(Sigma[train, train, drop = FALSE]))
    X_train <- X[train, , drop = FALSE]
    cov_beta <- chol2inv(chol(crossprod(X_train, Sigma_inv %*% X_train)))
    beta <- cov_beta %*% crossprod(X_train, Sigma_inv %*% y_krige[train])
    held_c <- Sigma[held, train, drop = FALSE]
    fit <- X[held, , drop = FALSE] %*% beta +
      held_c %*% Sigma_inv %*% (y_krige[train] - X_train %*% beta) + offset[held]
    Q <- X[held, , drop = FALSE] - held_c %*% Sigma_inv %*% X_train
    variance <- Sigma[held, held, drop = FALSE] -
      held_c %*% Sigma_inv %*% t(held_c) + Q %*% cov_beta %*% t(Q)
    list(fit = as.numeric(fit), se.fit = sqrt(diag(variance)), beta = as.numeric(beta))
  })
}

test_that("resolve_cv_auto_local() escalates only when local is omitted and n exceeds 5000, matching spmodel", {
  expect_message(
    escalated <- resolve_cv_auto_local(NULL, 5001, "kcv"),
    "rerun kcv\\(\\) with local = FALSE"
  )
  expect_true(escalated)

  expect_no_message(small <- resolve_cv_auto_local(NULL, 5000, "kcv"))
  expect_false(small)

  # an explicitly supplied local (including FALSE) is never overridden,
  # regardless of n -- only an omitted (NULL) local auto-escalates
  expect_no_message(explicit_false <- resolve_cv_auto_local(FALSE, 10000, "kcv"))
  expect_false(explicit_false)
  expect_no_message(explicit_true <- resolve_cv_auto_local(TRUE, 10, "kcv"))
  expect_true(explicit_true)
  local_spec <- list(method = "all")
  expect_no_message(explicit_list <- resolve_cv_auto_local(local_spec, 10000, "kcv"))
  expect_identical(explicit_list, local_spec)

  expect_message(
    resolve_cv_auto_local(NULL, 5001, "loocv"),
    "rerun loocv\\(\\) with local = FALSE"
  )
})

test_that("kcv()/loocv() omitting local behave exactly as local = FALSE for small samples (no auto-escalation message)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    tailup_initial = tailup_initial("exponential", de = 1, range = 1000, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.1, known = "given")
  )
  expect_true(fit$n < 5000)

  set.seed(2)
  expect_no_message(cv_omitted <- kcv(fit, k = 3))
  set.seed(2)
  expect_no_message(cv_explicit <- kcv(fit, k = 3, local = FALSE))
  expect_equal(cv_omitted, cv_explicit)

  expect_no_message(loocv_omitted <- loocv(fit))
  expect_no_message(loocv_explicit <- loocv(fit, local = FALSE))
  expect_equal(loocv_omitted, loocv_explicit)
})

test_that("Gaussian kcv maps original fold rows and matches direct held-fold GLS", {
  ssn <- mf04p
  n_original <- nrow(ssn$obs)
  ssn$obs$kcv_offset <- seq_len(n_original) / 15
  ssn$obs$kcv_response <- ssn$obs$Summer_mn + ssn$obs$kcv_offset
  ssn$obs$kcv_group <- factor(rep(1:4, length.out = n_original))
  ssn$obs$kcv_partition <- factor(rep(1:2, length.out = n_original))
  missing <- c(4, 19, 37)
  ssn$obs$kcv_response[missing] <- NA
  initial <- make_kcv_initials()
  fit <- ssn_lm(
    kcv_response ~ ELEV_DEM + offset(kcv_offset), ssn,
    tailup_initial = initial$tailup, taildown_initial = initial$taildown,
    nugget_initial = initial$nugget, additive = "afvArea",
    random = ~kcv_group,
    randcov_initial = randcov_initial(kcv_group = 1, known = "given"),
    partition_factor = ~kcv_partition
  )

  # Repeated visits share each site label. The full original-row index keeps
  # assignments for response-missing rows, which kcv() must ignore safely.
  site_fold <- paste0("site_", rep(seq_len(ceiling(n_original / 3)), each = 3)[seq_len(n_original)] %% 3)
  cv <- kcv(
    fit, folds_index = site_fold, cv_predict = TRUE, se.fit = TRUE,
    interval = "prediction", level = 0.73
  )
  expect_equal(kcv(
    fit, folds_index = site_fold, cv_predict = TRUE, se.fit = TRUE,
    interval = "prediction", level = 0.73,
    local = list(method = "all", parallel = TRUE, ncores = 2)
  ), cv, tolerance = 1e-8)
  folds <- split(seq_len(fit$n), site_fold[fit$observed_index])
  frame <- model.frame(fit)
  direct <- direct_kcv_lm(
    covmatrix(fit), model.matrix(fit), model.response(frame), model.offset(frame), folds
  )
  direct_fit <- numeric(fit$n)
  direct_se <- numeric(fit$n)
  for (i in seq_along(folds)) {
    direct_fit[folds[[i]]] <- direct[[i]]$fit
    direct_se[folds[[i]]] <- direct[[i]]$se.fit
  }
  expect_equal(unname(cv$cv_predict), direct_fit, tolerance = 1e-8)
  expect_equal(unname(cv$se.fit), direct_se, tolerance = 1e-8)
  expect_equal(cv$stats$MSPE, mean((direct_fit - model.response(frame))^2), tolerance = 1e-12)
  expect_identical(
    names(cv$stats),
    c("bias", "std.bias", "MSPE", "RMSPE", "std.MSPE", "RAV", "cor2", "cover.80", "cover.90", "cover.95", "cover.73")
  )
  expect_true(any(vapply(
    direct,
    function(value) max(abs(value$beta - coef(fit))) > 1e-7,
    logical(1)
  )))

  singleton_fold <- rep(NA_character_, n_original)
  singleton_fold[fit$observed_index] <- paste0("row_", seq_len(fit$n))
  cv_singleton <- kcv(fit, folds_index = singleton_fold, cv_predict = TRUE, se.fit = TRUE)
  cv_loocv <- loocv(fit, cv_predict = TRUE, se.fit = TRUE)
  expect_equal(cv_singleton$cv_predict, cv_loocv$cv_predict, tolerance = 1e-8)
  expect_equal(cv_singleton$se.fit, cv_loocv$se.fit, tolerance = 1e-8)

  bad_fold <- site_fold
  bad_fold[fit$observed_index[1]] <- NA_character_
  expect_error(kcv(fit, folds_index = bad_fold), "every fitted observation")
  expect_error(kcv(fit, folds_index = rep("one", n_original)), "at least two")
  expect_error(kcv(fit, k = fit$n + 1), "cannot exceed")
  expect_error(kcv(fit, k = Inf), "whole number")
})

test_that("default Gaussian k-fold statistics always use internal standard errors", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    tailup_initial = tailup_initial("exponential", de = 1, range = 1000, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.1, known = "given")
  )
  folds <- rep(1:3, length.out = nobs(fit))
  default <- kcv(fit, folds_index = folds)
  returned <- kcv(fit, folds_index = folds, se.fit = TRUE)
  expect_true(all(is.finite(unlist(default))))
  expect_equal(default, returned$stats)
})

test_that("Gaussian kcv is invariant to original observation order and rejects rank-deficient folds", {
  initial <- make_kcv_initials()
  make_fit <- function(ssn) {
    ssn_lm(
      kcv_response ~ ELEV_DEM + offset(kcv_offset), ssn,
      tailup_initial = initial$tailup, taildown_initial = initial$taildown,
      nugget_initial = initial$nugget, additive = "afvArea"
    )
  }
  ssn <- mf04p
  ssn$obs$kcv_offset <- seq_len(nrow(ssn$obs)) / 19
  ssn$obs$kcv_response <- ssn$obs$Summer_mn + ssn$obs$kcv_offset
  ssn$obs$kcv_response[c(5, 23)] <- NA
  pid <- ssn_get_netgeom(ssn$obs, "pid")$pid
  fold_by_pid <- setNames(paste0("reach_", seq_along(pid) %% 3), pid)
  fit <- make_fit(ssn)
  cv <- kcv(
    fit,
    folds_index = unname(fold_by_pid[as.character(pid)]),
    cv_predict = TRUE, se.fit = TRUE
  )

  set.seed(2)
  ssn_permuted <- ssn
  ssn_permuted$obs <- ssn_permuted$obs[sample(nrow(ssn_permuted$obs)), ]
  pid_permuted_all <- ssn_get_netgeom(ssn_permuted$obs, "pid")$pid
  fit_permuted <- make_fit(ssn_permuted)
  cv_permuted <- kcv(
    fit_permuted,
    folds_index = unname(fold_by_pid[as.character(pid_permuted_all)]),
    cv_predict = TRUE, se.fit = TRUE
  )
  pid_fit <- ssn_get_netgeom(fit$ssn.object$obs, "pid")$pid
  pid_fit_permuted <- ssn_get_netgeom(fit_permuted$ssn.object$obs, "pid")$pid
  expect_equal(
    cv$cv_predict[order(pid_fit)], cv_permuted$cv_predict[order(pid_fit_permuted)],
    tolerance = 1e-8
  )
  expect_equal(
    cv$se.fit[order(pid_fit)], cv_permuted$se.fit[order(pid_fit_permuted)],
    tolerance = 1e-8
  )

  ssn_rank <- mf04p
  ssn_rank$obs$kcv_rank_group <- factor(rep(c("A", "B"), each = ceiling(nrow(ssn_rank$obs) / 2))[seq_len(nrow(ssn_rank$obs))])
  rank_fit <- ssn_lm(
    Summer_mn ~ kcv_rank_group, ssn_rank,
    tailup_initial = initial$tailup, taildown_initial = initial$taildown,
    nugget_initial = initial$nugget, additive = "afvArea"
  )
  expect_error(
    kcv(rank_fit, folds_index = as.character(rank_fit$ssn.object$obs$kcv_rank_group)),
    "rank-deficient"
  )
})

test_that("GLM kcv uses offset-free latent updates and binomial trial alignment", {
  ssn <- mf04p
  n_original <- nrow(ssn$obs)
  ssn$obs$kcv_offset <- seq_len(n_original) / 35
  ssn$obs$kcv_trials <- rep(25, n_original)
  ssn$obs$kcv_success <- pmin(ssn$obs$kcv_trials, pmax(0, round(ssn$obs$Summer_mn)))
  ssn$obs$kcv_failure <- ssn$obs$kcv_trials - ssn$obs$kcv_success
  missing <- c(7, 31)
  ssn$obs$kcv_success[missing] <- NA
  ssn$obs$kcv_failure[missing] <- NA
  initial <- make_kcv_initials()
  fit <- ssn_glm(
    cbind(kcv_success, kcv_failure) ~ ELEV_DEM + offset(kcv_offset), ssn,
    family = "binomial", tailup_initial = initial$tailup,
    taildown_initial = initial$taildown, nugget_initial = initial$nugget,
    additive = "afvArea"
  )
  folds_original <- paste0("network_", rep(seq_len(3), length.out = n_original))
  cv <- kcv(
    fit, folds_index = folds_original, cv_predict = TRUE, se.fit = TRUE,
    type = "response", delta = TRUE
  )

  Sigma <- covmatrix(fit)
  X <- model.matrix(fit)
  y <- fit$y
  lowchol <- get_cholprods_glm(Sigma, X, y)$Sig_lowchol
  Sigma_inv <- chol2inv(lowchol)
  Sigma_inv_X <- backsolve(t(lowchol), forwardsolve(lowchol, X))
  cov_beta <- chol2inv(chol(crossprod(X, Sigma_inv_X)))
  offset <- model.offset(model.frame(fit))
  w_link <- fitted(fit, type = "link")
  w <- w_link - offset
  Ptheta <- Sigma_inv - Sigma_inv_X %*% tcrossprod(cov_beta, Sigma_inv_X)
  H <- get_D(fit$family, w_link, y, fit$size, as.vector(coef(fit, type = "dispersion"))) - Ptheta
  mHinv <- solve(-H)
  folds <- split(seq_len(fit$n), folds_original[fit$observed_index])
  direct_link <- numeric(fit$n)
  direct_se <- numeric(fit$n)
  beta_change <- logical(length(folds))
  full_beta <- cov_beta %*% crossprod(X, Sigma_inv %*% w)
  for (i in seq_along(folds)) {
    held <- folds[[i]]
    train <- setdiff(seq_len(fit$n), held)
    train_inv <- chol2inv(chol(Sigma[train, train, drop = FALSE]))
    X_train <- X[train, , drop = FALSE]
    train_cov_beta <- chol2inv(chol(crossprod(X_train, train_inv %*% X_train)))
    train_beta <- train_cov_beta %*% crossprod(X_train, train_inv %*% w[train])
    held_c <- Sigma[held, train, drop = FALSE]
    held_c_inv <- held_c %*% train_inv
    held_c_inv_X <- held_c_inv %*% X_train
    weights <- X[held, , drop = FALSE] %*% train_cov_beta %*% crossprod(X_train, train_inv) +
      held_c_inv - held_c_inv_X %*% train_cov_beta %*% crossprod(X_train, train_inv)
    direct_link[held] <- as.numeric(weights %*% matrix(w[train], ncol = 1) + offset[held])
    Q <- X[held, , drop = FALSE] - held_c_inv_X
    variance <- Sigma[held, held, drop = FALSE] - tcrossprod(held_c_inv, held_c) +
      Q %*% train_cov_beta %*% t(Q)
    train_mHinv <- mHinv[train, train, drop = FALSE] -
      mHinv[train, held, drop = FALSE] %*% solve(
        mHinv[held, held, drop = FALSE], mHinv[held, train, drop = FALSE]
      )
    direct_se[held] <- sqrt(diag(variance + weights %*% train_mHinv %*% t(weights)))
    beta_change[i] <- max(abs(train_beta - full_beta)) > 1e-7
  }
  expect_true(any(beta_change))
  expect_equal(
    cv$cv_predict,
    invlink(direct_link, fit$family, fit$size),
    tolerance = 1e-8
  )
  expect_equal(
    cv$se.fit,
    get_delta_se(direct_link, direct_se, fit$family, fit$size),
    tolerance = 1e-8
  )
  expect_true(is.finite(kcv(fit, folds_index = folds_original)$RAV))
  expect_error(kcv(fit, delta = NA), "delta must be TRUE or FALSE")
  expect_equal(kcv(
    fit, folds_index = folds_original, cv_predict = TRUE, se.fit = TRUE,
    type = "response", delta = TRUE,
    local = list(method = "all", parallel = TRUE, ncores = 2)
  ), cv, tolerance = 1e-8)
})

test_that("local kcv uses fixed-parameter fold refits through big distances", {
  ssn <- mf04p
  ssn_create_bigdist(ssn, overwrite = TRUE)
  n_original <- nrow(ssn$obs)
  ssn$obs$kcv_offset <- seq_len(n_original) / 40
  ssn$obs$kcv_response <- ssn$obs$Summer_mn + ssn$obs$kcv_offset
  ssn$obs$kcv_count <- round(ssn$obs$Summer_mn)
  initial <- make_kcv_initials()
  fit_lm <- ssn_lm(
    kcv_response ~ ELEV_DEM + offset(kcv_offset), ssn,
    tailup_initial = initial$tailup, taildown_initial = initial$taildown,
    nugget_initial = initial$nugget, additive = "afvArea"
  )
  folds <- rep(c("upstream", "downstream"), length.out = fit_lm$n)
  local <- list(method = "covariance", size = 1000)
  cv_lm <- kcv(fit_lm, folds_index = folds, cv_predict = TRUE, se.fit = TRUE, local = local)
  held <- which(folds == "upstream")

  # This independently assembles one known-parameter training fold. A single
  # local estimation group and a neighborhood larger than the training sample
  # make it the all-neighbor target of the local CV calculation.
  train_ssn <- fit_lm$ssn.object
  train_ssn$obs$kcv_response[fit_lm$observed_index[held]] <- NA
  refit_lm <- ssn_lm(
    kcv_response ~ ELEV_DEM + offset(kcv_offset), train_ssn,
    tailup_initial = initial$tailup, taildown_initial = initial$taildown,
    nugget_initial = initial$nugget, additive = "afvArea",
    local = list(method = "kmeans", size = 1000, var_adjust = "none", parallel = FALSE)
  )
  direct_lm <- predict(refit_lm, ".missing", se.fit = TRUE, interval = "none", local = local)
  direct_position <- match(as.character(fit_lm$observed_index[held]), names(direct_lm$fit))
  expect_true(all(unlist(refit_lm$is_known)))
  expect_equal(cv_lm$cv_predict[held], unname(direct_lm$fit[direct_position]), tolerance = 1e-8)
  expect_equal(cv_lm$se.fit[held], unname(direct_lm$se.fit[direct_position]), tolerance = 1e-8)

  fit_glm <- ssn_glm(
    kcv_count ~ ELEV_DEM + offset(kcv_offset), ssn, family = "poisson",
    tailup_initial = initial$tailup, taildown_initial = initial$taildown,
    nugget_initial = initial$nugget, additive = "afvArea"
  )
  cv_glm <- kcv(fit_glm, folds_index = folds, cv_predict = TRUE, se.fit = TRUE, local = local)
  train_ssn <- fit_glm$ssn.object
  train_ssn$obs$kcv_count[fit_glm$observed_index[held]] <- NA
  refit_glm <- ssn_glm(
    kcv_count ~ ELEV_DEM + offset(kcv_offset), train_ssn, family = "poisson",
    tailup_initial = initial$tailup, taildown_initial = initial$taildown,
    nugget_initial = initial$nugget,
    dispersion_initial = dispersion_initial("poisson", dispersion = 1, known = "given"),
    additive = "afvArea",
    local = list(method = "kmeans", size = 1000, var_adjust = "none", parallel = FALSE)
  )
  direct_glm <- predict(refit_glm, ".missing", type = "link", se.fit = TRUE, interval = "none", local = local)
  direct_position <- match(as.character(fit_glm$observed_index[held]), names(direct_glm$fit))
  expect_true(all(unlist(refit_glm$is_known)))
  expect_equal(cv_glm$cv_predict[held], unname(direct_glm$fit[direct_position]), tolerance = 1e-8)
  expect_equal(cv_glm$se.fit[held], unname(direct_glm$se.fit[direct_position]), tolerance = 1e-8)
  parallel_local <- c(local, list(parallel = TRUE, ncores = 2))
  expect_equal(kcv(
    fit_lm, folds_index = folds, cv_predict = TRUE, se.fit = TRUE,
    local = parallel_local
  ), cv_lm, tolerance = 1e-8)
  expect_equal(kcv(
    fit_glm, folds_index = folds, cv_predict = TRUE, se.fit = TRUE,
    local = parallel_local
  ), cv_glm, tolerance = 1e-8)
})

# loocv ddf
test_that("the automatic df threshold includes n = 500 and respects overrides", {
  expect_identical(determine_ddf(NULL, 499), "satterthwaite")
  expect_identical(determine_ddf(NULL, 500), "satterthwaite")
  expect_identical(determine_ddf(NULL, 501), "asymptotic")
  expect_identical(determine_ddf("asymptotic", 45), "asymptotic")
  expect_identical(determine_ddf("satterthwaite", 501), "satterthwaite")
  for (bad in list(NA_character_, character(), c("asymptotic", "satterthwaite"), "bad")) {
    expect_error(determine_ddf(bad, 45), "ddf must be")
  }
  local_mocked_bindings(compute_satterthwaite_fit_time = function(object) {
    list(ddf = c(intercept = 20), vcov_cov = matrix(1))
  })
  expect_equal(get_fit_ddf(list(n = 500, p = 1), NULL)$ddf, c(intercept = 20))
  expect_null(get_fit_ddf(list(n = 501, p = 1), NULL)$ddf)
})

test_that("default numerical df propagate to inference without changing fitted estimates", {
  args <- list(formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea")
  fit <- do.call(ssn_lm, args)
  asymptotic <- do.call(ssn_lm, c(args, list(ddf = "asymptotic")))
  satterthwaite_fit <- do.call(ssn_lm, c(args, list(ddf = "satterthwaite")))
  expect_true(all(is.finite(fit$ddf)))
  expect_length(fit$ddf, fit$p)
  expect_equal(coef(fit), coef(asymptotic))
  expect_equal(vcov(fit), vcov(asymptotic))
  # satterthwaite() retains its own method argument and recomputes on demand
  expect_equal(fit$ddf, satterthwaite(asymptotic, method = "numeric"))
  # summary()/confint()/anova()/emmeans() have no method/ddf argument of
  # their own (except anova()'s independent ddf argument) and always reflect
  # object$ddf, so the automatic n <= 500 default (fit) matches an explicit
  # ddf = "satterthwaite" fit (satterthwaite_fit)
  expect_equal(summary(fit)$coefficients$fixed, summary(satterthwaite_fit)$coefficients$fixed)
  expect_equal(confint(fit), confint(satterthwaite_fit))
  expect_equal(as.data.frame(anova(fit)), as.data.frame(anova(satterthwaite_fit, ddf = "satterthwaite")))
  if (requireNamespace("emmeans", quietly = TRUE)) {
    a <- as.data.frame(emmeans::emmeans(fit, ~ 1, data = mf04p$obs))
    b <- as.data.frame(emmeans::emmeans(satterthwaite_fit, ~ 1, data = mf04p$obs))
    expect_equal(a$df, b$df)
    expect_true(all(is.finite(a$df)))
  }
})

test_that("Satterthwaite failures never prevent a fitted model from being returned", {
  args <- list(formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_initial = tailup_initial("exponential", 2, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.5, known = "given"), additive = "afvArea")
  baseline <- do.call(ssn_lm, c(args, list(ddf = "asymptotic")))
  failures <- list(
    function(object) stop("failed numerical derivative"),
    function(object) {
      warning("The covariance matrix of the covariance parameters is not numerically positive definite.")
      list(ddf = rep(NA_real_, object$p), vcov_cov = NULL)
    },
    function(object) list(ddf = rep(NA_real_, object$p), vcov_cov = diag(object$p)),
    function(object) list(ddf = rep(-1, object$p), vcov_cov = diag(object$p))
  )
  for (failure in failures) local({
    local_mocked_bindings(compute_satterthwaite_fit_time = failure)
    expect_no_warning(fit <- do.call(ssn_lm, args))
    expect_s3_class(fit, "ssn_lm")
    expect_null(fit$ddf)
    expect_null(fit$vcov$cov)
    expect_equal(coef(fit), coef(baseline))
    expect_equal(summary(fit)$coefficients$fixed, summary(baseline)$coefficients$fixed)
    expect_equal(confint(fit), confint(baseline))
    expect_equal(as.data.frame(anova(fit)), as.data.frame(anova(baseline)))
  })
})

test_that("unsupported default df computation falls back without warnings", {
  args <- list(formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    euclid_initial = euclid_initial("exponential", 2, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.5, known = "given"))
  expect_no_warning(known <- do.call(ssn_lm, args))
  expect_null(known$ddf)
  expect_no_warning(local_fit <- do.call(ssn_lm, c(args,
    list(local = list(index = rep(1:3, length.out = NROW(mf04p$obs)), var_adjust = "none", parallel = FALSE)))))
  expect_null(local_fit$ddf)
  expect_true(all(is.finite(predict(local_fit, "CapeHorn"))))
})

# loocv iid
test_that("IID LOOCV is exact for all local requests with fixed covariance parameters", {
  s <- mf04p
  s$obs$x <- as.numeric(scale(s$obs$ELEV_DEM))
  s$obs$off <- sin(seq_len(nrow(s$obs))) / 3
  s$obs$z <- 1 + 0.4 * s$obs$x + s$obs$off + cos(seq_len(nrow(s$obs))) / 2
  s$obs$z[c(3, 17)] <- NA
  s$obs <- s$obs[rev(seq_len(nrow(s$obs))), ]
  for (known in c("given", "none")) {
    for (use_offset in c(FALSE, TRUE)) {
      for (method in c("ml", "reml")) {
        form <- if (use_offset) z ~ x + offset(off) else z ~ x
        fit <- ssn_lm(form, s, nugget_initial = nugget_initial("nugget", 0.2, known = known),
          estmethod = method, ddf = "asymptotic")
        frame <- model.frame(fit)
        y <- model.response(frame)
        off <- if (use_offset) model.offset(frame) else rep(0, length(y))
        X <- model.matrix(fit)
        variance <- unname(coef(fit, "nugget")[["nugget"]])
        ref <- vapply(seq_along(y), function(i) {
          Xi <- X[-i, , drop = FALSE]; xi <- X[i, , drop = FALSE]
          B <- variance * solve(crossprod(Xi))
          beta <- solve(crossprod(Xi), crossprod(Xi, (y - off)[-i]))
          c(fit = as.numeric(xi %*% beta) + off[i],
            variance = variance + as.numeric(xi %*% B %*% t(xi)))
        }, numeric(2))
        for (loc in list(FALSE, TRUE, list(size = 1),
                         list(size = 8, parallel = TRUE, ncores = 2))) {
          got <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, local = loc,
            interval = "prediction", level = 0.8)
          expect_equal(got$cv_predict, unname(ref["fit", ]), tolerance = 1e-10)
          expect_equal(got$se.fit^2, unname(ref["variance", ]), tolerance = 1e-10)
          expect_equal(got$stats$cover.8,
            mean(abs(y - ref["fit", ]) < qnorm(0.9) * sqrt(ref["variance", ])))
        }
      }
    }
  }
})

test_that("IID LOOCV avoids dense covariance and respects local-fit and floor settings", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    nugget_initial = nugget_initial("nugget", 0.2, known = "given"),
    local = list(index = rep(1:3, length.out = nrow(mf04p$obs)),
      var_adjust = "theoretical", parallel = FALSE), ddf = "asymptotic")
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  expect_equal(unname(vcov(fit)), unname(0.2 * solve(crossprod(X))), tolerance = 1e-8)
  fit$diagtol <- 0.3
  expected <- vapply(seq_along(y), function(i) {
    xi <- X[i, , drop = FALSE]; Xi <- X[-i, , drop = FALSE]
    c(fit = as.numeric(xi %*% solve(crossprod(Xi), crossprod(Xi, y[-i]))),
      variance = 0.3 * (1 + as.numeric(xi %*% solve(crossprod(Xi), t(xi)))))
  }, numeric(2))
  local_mocked_bindings(covmatrix = function(...) stop("Dense covariance requested"))
  result <- loocv(fit, local = TRUE, cv_predict = TRUE, se.fit = TRUE)
  expect_equal(result$cv_predict, unname(expected["fit", ]), tolerance = 1e-8)
  expect_equal(result$se.fit^2, unname(expected["variance", ]), tolerance = 1e-8)
  expect_equal(loocv(fit, local = list(size = 1)), result$stats)
  expect_equal(loocv(fit, local = TRUE, cv_predict = TRUE)$cv_predict, result$cv_predict)
  expect_equal(loocv(fit, local = TRUE, se.fit = TRUE)$se.fit, result$se.fit)
})

test_that("random effects and spatial components retain general local LOOCV", {
  s <- mf04p
  s$obs$grp <- factor(rep(1:4, length.out = nrow(s$obs)))
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s, random = ~grp,
    randcov_initial = randcov_initial(grp = 1, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.2, known = "given"), ddf = "asymptotic")
  expect_false(is_loocv_iid(fit))
  local <- loocv(fit, local = list(size = 8), cv_predict = TRUE, se.fit = TRUE)
  exact <- loocv(fit, local = FALSE, cv_predict = TRUE, se.fit = TRUE)
  expect_false(isTRUE(all.equal(local$cv_predict, exact$cv_predict)))
  expect_equal(loocv(fit, cv_predict = TRUE, se.fit = TRUE,
    local = list(size = 8, parallel = TRUE, ncores = 2)), local, tolerance = 1e-8)
  for (type in c("tailup", "taildown", "euclid")) {
    other <- fit
    other$random <- NULL
    class(other$coefficients$params_object[[type]]) <- paste0(type, "exponential")
    expect_false(is_loocv_iid(other))
  }
})

# loocv offset
test_that("Gaussian LOOCV matches direct offset-aware leave-one-out equations", {
  ssn <- mf04p
  ssn$obs$loocv_offset <- seq_len(nrow(ssn$obs)) / 10
  ssn$obs$Summer_mn_offset <- ssn$obs$Summer_mn + ssn$obs$loocv_offset

  tailup <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  taildown <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  nugget <- nugget_initial("nugget", nugget = 0.5, known = "nugget")
  fit <- ssn_lm(
    Summer_mn_offset ~ ELEV_DEM + offset(loocv_offset), ssn,
    tailup_initial = tailup, taildown_initial = taildown,
    nugget_initial = nugget, additive = "afvArea"
  )

  cv <- loocv(fit, cv_predict = TRUE, se.fit = TRUE)
  frame <- model.frame(fit)
  y <- model.response(frame)
  offset <- model.offset(frame)
  X <- model.matrix(fit)
  Sigma <- covmatrix(fit)
  y_krige <- y - offset

  direct <- lapply(seq_along(y), function(i) {
    retain <- -i
    Sigma_inv <- chol2inv(chol(Sigma[retain, retain, drop = FALSE]))
    X_retain <- X[retain, , drop = FALSE]
    y_retain <- y_krige[retain]
    cov_beta <- chol2inv(chol(crossprod(X_retain, Sigma_inv %*% X_retain)))
    beta <- cov_beta %*% crossprod(X_retain, Sigma_inv %*% y_retain)
    c_i <- Sigma[i, retain, drop = FALSE]
    residual <- y_retain - X_retain %*% beta
    fit_krige <- X[i, , drop = FALSE] %*% beta + c_i %*% Sigma_inv %*% residual
    Q <- X[i, , drop = FALSE] - c_i %*% Sigma_inv %*% X_retain
    variance <- Sigma[i, i] - c_i %*% Sigma_inv %*% t(c_i) + Q %*% cov_beta %*% t(Q)
    list(fit = as.numeric(fit_krige + offset[i]), se.fit = sqrt(as.numeric(variance)))
  })

  direct_fit <- vapply(direct, `[[`, numeric(1), "fit")
  direct_se <- vapply(direct, `[[`, numeric(1), "se.fit")
  expect_equal(unname(cv$cv_predict), direct_fit, tolerance = 1e-8)
  expect_equal(unname(cv$se.fit), direct_se, tolerance = 1e-8)

  errors <- y - cv$cv_predict
  expect_equal(cv$stats$std.bias, mean(errors / cv$se.fit), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(cv$stats$std.bias, mean(errors / sqrt(cv$se.fit)))))

  cv_level <- loocv(fit, interval = "prediction", level = 0.73)
  expect_identical(
    names(cv_level),
    c("bias", "std.bias", "MSPE", "RMSPE", "std.MSPE", "RAV", "cor2", "cover.80", "cover.90", "cover.95", "cover.73")
  )
  expect_equal(cv_level$cover.73, mean(abs(errors / cv$se.fit) < qnorm(1 - (1 - 0.73) / 2)))
  expect_error(loocv(fit, interval = "prediction", level = 1), "strictly between 0 and 1")
  expect_error(validate_loocv_se(c(1, 0)), "finite, positive prediction standard errors")

  cv_local_all <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, local = TRUE)
  expect_false(isTRUE(all.equal(cv_local_all$cv_predict, cv$cv_predict)))
  expect_true(all(is.finite(cv_local_all$se.fit)))

  local_size <- 20
  cv_local <- loocv(
    fit, cv_predict = TRUE, se.fit = TRUE,
    local = list(method = "covariance", size = local_size)
  )
  direct_local <- lapply(seq_along(y), function(i) {
    candidates <- seq_along(y)[-i]
    retain <- candidates[rev(order(abs(Sigma[i, candidates])))][seq_len(local_size)]
    Sigma_inv <- chol2inv(chol(Sigma[retain, retain, drop = FALSE]))
    X_retain <- X[retain, , drop = FALSE]
    y_retain <- y_krige[retain]
    cov_beta <- vcov(fit)
    beta <- coef(fit)
    c_i <- Sigma[i, retain, drop = FALSE]
    fit_krige <- X[i, , drop = FALSE] %*% beta + c_i %*% Sigma_inv %*% (y_retain - X_retain %*% beta)
    Q <- X[i, , drop = FALSE] - c_i %*% Sigma_inv %*% X_retain
    variance <- Sigma[i, i] - c_i %*% Sigma_inv %*% t(c_i) + Q %*% cov_beta %*% t(Q)
    list(fit = as.numeric(fit_krige + offset[i]), se.fit = sqrt(as.numeric(variance)))
  })
  expect_equal(unname(cv_local$cv_predict), vapply(direct_local, `[[`, numeric(1), "fit"), tolerance = 1e-8)
  expect_equal(unname(cv_local$se.fit), vapply(direct_local, `[[`, numeric(1), "se.fit"), tolerance = 1e-8)
  expect_true(all(is.finite(unlist(loocv(fit, local = list(size = 1))))))
  expect_equal(loocv(fit, cv_predict = TRUE, se.fit = TRUE,
    local = list(size = local_size, parallel = TRUE, ncores = 2)), cv_local, tolerance = 1e-8)
})

test_that("GLM LOOCV separates the latent predictor from its offset", {
  ssn <- mf04p
  ssn$obs$loocv_count <- round(ssn$obs$Summer_mn)
  ssn$obs$loocv_offset <- seq_len(nrow(ssn$obs)) / 25

  tailup <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  taildown <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  nugget <- nugget_initial("nugget", nugget = 0.5, known = "nugget")
  fit <- ssn_glm(
    loocv_count ~ ELEV_DEM + offset(loocv_offset), ssn, family = "poisson",
    tailup_initial = tailup, taildown_initial = taildown,
    nugget_initial = nugget, additive = "afvArea"
  )

  cv <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, type = "link")
  Sigma <- covmatrix(fit)
  X <- model.matrix(fit)
  y <- fit$y
  Sigma_lowchol <- get_cholprods_glm(Sigma, X, y)$Sig_lowchol
  Sigma_inv <- chol2inv(Sigma_lowchol)
  Sigma_inv_X <- backsolve(t(Sigma_lowchol), forwardsolve(Sigma_lowchol, X))
  cov_beta <- chol2inv(chol(crossprod(X, Sigma_inv_X)))
  Ptheta <- Sigma_inv - Sigma_inv_X %*% tcrossprod(cov_beta, Sigma_inv_X)
  offset <- model.offset(model.frame(fit))
  w_link <- fitted(fit, type = "link")
  w_krige <- w_link - offset
  dispersion <- as.vector(coef(fit, type = "dispersion"))
  H <- get_D(fit$family, w_link, y, fit$size, dispersion) - Ptheta
  mHinv <- solve(-H)
  wX <- cbind(w_krige, X)
  Sigma_inv_wX <- Sigma_inv %*% wX

  direct <- lapply(seq_len(fit$n), get_loocv_glm,
    Sig = Sigma, SigInv = Sigma_inv, Xmat = X,
    w = matrix(w_krige, ncol = 1), wX = wX,
    SigInv_wX = Sigma_inv_wX, mHinv = mHinv, se.fit = TRUE
  )
  direct_link <- vapply(direct, `[[`, numeric(1), "pred") + offset
  direct_se <- vapply(direct, `[[`, numeric(1), "se.fit")

  expect_equal(unname(cv$cv_predict), direct_link, tolerance = 1e-8)
  expect_equal(unname(cv$se.fit), direct_se, tolerance = 1e-8)

  cv_response <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE)
  expect_equal(unname(cv_response$se.fit), get_delta_se(direct_link, direct_se, fit$family, fit$size), tolerance = 1e-8)
  expect_error(loocv(fit, delta = NA), "delta must be TRUE or FALSE")
  expect_equal(loocv(fit, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE,
    local = list(method = "all", parallel = TRUE, ncores = 2)), cv_response, tolerance = 1e-8)

  local_size <- 20
  cv_local <- loocv(
    fit, cv_predict = TRUE, se.fit = TRUE, type = "link",
    local = list(method = "covariance", size = local_size)
  )
  expect_equal(loocv(fit, cv_predict = TRUE, se.fit = TRUE, type = "link",
    local = list(size = local_size, parallel = TRUE, ncores = 2)), cv_local, tolerance = 1e-8)
  obs <- 1
  candidates <- seq_len(fit$n)[-obs]
  retain <- candidates[rev(order(abs(Sigma[obs, candidates])))][seq_len(local_size)]
  Sigma_upchol <- chol(Sigma[retain, retain, drop = FALSE])
  Sigma_lowchol <- t(Sigma_upchol)
  Sigma_inv <- chol2inv(Sigma_upchol)
  X_retain <- X[retain, , drop = FALSE]
  Sigma_inv_X <- backsolve(t(Sigma_lowchol), forwardsolve(Sigma_lowchol, X_retain))
  cov_beta <- chol2inv(chol(crossprod(X_retain, Sigma_inv_X)))
  local_precision <- Sigma_inv - Sigma_inv_X %*% tcrossprod(cov_beta, Sigma_inv_X)
  local_curvature <- get_D(fit$family, w_link[retain], y[retain], NULL, dispersion)
  local_latent_cov <- solve(local_precision - local_curvature)
  c_obs <- Sigma[obs, retain, drop = FALSE]
  weights_beta <- tcrossprod(cov_beta, Sigma_inv_X)
  c_Sigma_inv <- c_obs %*% Sigma_inv
  c_Sigma_inv_X <- c_obs %*% Sigma_inv_X
  weights <- X[obs, , drop = FALSE] %*% weights_beta + c_Sigma_inv - c_Sigma_inv_X %*% weights_beta
  direct_local_link <- as.numeric(X[obs, , drop = FALSE] %*% coef(fit) +
    c_Sigma_inv %*% (w_krige[retain] - X_retain %*% coef(fit)) + offset[obs])
  Q <- X[obs, , drop = FALSE] - c_Sigma_inv_X
  direct_local_variance <- Sigma[obs, obs] - tcrossprod(c_Sigma_inv, c_obs) +
    Q %*% tcrossprod(vcov(fit, var_correct = FALSE), Q) + weights %*% tcrossprod(local_latent_cov, weights)
  expect_equal(unname(cv_local$cv_predict[obs]), direct_local_link, tolerance = 1e-8)
  expect_equal(unname(cv_local$se.fit[obs]), sqrt(as.numeric(direct_local_variance)), tolerance = 1e-8)
  cv_local_response <- loocv(
    fit, cv_predict = TRUE, se.fit = TRUE, type = "response", delta = TRUE,
    local = list(method = "covariance", size = local_size)
  )
  expect_equal(cv_local_response$cv_predict, invlink(cv_local$cv_predict, fit$family, fit$size), tolerance = 1e-8)
  expect_equal(cv_local_response$se.fit, get_delta_se(cv_local$cv_predict, cv_local$se.fit, fit$family, fit$size), tolerance = 1e-8)
})

test_that("Gaussian local LOOCV retains covariance structure across missing and reordered rows", {
  make_fit <- function(ssn) {
    tailup <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
    taildown <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
    nugget <- nugget_initial("nugget", nugget = 0.5, known = "nugget")
    ssn_lm(
      Summer_mn_offset ~ ELEV_DEM + offset(loocv_offset), ssn,
      tailup_initial = tailup, taildown_initial = taildown,
      nugget_initial = nugget, additive = "afvArea",
      random = ~loocv_random,
      randcov_initial = randcov_initial(loocv_random = 1, known = "given"),
      partition_factor = ~loocv_partition
    )
  }

  set.seed(2)
  ssn <- mf04p
  ssn$obs$loocv_offset <- seq_len(nrow(ssn$obs)) / 15
  ssn$obs$Summer_mn_offset <- ssn$obs$Summer_mn + ssn$obs$loocv_offset
  ssn$obs$Summer_mn_offset[c(3, 17, 36)] <- NA
  ssn$obs$loocv_random <- factor(rep(1:4, length.out = nrow(ssn$obs)))
  ssn$obs$loocv_partition <- factor(rep(1:2, length.out = nrow(ssn$obs)))

  fit <- make_fit(ssn)
  cv_exact <- loocv(fit, cv_predict = TRUE, se.fit = TRUE)
  cv_local <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, local = TRUE)
  pid <- ssn_get_netgeom(fit$ssn.object$obs, "pid")$pid
  expect_false(isTRUE(all.equal(cv_local$cv_predict, cv_exact$cv_predict)))
  expect_true(all(is.finite(cv_local$se.fit)))
  frame <- model.frame(fit)
  y <- model.response(frame)
  offset <- model.offset(frame)
  X <- model.matrix(fit)
  Sigma <- covmatrix(fit)
  obs <- 1
  retain <- seq_along(y)[-obs]
  Sigma_inv <- chol2inv(chol(Sigma[retain, retain, drop = FALSE]))
  X_retain <- X[retain, , drop = FALSE]
  y_retain <- y[retain] - offset[retain]
  cov_beta <- chol2inv(chol(crossprod(X_retain, Sigma_inv %*% X_retain)))
  beta <- cov_beta %*% crossprod(X_retain, Sigma_inv %*% y_retain)
  c_obs <- Sigma[obs, retain, drop = FALSE]
  direct_fit <- X[obs, , drop = FALSE] %*% beta +
    c_obs %*% Sigma_inv %*% (y_retain - X_retain %*% beta) + offset[obs]
  Q <- X[obs, , drop = FALSE] - c_obs %*% Sigma_inv %*% X_retain
  direct_variance <- Sigma[obs, obs] - c_obs %*% Sigma_inv %*% t(c_obs) +
    Q %*% cov_beta %*% t(Q)
  expect_equal(unname(cv_exact$cv_predict[obs]), as.numeric(direct_fit), tolerance = 1e-8)
  expect_equal(unname(cv_exact$se.fit[obs]), sqrt(as.numeric(direct_variance)), tolerance = 1e-8)

  ssn_permuted <- ssn
  ssn_permuted$obs <- ssn_permuted$obs[sample(nrow(ssn_permuted$obs)), ]
  fit_permuted <- make_fit(ssn_permuted)
  cv_permuted <- loocv(fit_permuted, cv_predict = TRUE, se.fit = TRUE, local = TRUE)
  pid_permuted <- ssn_get_netgeom(fit_permuted$ssn.object$obs, "pid")$pid
  expect_equal(
    cv_local$cv_predict[order(pid)],
    cv_permuted$cv_predict[order(pid_permuted)],
    tolerance = 1e-8
  )
  expect_equal(
    cv_local$se.fit[order(pid)],
    cv_permuted$se.fit[order(pid_permuted)],
    tolerance = 1e-8
  )
})

