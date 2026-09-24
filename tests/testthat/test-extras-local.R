skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# large prediction
block_reference <- function(object, newdata_name) {
  covariance <- covmatrix(object)
  cross_covariance <- covmatrix(object, newdata_name)
  block_covariance <- covmatrix(object, newdata_name, cov_type = "pred.pred")
  model <- SSN2:::get_newdata_model_matrix(
    object, object$ssn.object$preds[[newdata_name]]
  )
  x0 <- matrix(colMeans(model$newdata_model), nrow = 1)
  response <- model.response(model.frame(object))
  offset <- model.offset(model.frame(object))
  if (!is.null(offset)) response <- response - offset
  lowchol <- t(chol(covariance))
  siginv_x <- forwardsolve(lowchol, model.matrix(object))
  siginv_c0 <- forwardsolve(lowchol, colMeans(cross_covariance))
  residuals <- forwardsolve(lowchol, response) - siginv_x %*% coefficients(object)
  fit <- as.numeric(x0 %*% coefficients(object) + Matrix::crossprod(siginv_c0, residuals))
  h <- x0 - Matrix::crossprod(siginv_c0, siginv_x)
  variance <- as.numeric(
    mean(block_covariance) - Matrix::crossprod(siginv_c0, siginv_c0) +
      h %*% Matrix::tcrossprod(vcov(object), h)
  )
  list(fit = fit, se.fit = sqrt(pmax(variance, 0)))
}

test_that("streamed block predictions retain the exact Gaussian limit", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", taildown_type = "exponential",
    euclid_type = "exponential", nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )

  reference <- block_reference(fit, "pred1km")
  streamed <- predict(
    fit, "pred1km", block = TRUE, se.fit = TRUE,
    interval = "prediction", local = FALSE
  )
  point_average <- mean(predict(fit, "pred1km", local = FALSE))
  expect_equal(streamed$fit[1, "fit"], reference$fit, tolerance = 1e-8)
  expect_equal(unname(streamed$se.fit), reference$se.fit, tolerance = 1e-8)
  expect_equal(streamed$fit[1, "fit"], point_average, tolerance = 1e-8)

  chunked <- predict(fit, "pred1km", block = TRUE, local = list(
    method = "all", method_new = "basis", size_new = Inf,
    ordering = "none", chunk_size = 17, parallel = FALSE
  ))
  expect_equal(unname(chunked), reference$fit, tolerance = 1e-8)
})

test_that("block basis and subset controls have explicit exact limits", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", taildown_type = "exponential",
    euclid_type = "exponential", nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  small <- fit
  small$ssn.object$preds$pred1km <- fit$ssn.object$preds$pred1km[1:12, , drop = FALSE]
  exact <- predict(small, "pred1km", block = TRUE, se.fit = TRUE, local = FALSE)
  basis <- predict(small, "pred1km", block = TRUE, se.fit = TRUE, local = list(
    method = "all", method_new = "basis", size_new = 4,
    ordering = "pid", chunk_size = 3, parallel = FALSE
  ))
  subset <- predict(small, "pred1km", block = TRUE, se.fit = TRUE, local = list(
    method = "all", method_new = "subset", size_new = 4,
    ordering = "pid", chunk_size = 3, parallel = FALSE
  ))
  all_nodes <- predict(small, "pred1km", block = TRUE, se.fit = TRUE, local = list(
    method = "all", method_new = "basis", size_new = 12,
    ordering = "random", chunk_size = 3, parallel = FALSE
  ))
  expect_equal(basis$fit, exact$fit, tolerance = 1e-8)
  expect_equal(all_nodes, exact, tolerance = 1e-8)
  expect_true(is.finite(subset$fit) && is.finite(subset$se.fit))

  parallel_basis <- predict(small, "pred1km", block = TRUE, se.fit = TRUE, local = list(
    method = "all", method_new = "basis", size_new = 4,
    ordering = "pid", chunk_size = 3, parallel = TRUE, ncores = 2
  ))
  expect_equal(parallel_basis, basis, tolerance = 1e-8)
})

test_that("block covariance chunks retain random and partition terms", {
  ssn <- mf04p
  group_levels <- c("a", "b", "c")
  partition_levels <- c("inside", "outside")
  ssn$obs$block_group <- factor(rep(group_levels, length.out = NROW(ssn$obs)), levels = group_levels)
  ssn$preds$pred1km$block_group <- factor(rep(group_levels, length.out = NROW(ssn$preds$pred1km)), levels = group_levels)
  ssn$obs$block_partition <- factor(rep(partition_levels, length.out = NROW(ssn$obs)), levels = partition_levels)
  ssn$preds$pred1km$block_partition <- factor(rep(partition_levels, length.out = NROW(ssn$preds$pred1km)), levels = partition_levels)
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "none", taildown_type = "none", euclid_type = "none",
    nugget_type = "nugget", random = ~block_group,
    partition_factor = ~block_partition
  )
  small <- fit
  small$ssn.object$preds$pred1km <- fit$ssn.object$preds$pred1km[1:12, , drop = FALSE]
  expected <- covmatrix(small, "pred1km", cov_type = "pred.pred")
  actual <- SSN2:::get_block_pred_covariance(small, "pred1km", 1:4, 5:9)
  expect_equal(unname(actual), unname(expected[1:4, 5:9]), tolerance = 1e-10)
})

test_that("block prediction extends partition levels from all prediction rows", {
  ssn <- mf04p
  ssn$obs$new_partition <- factor(rep(c("a", "b"), length.out = NROW(ssn$obs)))
  ssn$preds$pred1km$new_partition <- factor(
    rep(c("c", "d"), length.out = NROW(ssn$preds$pred1km))
  )
  fit <- ssn_lm(
    Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "none", taildown_type = "none", euclid_type = "none",
    nugget_initial = nugget_initial("nugget", nugget = 1, known = "given"),
    partition_factor = ~new_partition
  )
  unpartitioned <- fit
  unpartitioned$partition_factor <- NULL
  unpartitioned$partition_xlev <- NULL
  rows <- 1:8
  columns <- 1:8
  base <- SSN2:::get_block_pred_covariance(
    unpartitioned, "pred1km", rows, columns
  )
  actual <- SSN2:::get_block_pred_covariance(
    fit, "pred1km", rows, columns
  )
  groups <- ssn$preds$pred1km$new_partition
  expected <- base * outer(groups[rows], groups[columns], `==`)

  expect_equal(unname(actual), unname(expected), tolerance = 1e-12)
  expect_true(any(expected != 0))
  expect_true(any(expected == 0))
  block <- predict(
    fit, "pred1km", block = TRUE, se.fit = TRUE,
    local = list(method = "all", chunk_size = 3, parallel = FALSE)
  )
  expect_true(all(is.finite(c(block$fit, block$se.fit))))
})

test_that("GLM block prediction errors explicitly", {
  expect_error(
    predict(structure(list(), class = "ssn_glm"), block = TRUE),
    "not supported for ssn_glm"
  )
})

# local alignment
ssn_create_bigdist(mf04p, predpts = "CapeHorn", overwrite = TRUE, no_cores = 1, verbose = FALSE)

test_that("bigdata Gaussian offset alignment holds across one/multi/singleton local groups", {
  set.seed(2)
  offset_ssn <- mf04p
  offset_ssn$obs <- offset_ssn$obs[sample(nrow(offset_ssn$obs)), ]
  offset_ssn$obs$audit_offset <- seq_len(nrow(offset_ssn$obs)) / 10
  offset_ssn$obs$Summer_mn_off <- offset_ssn$obs$Summer_mn + offset_ssn$obs$audit_offset

  form <- Summer_mn_off ~ ELEV_DEM + offset(audit_offset)
  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  td <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  ng <- nugget_initial("nugget", nugget = 0.5, known = "nugget")

  n <- nrow(offset_ssn$obs)
  y_by_pid <- setNames(offset_ssn$obs$Summer_mn_off, ssn_get_netgeom(offset_ssn$obs, "pid")$pid)

  group_specs <- list(one = rep(1, n), multi = rep(c(2, 7, 10), length.out = n), singleton = seq_len(n))
  for (gi in group_specs) {
    fit_local <- ssn_lm(form, offset_ssn,
      tailup_initial = tu, taildown_initial = td, nugget_initial = ng,
      additive = "afvArea", local = list(index = gi)
    )
    reconstructed <- fitted(fit_local) + residuals(fit_local)
    expect_equal(reconstructed, y_by_pid[names(reconstructed)], tolerance = 1e-6, ignore_attr = TRUE)
  }

  # one local group must exactly reproduce the dense (exact) fit
  fit_dense <- ssn_lm(form, offset_ssn, tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea")
  fit_local_one <- ssn_lm(form, offset_ssn,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng,
    additive = "afvArea", local = list(index = rep(1, n))
  )
  expect_equal(coef(fit_local_one), coef(fit_dense), tolerance = 1e-8)
  fl <- fitted(fit_local_one)
  fd <- fitted(fit_dense)
  expect_equal(unname(fl), unname(fd[names(fl)]), tolerance = 1e-6)
})

test_that("bigdata offset alignment holds with missing responses (nonconsecutive observed pids)", {
  set.seed(2)
  miss_ssn <- mf04p
  miss_ssn$obs <- miss_ssn$obs[sample(nrow(miss_ssn$obs)), ]
  miss_ssn$obs$audit_offset <- seq_len(nrow(miss_ssn$obs)) / 10
  miss_ssn$obs$Summer_mn_off <- miss_ssn$obs$Summer_mn + miss_ssn$obs$audit_offset
  miss_ssn$obs$Summer_mn_off[c(2, 9, 17)] <- NA

  form <- Summer_mn_off ~ ELEV_DEM + offset(audit_offset)
  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  td <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  ng <- nugget_initial("nugget", nugget = 0.5, known = "nugget")

  fit_local <- ssn_lm(form, miss_ssn,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng,
    additive = "afvArea", local = list(index = rep(c(2, 7, 10), length.out = nrow(miss_ssn$obs)))
  )
  observed <- !is.na(miss_ssn$obs$Summer_mn_off)
  y_by_pid <- setNames(miss_ssn$obs$Summer_mn_off[observed], ssn_get_netgeom(miss_ssn$obs[observed, ], "pid")$pid)
  reconstructed <- fitted(fit_local) + residuals(fit_local)
  expect_equal(reconstructed, y_by_pid[names(reconstructed)], tolerance = 1e-6, ignore_attr = TRUE)
  expect_length(predict(fit_local, ".missing"), 3)
})

test_that("bigdata offset alignment holds with a random intercept or a partition factor", {
  set.seed(2)
  s <- mf04p
  s$obs <- s$obs[sample(nrow(s$obs)), ]
  s$obs$audit_offset <- seq_len(nrow(s$obs)) / 10
  s$obs$Summer_mn_off <- s$obs$Summer_mn + s$obs$audit_offset
  s$obs$audit_group <- factor(rep(1:4, length.out = nrow(s$obs)))

  form <- Summer_mn_off ~ ELEV_DEM + offset(audit_offset)
  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  td <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  ng <- nugget_initial("nugget", nugget = 0.5, known = "nugget")
  y_by_pid <- setNames(s$obs$Summer_mn_off, ssn_get_netgeom(s$obs, "pid")$pid)

  fit_rand_one <- ssn_lm(form, s,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea",
    local = list(index = rep(1, nrow(s$obs))),
    random = ~audit_group, randcov_initial = randcov_initial(audit_group = 1, known = "given")
  )
  reconstructed <- fitted(fit_rand_one) + residuals(fit_rand_one)
  expect_equal(reconstructed, y_by_pid[names(reconstructed)], tolerance = 1e-6, ignore_attr = TRUE)

  fit_rand_multi <- ssn_lm(form, s,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea",
    local = list(index = rep(c(2, 7, 10), length.out = nrow(s$obs)), var_adjust = "none"),
    random = ~audit_group, randcov_initial = randcov_initial(audit_group = 1, known = "given")
  )
  reconstructed <- fitted(fit_rand_multi) + residuals(fit_rand_multi)
  expect_equal(reconstructed, y_by_pid[names(reconstructed)], tolerance = 1e-6, ignore_attr = TRUE)

  fit_partition <- ssn_lm(form, s,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea",
    local = list(index = rep(c(2, 7, 10), length.out = nrow(s$obs))),
    partition_factor = ~audit_group
  )
  reconstructed <- fitted(fit_partition) + residuals(fit_partition)
  expect_equal(reconstructed, y_by_pid[names(reconstructed)], tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("bigdata GLM builder aligns binomial trial totals with reordered rows", {
  set.seed(2)
  binom_ssn <- mf04p
  binom_ssn$obs <- binom_ssn$obs[sample(nrow(binom_ssn$obs)), ]
  binom_ssn$obs$success <- rep(1:3, length.out = nrow(binom_ssn$obs))
  binom_ssn$obs$failure <- seq_len(nrow(binom_ssn$obs)) + 4
  initial_glm <- get_initial_object_glm("exponential", "none", "none", "nugget", NULL, NULL, NULL, NULL, "binomial", NULL)

  for (gi in list(rep(1, nrow(binom_ssn$obs)), rep(c(2, 7, 10), length.out = nrow(binom_ssn$obs)))) {
    bd <- get_data_object_bigdata_glm(
      cbind(success, failure) ~ ELEV_DEM, binom_ssn, "binomial", "afvArea", FALSE,
      initial_glm, NULL, NULL, NULL, list(index = gi)
    )
    expect_equal(bd$size, (binom_ssn$obs$success + binom_ssn$obs$failure)[bd$order_bigdata], ignore_attr = TRUE)
  }

  fit_binom_local <- ssn_glm(cbind(success, failure) ~ ELEV_DEM, binom_ssn,
    family = "binomial", tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    local = list(index = rep(c(2, 7, 10), length.out = nrow(binom_ssn$obs)))
  )
  expect_s3_class(fit_binom_local, "ssn_glm")
  expect_true(is.finite(as.numeric(logLik(fit_binom_local))))
})

# local parity
test_that("Gaussian local predictions use each target's absolute-covariance neighbors", {
  network <- mf04p
  network$obs$x <- seq(-2, 2, length.out = nrow(network$obs))
  network$obs$off <- network$obs$x / 10
  network$obs$grp <- factor(rep(c("a", "b"), length.out = nrow(network$obs)))
  network$preds$CapeHorn <- network$preds$CapeHorn[c(1, 15, 42), ]
  network$preds$CapeHorn$x <- c(-2, 0.5, 2)
  network$preds$CapeHorn$off <- c(0.1, 0.2, -0.1)
  network$preds$CapeHorn$grp <- factor(c("a", "b", "a"))
  for (partition in list(NULL, ~grp)) {
    fit <- ssn_lm(Summer_mn ~ x + offset(off), network,
      tailup_initial = tailup_initial("exponential", 0.3, 10000, known = "given"),
      euclid_initial = euclid_initial("wave", 0.4, 4000, known = "given"),
      nugget_initial = nugget_initial("nugget", 0.2, known = "given"),
      random = ~(x | grp),
      randcov_initial = randcov_initial(`1 | grp` = 0.2, `x | grp` = 2, known = "given"),
      partition_factor = partition, additive = "afvArea", ddf = "asymptotic")
    local <- list(method = "covariance", size = 8)
    batch <- predict(fit, "CapeHorn", local = local, se.fit = TRUE, interval = "prediction")
    weights <- predict(fit, "CapeHorn", local = local, type = "weight")
    V <- as.matrix(covmatrix(fit))
    C <- as.matrix(covmatrix(fit, "CapeHorn"))
    X <- model.matrix(fit)
    X0 <- model.matrix(delete.response(terms(fit)), network$preds$CapeHorn)
    y <- model.response(model.frame(fit)) - model.offset(model.frame(fit))
    B <- vcov(fit)
    expect_true(any(C < -0.01))
    for (i in seq_len(nrow(C))) {
      keep <- order(-abs(C[i, ]), -seq_len(ncol(C)))[seq_len(local$size)]
      c0 <- C[i, keep, drop = FALSE]
      A <- solve(V[keep, keep])
      H <- X0[i, , drop = FALSE] - c0 %*% A %*% X[keep, , drop = FALSE]
      mu <- X0[i, , drop = FALSE] %*% coef(fit) +
        c0 %*% A %*% (y[keep] - X[keep, , drop = FALSE] %*% coef(fit)) +
        network$preds$CapeHorn$off[i]
      target_variance <- 0.3 + 0.4 + 0.2 + 0.2 + 2 * network$preds$CapeHorn$x[i]^2
      variance <- target_variance - c0 %*% A %*% t(c0) + H %*% B %*% t(H)
      expect_equal(unname(batch$fit[i, "fit"]), as.numeric(mu), tolerance = 1e-8)
      expect_equal(unname(batch$se.fit[i]^2), as.numeric(variance), tolerance = 1e-8)
      expect_equal(unname(weights[i, -keep]), rep(0, fit$n - length(keep)))
      expected_weight <- c0 %*% A + H %*% B %*% t(X[keep, , drop = FALSE]) %*% A
      expect_equal(unname(weights[i, keep]), as.numeric(expected_weight), tolerance = 1e-8)
      single <- fit
      single$ssn.object$preds$CapeHorn <- network$preds$CapeHorn[i, ]
      one <- predict(single, "CapeHorn", local = local, se.fit = TRUE, interval = "prediction")
      expect_equal(unname(batch$fit[i, ]), as.numeric(one$fit), tolerance = 1e-8)
      expect_equal(unname(batch$se.fit[i]), unname(one$se.fit), tolerance = 1e-8)
    }
    expect_equal(predict(fit, "CapeHorn", local = list(size = fit$n), se.fit = TRUE),
      predict(fit, "CapeHorn", local = FALSE, se.fit = TRUE), tolerance = 1e-8)
    expect_equal(predict(fit, "CapeHorn", local = c(local, list(parallel = TRUE, ncores = 2)),
      se.fit = TRUE, interval = "prediction"), batch, tolerance = 1e-8)
  }
})

test_that("Gaussian local covariance ties prefer later observations", {
  weights <- get_pred(list(c0 = c(1, -2, 2, 0), x0 = matrix(1, 1, 1)),
    se.fit = FALSE, interval = "none", formula = ~1, obdata = NULL,
    cov_matrix_val = diag(4), spatial_nugget_var = 1, randcov_params = NULL,
    cov_lowchol = NULL, Xmat = matrix(1, 4, 1), y = 1:4, offset = NULL,
    betahat = 0, cov_betahat = matrix(0.25), contrasts = NULL,
    local = list(method = "covariance", size = 1), xlevels = NULL, type = "weight")$fit
  expect_equal(which(as.numeric(weights) != 0), 3L)
})

test_that("local CV preserves fitting groups and theoretical fold covariance with missing rows", {
  network <- mf04p
  set.seed(2)
  network$obs <- network$obs[sample(nrow(network$obs)), ]
  network$obs$x <- as.numeric(scale(network$obs$ELEV_DEM))
  network$obs$off <- seq(-0.2, 0.2, length.out = nrow(network$obs))
  network$obs$y <- network$obs$Summer_mn
  network$obs$count <- round(network$obs$Summer_mn)
  network$obs$y[c(2, 9, 17)] <- NA
  network$obs$count[c(2, 9, 17)] <- NA
  groups <- rep(c(3, 7, 11), length.out = nrow(network$obs))
  ssn_create_bigdist(network, overwrite = TRUE, no_cores = 1, verbose = FALSE)
  common <- list(ssn.object = network,
    tailup_initial = tailup_initial("exponential", 0.3, 10000, known = "given"),
    euclid_initial = euclid_initial("exponential", 0.4, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.2, known = "given"), additive = "afvArea",
    local = list(index = groups, var_adjust = "none"))
  for (family in c("Gaussian", "poisson")) {
    fit <- if (family == "Gaussian") {
      do.call(ssn_lm, c(list(formula = y ~ x + offset(off), ddf = "asymptotic"), common))
    } else {
      do.call(ssn_glm, c(list(formula = count ~ x + offset(off), family = family), common))
    }
    folds <- rep(c("a", "b"), length.out = fit$n)
    held <- which(folds == "a")
    local <- list(method = "covariance", size = 8)
    refit <- if (family == "Gaussian") get_kcv_local_lm_refit(fit, held, local) else
      get_kcv_local_glm_refit(fit, held, local)
    expect_equal(unname(refit$local_index), unname(groups[fit$observed_index[-held]]))
    expect_true(all(unlist(refit$is_known)))
    if (family == "Gaussian") {
      V <- as.matrix(covmatrix(fit))[-held, -held]
      X <- model.matrix(fit)[-held, , drop = FALSE]
      W <- matrix(0, nrow(V), ncol(V))
      for (g in unique(refit$local_index)) {
        rows <- which(refit$local_index == g)
        W[rows, rows] <- solve(V[rows, rows, drop = FALSE])
      }
      naive <- solve(crossprod(X, W %*% X))
      weights <- naive %*% t(X) %*% W
      theoretical <- weights %*% V %*% t(weights)
      expect_gt(max(abs(theoretical - naive)), 1e-5)
      expect_equal(unname(vcov(refit)), unname(theoretical), tolerance = 1e-8)
    }
    cv <- kcv(fit, folds_index = folds, local = local, cv_predict = TRUE, se.fit = TRUE,
      type = "link")
    direct <- predict(refit, ".missing", local = local, se.fit = TRUE,
      type = if (family == "Gaussian") "response" else "link")
    positions <- match(as.character(fit$observed_index[held]), names(direct$fit))
    expect_equal(cv$cv_predict[held], unname(direct$fit[positions]), tolerance = 1e-8)
    expect_equal(cv$se.fit[held], unname(direct$se.fit[positions]), tolerance = 1e-8)
  }
})

# behavior parity
test_that("GLM local prediction is independent of other requested locations", {
  network <- mf04p
  network$obs$off <- seq(-0.2, 0.2, length.out = NROW(network$obs))
  network$preds$CapeHorn <- network$preds$CapeHorn[c(1, 15, 42), ]
  network$preds$CapeHorn$off <- c(0.1, 0.2, 0.3)
  fit <- ssn_glm(C16 ~ ELEV_DEM + offset(off), network, family = "poisson",
    tailup_initial = tailup_initial("exponential", 0.2, 10000, known = "given"),
    euclid_initial = euclid_initial("exponential", 0.3, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.2, known = "given"), additive = "afvArea")
  local <- list(method = "covariance", size = 12, parallel = FALSE)
  for (correct in c(FALSE, TRUE)) {
    batch <- predict(fit, "CapeHorn", local = local, se.fit = TRUE, var_correct = correct)
    for (i in seq_len(3)) {
      single <- fit
      single$ssn.object$preds$CapeHorn <- fit$ssn.object$preds$CapeHorn[i, ]
      individual <- predict(single, "CapeHorn", local = local, se.fit = TRUE, var_correct = correct)
      expect_equal(batch$fit[i], individual$fit, ignore_attr = TRUE, tolerance = 1e-8)
      expect_equal(batch$se.fit[i], individual$se.fit, ignore_attr = TRUE, tolerance = 1e-8)
    }
  }
  expect_equal(predict(fit, "CapeHorn", local = list(size = fit$n), se.fit = TRUE),
    predict(fit, "CapeHorn", local = FALSE, se.fit = TRUE), tolerance = 1e-8)
  expect_equal(dim(predict(fit, "CapeHorn", local = local, type = "weight")), c(3L, fit$n))
})

test_that("local GLM prediction matches spmodel with offsets and intervals", {
  fit_tolerance <- 0.1
  se_tolerance <- 2e-4
  beta_tolerance <- 0.03
  network <- mf04p
  network$preds$CapeHorn <- network$preds$CapeHorn[c(1, 15, 42), ]
  network$obs$off <- seq(-0.2, 0.2, length.out = NROW(network$obs))
  network$preds$CapeHorn$off <- c(0.1, 0.2, 0.3)
  observed <- sf::st_drop_geometry(network$obs)
  observed$cx <- sf::st_coordinates(network$obs)[, 1]
  observed$cy <- sf::st_coordinates(network$obs)[, 2]
  newdata <- sf::st_drop_geometry(network$preds$CapeHorn)
  newdata$cx <- sf::st_coordinates(network$preds$CapeHorn)[, 1]
  newdata$cy <- sf::st_coordinates(network$preds$CapeHorn)[, 2]
  a <- ssn_glm(C16 ~ ELEV_DEM + offset(off), network, family = "poisson",
    euclid_initial = euclid_initial("exponential", 0.3, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.2, known = "given"))
  b <- spmodel::spglm(C16 ~ ELEV_DEM + offset(off), observed, family = "poisson",
    spcov_initial = spmodel::spcov_initial("exponential", de = 0.3, ie = 0.2, range = 12000, known = "given"),
    xcoord = "cx", ycoord = "cy")
  local <- list(method = "covariance", size = 12, parallel = FALSE)
  for (correct in c(FALSE, TRUE)) for (type in c("link", "response")) {
    actual <- predict(a, "CapeHorn", local = local, type = type, interval = "prediction", se.fit = TRUE, var_correct = correct)
    expected <- predict(b, newdata, local = local, type = type, interval = "prediction", se.fit = TRUE, var_correct = correct)
    expect_equal(unname(actual$fit), unname(expected$fit), tolerance = fit_tolerance)
    expect_equal(unname(actual$se.fit), unname(expected$se.fit), tolerance = se_tolerance)
  }
  set.seed(2)
  beta_a <- conditional(a, "CapeHorn", output = "beta", samples = 200)
  set.seed(2)
  beta_b <- spmodel::conditional(b, newdata, output = "beta", samples = 200)
  expect_equal(unname(beta_a), unname(beta_b), tolerance = beta_tolerance)
})

test_that("decorrelation is invariant to common covariance scaling", {
  for (local in list(FALSE, list(method = "covariance", size = 12))) {
    transformed <- lapply(c(0.2, 3), function(scale) {
      ssn_decorrelate_data(Summer_mn ~ ELEV_DEM, mf04p, additive = "afvArea",
        tailup_params = tailup_params("exponential", scale, 10000),
        taildown_params = taildown_params("exponential", 2 * scale, 10000),
        euclid_params = euclid_params("exponential", scale, 12000),
        nugget_params = nugget_params("nugget", scale), local = local)
    })
    expect_equal(transformed[[1]]$tX, transformed[[2]]$tX, tolerance = 1e-9)
    expect_equal(transformed[[1]]$ty, transformed[[2]]$ty, tolerance = 1e-9)
    newdata <- lapply(transformed, ssn_decorrelate_newdata, "CapeHorn")
    expect_equal(newdata[[1]]$tX_newdata, newdata[[2]]$tX_newdata, tolerance = 1e-9)
    expect_equal(newdata[[1]]$yscale, newdata[[2]]$yscale, tolerance = 1e-9)
  }
})

# kmeans argument forwarding
test_that("local kmeans control values are forwarded to kmeans() by value, not by name", {
  set.seed(2)
  coords <- cbind(x = runif(24), y = runif(24))
  d_sf <- sf::st_as_sf(data.frame(x1 = coords[, 1], x2 = coords[, 2]), coords = c("x1", "x2"))

  # capture every argument actually reaching kmeans() (including x/centers/
  # iter.max) without changing its behavior; a formal signature naming
  # iter.max explicitly would swallow it out of "..." and hide a duplicate
  # or missing forward, so capture everything through "..." instead
  capture_kmeans_call <- function(local) {
    captured_args <- NULL
    real_kmeans <- kmeans
    local_mocked_bindings(kmeans = function(...) {
      captured_args <<- list(...)
      do.call(real_kmeans, captured_args)
    })
    index <- SSN2:::get_local_estimation_index(local, d_sf, 24)
    list(captured = captured_args, index = index)
  }

  # a nondefault nstart must fail before the fix (extra names, not values,
  # were previously spliced into the do.call()) and succeed after it
  local_nstart <- list(method = "kmeans", groups = 3, size = 8, nstart = 5)
  result <- capture_kmeans_call(local_nstart)
  expect_equal(result$captured$nstart, 5)
  expect_equal(sum(names(result$captured) == "nstart"), 1)

  # a nondefault algorithm is forwarded the same way
  local_alg <- list(method = "kmeans", groups = 3, size = 8, algorithm = "Lloyd")
  result_alg <- capture_kmeans_call(local_alg)
  expect_equal(result_alg$captured$algorithm, "Lloyd")

  # reserved local-control names are never forwarded as kmeans args
  local_reserved <- list(method = "kmeans", groups = 3, size = 8, parallel = FALSE, ncores = 2, var_adjust = "none")
  result_reserved <- capture_kmeans_call(local_reserved)
  expect_false(any(c("parallel", "ncores", "var_adjust") %in% names(result_reserved$captured)))

  # default (no extra controls): iter.max appears exactly once, at its default
  result_default <- capture_kmeans_call(list(method = "kmeans", groups = 3, size = 8))
  expect_equal(sum(names(result_default$captured) == "iter.max"), 1)
  expect_equal(result_default$captured$iter.max, 30)

  # an explicit iter.max overrides the default without a duplicate-argument
  # error -- this reproduces a real error (kmeans()'s own "formal argument
  # \"iter.max\" matched by multiple actual arguments") if the extra control
  # is simply appended alongside the hardcoded default instead of replacing it
  result_itermax <- capture_kmeans_call(list(method = "kmeans", groups = 3, size = 8, iter.max = 5))
  expect_equal(sum(names(result_itermax$captured) == "iter.max"), 1)
  expect_equal(result_itermax$captured$iter.max, 5)
})

test_that("forwarded kmeans controls match a direct seeded stats::kmeans() call", {
  set.seed(2)
  coords <- cbind(x = runif(30), y = runif(30))
  d_sf <- sf::st_as_sf(data.frame(x1 = coords[, 1], x2 = coords[, 2]), coords = c("x1", "x2"))
  x <- sf::st_coordinates(d_sf)

  for (local in list(
    list(method = "kmeans", groups = 4, size = 8),
    list(method = "kmeans", groups = 4, size = 8, nstart = 6),
    list(method = "kmeans", groups = 4, size = 8, nstart = 3, algorithm = "MacQueen"),
    list(method = "kmeans", groups = 4, size = 8, iter.max = 5)
  )) {
    set.seed(2)
    actual <- SSN2:::get_local_estimation_index(local, d_sf, 30)
    extra <- local[setdiff(names(local), c("size", "groups", "method", "index", "parallel", "ncores", "var_adjust"))]
    set.seed(2)
    expected <- do.call(kmeans, modifyList(list(x = x, centers = local$groups, iter.max = 30), extra))$cluster
    expect_equal(as.integer(actual), as.integer(expected))
  }
})

test_that("local$index still bypasses kmeans even when extra kmeans-like controls are present", {
  set.seed(2)
  coords <- cbind(x = runif(10), y = runif(10))
  d_sf <- sf::st_as_sf(data.frame(x1 = coords[, 1], x2 = coords[, 2]), coords = c("x1", "x2"))
  local <- list(index = rep(c(1, 2), 5), nstart = 5, method = "kmeans")
  called <- FALSE
  local_mocked_bindings(kmeans = function(...) {
    called <<- TRUE
    stop("kmeans() should not be called when local$index is supplied")
  })
  built <- SSN2:::get_local_list_estimation(local, d_sf, 10, NULL)
  expect_false(called)
  expect_equal(built$index, rep(c(1, 2), 5))
})

test_that("public Gaussian and GLM local fitting accept supplied kmeans controls", {
  ssn_create_bigdist(mf04p, predpts = "CapeHorn", overwrite = TRUE, no_cores = 1, verbose = FALSE)
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    local = list(size = 30, nstart = 3)
  )
  expect_s3_class(fit, "ssn_lm")
  expect_true(all(is.finite(coef(fit))))

  fit_glm <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p,
    family = "Gamma", tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    local = list(size = 30, nstart = 3)
  )
  expect_s3_class(fit_glm, "ssn_glm")
  expect_true(all(is.finite(coef(fit_glm))))
})

test_that("CV bias reports observed minus predicted", {
  predicted <- c(1, 2, 5)
  observed <- c(2, 4, 6)
  se <- c(1, 2, 1)
  gaussian <- get_kcv_lm_stats(predicted, observed, se, "none", 0.95)
  expect_equal(gaussian$bias, 4/3)
  expect_equal(gaussian$std.bias, 1)
  expect_equal(get_kcv_glm_stats(predicted, observed, se)$bias, 4/3)
})

