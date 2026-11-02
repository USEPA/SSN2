skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")
skip_if_not_installed("pracma")

stream_wave_network <- function(multiple = FALSE, prediction_only = FALSE) {
  network <- mf04p
  network$obs <- network$obs[rev(round(seq(1, NROW(network$obs), length.out = 12))), ]
  network$preds <- list(pred1km = network$preds$pred1km[c(5, 1, 8, 3), ])
  if (multiple || prediction_only) {
    pred <- mf04p$preds$pred1km
    rows <- unlist(lapply(split(seq_len(NROW(pred)), pred$netID), head, 2), use.names = FALSE)
    network$preds$pred1km <- pred[rev(rows), ]
    stopifnot(NROW(network$preds$pred1km) == 4)
  }
  if (prediction_only) {
    id <- names(which.max(table(mf04p$obs$netID)))
    rows <- which(as.character(mf04p$obs$netID) == id)
    network$obs <- mf04p$obs[rev(rows[round(seq(1, length(rows), length.out = 12))]), ]
  }
  network$obs$x <- seq(-1, 1, length.out = 12)
  network$obs$off <- seq(-.2, .2, length.out = 12)
  network$obs$y <- c(1, 3, 2, 4, 2, 5, 4, 2, 5, 3, 6, 4)
  network$obs$g <- factor(rep(c("a", "b"), 6))
  network$obs$part <- factor(rep(c("one", "two"), each = 6))
  network$preds$pred1km$x <- c(-.5, .2, .8, 0)
  network$preds$pred1km$off <- c(.1, -.3, .2, .4)
  network$preds$pred1km$g <- factor(c("a", "b", "new", "a"))
  network$preds$pred1km$part <- factor(c("one", "two", "one", "two"))
  network$path <- withr::local_tempdir(.local_envir = parent.frame())
  stopifnot(file.copy(file.path(mf04p$path, "binaryID.db"), network$path))
  ssn_create_distmat(network, predpts = "pred1km", among_predpts = TRUE, overwrite = TRUE)
  network
}

stream_wave_fit <- function(network, component = "both", glm = FALSE,
                            estimated = FALSE, effects = FALSE, ...) {
  known <- if (estimated) "range" else "given"
  args <- list(formula = y ~ x + offset(off), ssn.object = network, additive = "afvArea",
    tailup_initial = if (component == "down") tailup_initial("none") else
      tailup_initial("wave", .4, 2500, known = known),
    taildown_initial = if (component == "up") taildown_initial("none") else
      taildown_initial("wave", .3, 3500, known = known),
    euclid_initial = euclid_initial("none"),
    nugget_initial = nugget_initial("nugget", .25, known = "given"))
  if (effects) {
    args$euclid_initial <- euclid_initial("exponential", .2, 12000,
      rotate = .3, scale = .6, known = "given")
    args$anisotropy <- TRUE
    args$random <- ~g
    args$randcov_initial <- randcov_initial(g = .15, known = "given")
    args$partition_factor <- ~part
  }
  if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
  do.call(if (glm) ssn_glm else ssn_lm, c(args, list(...)))
}

stream_wave_prediction_reference <- function(fit, rows = 1:4, keep = seq_len(fit$n)) {
  X <- model.matrix(fit)[keep, , drop = FALSE]
  pred <- fit$ssn.object$preds$pred1km[rows, ]
  Xp <- cbind(1, pred$x)
  S <- covmatrix(fit)[keep, keep, drop = FALSE]
  C <- covmatrix(fit, "pred1km")[rows, keep, drop = FALSE]
  P <- covmatrix(fit, "pred1km", cov_type = "pred.pred")[rows, rows, drop = FALSE]
  W <- C %*% solve(S)
  H <- Xp - W %*% X
  y <- (model.response(model.frame(fit)) - model.offset(model.frame(fit)))[keep]
  list(mean = as.vector(Xp %*% coef(fit) + W %*% (y - X %*% coef(fit)) + pred$off),
       covariance = P - W %*% t(C) + H %*% vcov(fit) %*% t(H))
}

test_that("stream wave fits recover GLS coefficients and covariance under ML and REML", {
  set.seed(2)
  network <- stream_wave_network()
  for (component in c("up", "down", "both")) for (method in c("ml", "reml")) {
    for (estimated in c(FALSE, TRUE)) {
      fit <- stream_wave_fit(network, component, estimated = estimated, estmethod = method)
      X <- model.matrix(fit)
      y <- model.response(model.frame(fit)) - model.offset(model.frame(fit))
      S <- covmatrix(fit)
      B <- solve(crossprod(X, solve(S, X)))
      beta <- B %*% crossprod(X, solve(S, y))
      expect_equal(as.vector(coef(fit)), as.vector(beta), tolerance = 1e-7)
      expect_equal(unname(vcov(fit)), unname(B), tolerance = 1e-7)
      expect_true(is.finite(as.numeric(logLik(fit))))
      expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
    }
  }
})

test_that("stream wave point and block predictions retain signed covariance algebra", {
  set.seed(2)
  network <- stream_wave_network()
  for (component in c("up", "down", "both")) {
    fit <- stream_wave_fit(network, component)
    reference <- stream_wave_prediction_reference(fit)
    expect_true(any(covmatrix(fit, "pred1km") < -1e-4))
    exact <- predict(fit, "pred1km", se.fit = TRUE, local = FALSE)
    expect_equal(as.vector(exact$fit), reference$mean, tolerance = 1e-8)
    expect_equal(as.vector(exact$se.fit^2), unname(diag(reference$covariance)), tolerance = 1e-8)
    expect_equal(predict(fit, "pred1km", se.fit = TRUE,
      local = list(method = "covariance", size = fit$n)), exact, tolerance = 1e-8)
    for (size_new in c(2, 4)) {
      block <- predict(fit, "pred1km", block = TRUE, se.fit = TRUE,
        local = list(method = "covariance", size = fit$n, method_new = "subset",
                     size_new = size_new, ordering = "none", chunk_size = 2))
      X <- model.matrix(fit)
      pred <- fit$ssn.object$preds$pred1km
      x0 <- matrix(c(1, mean(pred$x)), 1)
      C <- covmatrix(fit, "pred1km")
      P <- covmatrix(fit, "pred1km", cov_type = "pred.pred")
      rows <- seq_len(size_new)
      c0 <- matrix(colMeans(C[rows, , drop = FALSE]), 1)
      W <- c0 %*% solve(covmatrix(fit))
      H <- x0 - W %*% X
      y <- model.response(model.frame(fit)) - model.offset(model.frame(fit))
      mu <- x0 %*% coef(fit) + W %*% (y - X %*% coef(fit)) + mean(pred$off)
      # Subsets approximate off-diagonal averages while retaining the full-grid diagonal.
      s0 <- mean(P[rows, rows]) - sum(diag(P)[rows]) / size_new^2 + sum(diag(P)) / 16
      variance <- s0 - W %*% t(c0) + H %*% vcov(fit) %*% t(H)
      expect_equal(as.numeric(block$fit), as.numeric(mu), tolerance = 1e-8)
      expect_equal(as.numeric(block$se.fit^2), as.numeric(variance), tolerance = 1e-8)
    }
    block <- predict(fit, "pred1km", block = TRUE, se.fit = TRUE, local = FALSE)
    expect_equal(as.numeric(block$fit), mean(reference$mean), tolerance = 1e-8)
    expect_equal(as.numeric(block$se.fit^2), mean(reference$covariance), tolerance = 1e-8)
    small <- predict(fit, "pred1km", se.fit = TRUE, local = list(size = 5))
    weights <- predict(fit, "pred1km", type = "weight", local = list(size = 5))
    C <- covmatrix(fit, "pred1km")
    for (i in 1:4) {
      keep <- order(-abs(C[i, ]), -seq_len(fit$n))[1:5]
      target <- stream_wave_prediction_reference(fit, i, keep)
      expect_equal(as.numeric(small$fit[i]), target$mean, tolerance = 1e-8)
      expect_equal(as.numeric(small$se.fit[i]^2), as.numeric(target$covariance), tolerance = 1e-8)
      expect_equal(as.vector(weights[i, -keep]), rep(0, fit$n - 5))
    }
  }
})

test_that("wave mixtures retain partitions, row order, and indexed distance agreement", {
  set.seed(2)
  network <- stream_wave_network(multiple = TRUE)
  fit <- stream_wave_fit(network, effects = TRUE)
  obs <- network$obs
  pred <- network$preds$pred1km
  S <- covmatrix(fit)
  C <- covmatrix(fit, "pred1km")
  P <- covmatrix(fit, "pred1km", cov_type = "pred.pred")
  expect_true(all(S[outer(obs$part, obs$part, "!=")] == 0))
  expect_true(all(C[outer(pred$part, obs$part, "!=")] == 0))
  expect_equal(unname(diag(S)), rep(.4 + .3 + .2 + .25 + .15, 12))
  reference <- stream_wave_prediction_reference(fit)
  point <- predict(fit, "pred1km", se.fit = TRUE, local = FALSE)
  expect_equal(as.vector(point$fit), reference$mean, tolerance = 1e-8)
  expect_equal(as.vector(point$se.fit^2), unname(diag(reference$covariance)), tolerance = 1e-8)

  ssn_create_bigdist(network, predpts = "pred1km", among_predpts = TRUE,
    overwrite = TRUE, no_cores = 1, verbose = FALSE)
  indexed <- stream_wave_fit(network, effects = TRUE, local = list(index = rep(1, 12)))
  expect_identical(SSN2:::get_decorrelate_observed_backend(indexed), "bigdata")
  expect_equal(covmatrix(indexed), S, tolerance = 1e-10)
  expect_equal(covmatrix(indexed, "pred1km"), C, tolerance = 1e-10)
  expect_equal(covmatrix(indexed, "pred1km", cov_type = "obs.pred"), t(C), tolerance = 1e-10)
  expect_equal(covmatrix(indexed, "pred1km", cov_type = "pred.pred"), P, tolerance = 1e-10)
  expect_equal(coef(indexed), coef(fit), tolerance = 1e-8)
  expect_equal(predict(indexed, "pred1km", se.fit = TRUE, local = list(size = 12)),
    point, tolerance = 1e-8)
  perm <- c(4, 2, 1, 3)
  indexed$ssn.object$preds$pred1km <- pred[perm, ]
  expect_equal(covmatrix(indexed, "pred1km"), C[perm, ], tolerance = 1e-10)
  expect_equal(covmatrix(indexed, "pred1km", cov_type = "pred.pred"), P[perm, perm], tolerance = 1e-10)
})

test_that("wave backends retain cross-network and prediction-only-network zeros", {
  set.seed(2)
  for (prediction_only in c(FALSE, TRUE)) {
    network <- stream_wave_network(multiple = TRUE, prediction_only = prediction_only)
    obs_id <- as.character(network$obs$netID)
    pred_id <- as.character(network$preds$pred1km$netID)
    expect_equal(length(unique(obs_id)), if (prediction_only) 1 else 2)
    expect_equal(length(unique(pred_id)), 2)
    references <- list()
    for (backend in c("dense", "bigdata")) {
      if (backend == "bigdata") ssn_create_bigdist(network, predpts = "pred1km",
        among_predpts = TRUE, overwrite = TRUE, no_cores = 1, verbose = FALSE)
      for (component in c("up", "down", "both")) {
        fit <- stream_wave_fit(network, component)
        S <- covmatrix(fit)
        C <- covmatrix(fit, "pred1km")
        P <- covmatrix(fit, "pred1km", cov_type = "pred.pred")
        expect_true(all(S[outer(obs_id, obs_id, "!=")] == 0))
        expect_true(all(C[outer(pred_id, obs_id, "!=")] == 0))
        expect_true(all(P[outer(pred_id, pred_id, "!=")] == 0))
        if (backend == "dense") references[[component]] <- list(S = S, C = C, P = P)
        else expect_equal(list(S = S, C = C, P = P), references[[component]], tolerance = 1e-10)
        if (prediction_only) {
          rows <- which(!pred_id %in% obs_id)
          expect_length(rows, 2)
          expect_true(all(C[rows, ] == 0))
          point <- predict(fit, "pred1km", se.fit = TRUE, local = list(size = 12))
          pred <- network$preds$pred1km[rows, ]
          Xp <- cbind(1, pred$x)
          expect_equal(as.vector(point$fit[rows]), as.vector(Xp %*% coef(fit) + pred$off), tolerance = 1e-8)
          variance <- diag(P)[rows] + diag(Xp %*% vcov(fit) %*% t(Xp))
          expect_equal(as.vector(point$se.fit[rows]^2), unname(variance), tolerance = 1e-8)
        }
      }
    }
  }
})

test_that("grouped wave fitting retains local GLS and Poisson full-block limits", {
  set.seed(2)
  network <- stream_wave_network()
  ssn_create_bigdist(network, predpts = "pred1km", among_predpts = TRUE,
    overwrite = TRUE, no_cores = 1, verbose = FALSE)
  groups <- rep(1:3, 4)
  fit <- stream_wave_fit(network, local = list(index = groups, var_adjust = "none"))
  S <- covmatrix(fit) * outer(groups, groups, "==")
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit)) - model.offset(model.frame(fit))
  B <- solve(crossprod(X, solve(S, X)))
  expect_equal(as.vector(coef(fit)), as.vector(B %*% crossprod(X, solve(S, y))), tolerance = 1e-8)
  expect_equal(unname(vcov(fit)), unname(B), tolerance = 1e-8)
  exact <- stream_wave_fit(network, glm = TRUE)
  full <- stream_wave_fit(network, glm = TRUE, local = list(index = rep(1, 12)))
  expect_equal(coef(full), coef(exact), tolerance = 1e-7)
  expect_equal(vcov(full), vcov(exact), tolerance = 1e-7)
  grouped <- stream_wave_fit(network, glm = TRUE, local = list(index = groups, var_adjust = "none"))
  expect_true(all(is.finite(c(coef(grouped), vcov(grouped)))))
  point <- predict(grouped, "pred1km", se.fit = TRUE, local = list(size = 5))
  expect_true(all(is.finite(c(point$fit, point$se.fit))))
})

test_that("Poisson wave fits support fixed and estimated variances and local prediction", {
  set.seed(2)
  network <- stream_wave_network()
  for (component in c("up", "down", "both")) for (method in c("ml", "reml")) {
    for (estimated in c(FALSE, TRUE)) {
      fit <- stream_wave_fit(network, component, glm = TRUE, estimated = estimated,
        estmethod = method)
      expect_true(all(is.finite(c(coef(fit), vcov(fit), as.numeric(logLik(fit))))))
      exact <- predict(fit, "pred1km", type = "link", se.fit = TRUE, local = FALSE)
      expect_true(all(is.finite(c(exact$fit, exact$se.fit))))
      expect_equal(predict(fit, "pred1km", type = "link", se.fit = TRUE,
        local = list(size = fit$n)), exact, tolerance = 1e-7)
      small <- predict(fit, "pred1km", type = "response", se.fit = TRUE,
        local = list(size = 5))
      expect_true(all(is.finite(c(small$fit, small$se.fit))))
      expect_true(all(small$fit > 0 & small$se.fit >= 0))
    }
  }
})

test_that("wave unconditional draws recover covariance through full local limits", {
  set.seed(2)
  network <- stream_wave_network()
  fit <- stream_wave_fit(network)
  S <- covmatrix(fit)
  args <- list(ssn.object = network, additive = "afvArea",
    tailup_params = tailup_params("wave", .4, 2500),
    taildown_params = taildown_params("wave", .3, 3500),
    euclid_params = euclid_params("none"), nugget_params = nugget_params("nugget", .25))
  for (local in list(FALSE,
    list(approximation = "vecchia", method = "covariance", size = 12, ordering = "none"),
    list(approximation = "low-rank", method_base = "all"))) {
    set.seed(2)
    draws <- do.call(ssn_rnorm, c(args, list(samples = 3000, local = local)))
    expect_equal(dim(draws), c(12L, 3000L))
    expect_true(all(is.finite(draws)))
    expect_lt(max(abs(rowMeans(draws)) / sqrt(diag(S) / 3000)), 6)
    mcse <- sqrt((outer(diag(S), diag(S)) + S^2) / 2999)
    expect_lt(max(abs(cov(t(draws)) - S) / mcse), 6)
  }
  for (local in list(list(approximation = "vecchia", size = 3, ordering = "none"),
    list(approximation = "low-rank", size_base = 5, size_new = 3,
         reorder_base = "none", kmeans_new = FALSE))) {
    set.seed(2)
    draws <- do.call(ssn_rnorm, c(args, list(samples = 7, local = local)))
    expect_equal(dim(draws), c(12L, 7L))
    expect_true(all(is.finite(draws)))
  }
})

test_that("wave CV and fixed-covariance decorrelation remain finite and agree with GLS", {
  set.seed(2)
  network <- stream_wave_network()
  for (glm in c(FALSE, TRUE)) {
    fit <- stream_wave_fit(network, glm = glm)
    for (local in list(FALSE, list(size = 5))) {
      cv <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, local = local)
      expect_true(all(is.finite(c(cv$cv_predict, cv$se.fit))))
      cv <- kcv(fit, folds_index = rep(1:3, 4), cv_predict = TRUE, se.fit = TRUE, local = local)
      expect_true(all(is.finite(c(cv$cv_predict, cv$se.fit))))
    }
  }
  args <- list(formula = y ~ x + offset(off), ssn.object = network, additive = "afvArea",
    tailup_params = tailup_params("wave", .4, 2500),
    taildown_params = taildown_params("wave", .3, 3500),
    euclid_params = euclid_params("none"), nugget_params = nugget_params("nugget", .25),
    ordering = "none")
  transformed <- do.call(ssn_decorrelate_data, c(args, list(local = FALSE)))
  S <- covmatrix(transformed$covariance_fit) / transformed$total_var
  X <- transformed$X
  y <- transformed$y
  expect_equal(unname(crossprod(transformed$tX)), unname(crossprod(X, solve(S, X))), tolerance = 1e-8)
  expect_equal(as.vector(crossprod(transformed$tX, transformed$ty)),
    as.vector(crossprod(X, solve(S, y))), tolerance = 1e-8)
  expect_equal(as.numeric(crossprod(transformed$ty)), as.numeric(crossprod(y, solve(S, y))), tolerance = 1e-8)
  expect_equal(do.call(ssn_decorrelate_data, c(args, list(local = list(size = 12))))$tX,
    transformed$tX, tolerance = 1e-8)
  pred <- ssn_decorrelate_newdata(transformed, "pred1km")
  expect_true(all(is.finite(c(pred$tX_newdata, pred$yscale, pred$yoffset))))
  for (dense in c(FALSE, TRUE)) {
    expect_warning(grid <- ssn_decorrelate_grid(y ~ x, network,
      tailup_type = "wave", taildown_type = "wave", additive = "afvArea", dense_grid = dense), NA)
    expect_true(any(grid$tailup_type == "wave" & grid$taildown_type == "wave"))
    expect_true(all(is.finite(grid$tailup_range[grid$tailup_type == "wave"])))
  }
})

test_that("installed wave workers match serial prediction and decorrelation", {
  skip_if(requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("SSN2"),
          "installed-package worker coverage")
  set.seed(2)
  network <- stream_wave_network()
  ssn_create_bigdist(network, predpts = "pred1km", among_predpts = TRUE,
    overwrite = TRUE, no_cores = 1, verbose = FALSE)
  fit <- stream_wave_fit(network)
  serial <- list(size = 5, parallel = FALSE, ncores = 2)
  parallel <- modifyList(serial, list(parallel = TRUE))
  expect_equal(predict(fit, "pred1km", se.fit = TRUE, local = parallel),
    predict(fit, "pred1km", se.fit = TRUE, local = serial), tolerance = 1e-10)
  for (block in c(FALSE, TRUE)) {
    expect_equal(predict(fit, "pred1km", block = block, se.fit = TRUE,
      local = c(parallel, list(size_new = 2, ordering = "none", chunk_size = 2))),
      predict(fit, "pred1km", block = block, se.fit = TRUE,
        local = c(serial, list(size_new = 2, ordering = "none", chunk_size = 2))), tolerance = 1e-10)
  }
  args <- list(formula = y ~ x, ssn.object = network, additive = "afvArea",
    tailup_params = tailup_params("wave", .4, 2500),
    taildown_params = taildown_params("wave", .3, 3500),
    nugget_params = nugget_params("nugget", .25), ordering = "none")
  a <- do.call(ssn_decorrelate_data, c(args, list(local = serial)))
  b <- do.call(ssn_decorrelate_data, c(args, list(local = parallel)))
  expect_equal(a$tX, b$tX, tolerance = 1e-10)
  expect_equal(a$ty, b$ty, tolerance = 1e-10)
  expect_equal(ssn_decorrelate_newdata(a, "pred1km")$tX_newdata,
    ssn_decorrelate_newdata(b, "pred1km")$tX_newdata, tolerance = 1e-10)
})
