skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_lowrank_fixture <- function(glm = FALSE, effects = FALSE, stream = FALSE, missing = FALSE) {
  network <- mf04p
  network$obs <- network$obs[round(seq(1, NROW(network$obs), length.out = 12)), ]
  network$preds <- list(pred1km = network$preds$pred1km[c(5, 1, 8, 3), ])
  network$obs$x <- seq(-1, 1, length.out = 12)
  network$obs$off <- seq(-.2, .2, length.out = 12)
  network$obs$y <- c(1, 3, 2, 4, 2, 5, 4, 2, 5, 3, 6, 4)
  network$obs$g <- factor(rep(c("a", "b", "c"), each = 4))
  network$obs$part <- factor(rep(c("one", "two"), each = 6))
  network$preds$pred1km$x <- c(-.5, .2, .8, 0)
  network$preds$pred1km$off <- c(.1, -.3, .2, .4)
  network$preds$pred1km$g <- factor(c("a", "c", "new", "new"))
  network$preds$pred1km$part <- factor(c("one", "two", "one", "one"))
  if (stream) {
    root <- tempfile("conditional-lowrank-")
    dir.create(root)
    stopifnot(file.copy(network$path, root, recursive = TRUE, copy.mode = FALSE))
    network$path <- file.path(root, basename(network$path))
    ssn_create_distmat(network, predpts = "pred1km", among_predpts = TRUE, overwrite = TRUE)
  }
  if (missing) network$obs$y[c(2, 7)] <- NA
  args <- list(formula = y ~ x + offset(off), ssn.object = network, additive = "afvArea",
    tailup_initial = if (stream) tailup_initial("exponential", .3, 10000, known = "given") else tailup_initial("none"),
    taildown_initial = if (stream) taildown_initial("spherical", .2, 10000, known = "given") else taildown_initial("none"),
    euclid_initial = euclid_initial("exponential", .4, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", .25, known = "given"))
  if (effects) {
    args$random <- ~g
    args$randcov_initial <- randcov_initial(g = .3, known = "given")
    args$partition_factor <- ~part
  }
  if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
  do.call(if (glm) ssn_glm else ssn_lm, args)
}

conditional_lowrank_settings <- function(size = 6, block = 2, chunk = 1000) {
  list(approximation = "low-rank", size_base = size, reorder_base = "none",
       size_new = block, reorder_new = "none", kmeans_new = FALSE, chunk_size = chunk)
}

conditional_lowrank_oracle <- function(fit, base, blocks) {
  glm <- inherits(fit, "ssn_glm")
  X <- model.matrix(fit)[base, , drop = FALSE]
  Xp <- cbind(1, fit$ssn.object$preds$pred1km$x)
  Sigma <- as.matrix(covmatrix(fit))[base, base]
  spatial <- conditional_spatial_object(fit)
  S <- as.matrix(covmatrix(spatial))[base, base]
  Sp <- as.matrix(covmatrix(spatial, "pred1km", cov_type = "pred.pred"))
  Spo <- t(as.matrix(covmatrix(spatial, "pred1km", cov_type = "obs.pred"))[base, , drop = FALSE])
  W <- Spo %*% solve(S)
  R <- matrix(0, NROW(Xp), NROW(Xp))
  for (rows in blocks) R[rows, rows] <- Sp[rows, rows] - W[rows, , drop = FALSE] %*% t(Spo[rows, , drop = FALSE])
  if (!is.null(fit$random)) {
    design <- get_conditional_random_design(fit, fit$ssn.object$preds$pred1km, base)
    Z <- as.matrix(design$observed)
    Zp <- as.matrix(design$newdata)
    H <- solve(diag(1 / design$variance) + t(Z) %*% solve(S, Z))
    A <- Zp - W %*% Z
    R <- R + A %*% H %*% t(A)
    W <- W + A %*% H %*% t(Z) %*% solve(S)
  }
  C <- as.matrix(vcov(fit))
  eta <- if (glm) fitted(fit, type = "link") else model.response(model.frame(fit))
  latent_mean <- if (glm) eta + fit$laplace$delta[names(eta)] else eta
  center <- (latent_mean - model.offset(model.frame(fit)))[base]
  V <- if (glm) solve(solve(Sigma) + diag(exp(eta[base]))) else matrix(0, length(base), length(base))
  M <- if (glm) V %*% solve(Sigma, X) else matrix(0, length(base), NCOL(X))
  B <- Xp + W %*% (M - X)
  beta <- as.vector(if (glm) fit$laplace$betahat else coef(fit))
  mean <- as.vector(Xp %*% beta + W %*% (center - X %*% beta) + fit$ssn.object$preds$pred1km$off)
  covariance <- R + W %*% V %*% t(W) + B %*% C %*% t(B)
  list(mean = c(beta, mean), covariance = rbind(cbind(C, C %*% t(B)), cbind(B %*% C, covariance)),
       V = V, M = M, Sigma = Sigma)
}

expect_lowrank_moments <- function(draws, mean, covariance) {
  n <- NCOL(draws)
  expect_lt(max(abs(rowMeans(draws) - mean) / sqrt(diag(covariance) / n)), 6)
  mcse <- sqrt((outer(diag(covariance), diag(covariance)) + covariance^2) / (n - 1))
  expect_lt(max(abs(cov(t(draws)) - covariance) / mcse), 6)
}

test_that("reduced low-rank draws recover base and shared-effect block moments", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) for (effects in c(FALSE, TRUE)) {
    fit <- conditional_lowrank_fixture(glm, effects, stream = TRUE)
    reference <- conditional_lowrank_oracle(fit, 1:6, list(1:2, 3:4))
    for (chunk in c(1000, 333)) {
      set.seed(2)
      value <- conditional(fit, "pred1km", output = c("beta", "newdata"),
                           samples = 6000, local = conditional_lowrank_settings(chunk = chunk))
      expect_lowrank_moments(rbind(value$beta, value$newdata), reference$mean, reference$covariance)
      h <- c(0, 0, rep(.25, 4))
      expect_lowrank_moments(matrix(colMeans(value$newdata), 1), sum(h * reference$mean),
                            matrix(as.numeric(t(h) %*% reference$covariance %*% h), 1))
    }
    if (effects) {
      design <- get_conditional_random_design(fit, fit$ssn.object$preds$pred1km, 1:6)
      expect_true(any(design$info$observed & !design$info$conditioned))
      expect_true(any(!design$info$observed))
    }
  }
})

test_that("whitened GLM base factors agree with natural-scale independent operators", {
  fit <- conditional_lowrank_fixture(TRUE, TRUE)
  index <- c(6, 1, 8, 3, 10, 12)
  oracle <- conditional_lowrank_oracle(fit, index, list(1:4))
  base <- get_conditional_glm_base(fit, index)
  L <- base$cov_lowchol_base
  Jinv <- solve(t(base$cov_lowchol_mH))
  expect_equal(tcrossprod(L), oracle$Sigma, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(tcrossprod(L %*% Jinv), oracle$V, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(L %*% base$wts_latent, oracle$M, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(tcrossprod(base$cov_betahat_lowchol), vcov(fit), ignore_attr = TRUE, tolerance = 1e-10)
  local_mocked_bindings(get_D = function(...) diag(1e6, length(index)))
  expect_error(get_conditional_glm_base(fit, index), "simulation base latent precision is not positive definite")
})

test_that("full-base single-block low-rank dispatch is exact", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) for (effects in c(FALSE, TRUE)) {
    fit <- conditional_lowrank_fixture(glm, effects)
    output <- c("all", if (glm) "latent", if (effects) "random")
    set.seed(2)
    exact <- conditional(fit, "pred1km", output = output, samples = 11, local = FALSE)
    set.seed(2)
    full <- conditional(fit, "pred1km", output = output, samples = 11,
                        local = conditional_lowrank_settings(size = 12, block = 4))
    expect_identical(full, exact)
  }
})

test_that("full-base block GLM simulation exposes latent draws with offsets once", {
  set.seed(2)
  fit <- conditional_lowrank_fixture(TRUE, TRUE)
  settings <- conditional_lowrank_settings(size = 12)
  set.seed(2)
  link <- conditional(fit, "pred1km", output = c("all", "latent"), samples = 17, local = settings)
  set.seed(2)
  response <- conditional(fit, "pred1km", output = c("all", "latent"),
                          samples = 17, local = settings, type = "response")
  expect_equal(response$newdata, exp(link$newdata))
  expect_identical(response$latent, link$latent)
  expect_identical(rownames(link$latent), rownames(model.matrix(fit)))
  expect_equal(dim(link$latent), c(12, 17))
  expect_equal(link$object, matrix(fitted(fit, type = "link"), 12, 17), ignore_attr = TRUE)
  expect_error(conditional(fit, "pred1km", output = c("newdata", "random"), samples = 1, local = settings),
               "Random-effect output is unavailable on low-rank paths", fixed = TRUE)
  set.seed(2)
  exact <- conditional(fit, output = "random", samples = 3)
  set.seed(2)
  expect_identical(conditional(fit, output = "random", samples = 3, local = settings), exact)
  expect_error(conditional(fit, "pred1km", output = c("newdata", "latent"), samples = 1,
                           local = conditional_lowrank_settings()), "Full latent output is unavailable")
  expect_error(conditional(fit, "pred1km", samples = 1,
                           local = conditional_lowrank_settings(chunk = Inf)), "positive integer")
})

test_that("serial and parallel block preparation uses the same parent RNG stream", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_lowrank_fixture(glm, TRUE)
    settings <- conditional_lowrank_settings(chunk = 3)
    set.seed(2)
    serial <- conditional(fit, "pred1km", output = "all", samples = 7, local = settings)
    set.seed(2)
    parallel <- conditional(fit, "pred1km", output = "all", samples = 7,
                            local = c(settings, list(parallel = TRUE, ncores = 2)))
    expect_identical(parallel, serial)
  }
})

test_that("low-rank preparation reads only the base and factors each block once", {
  set.seed(2)
  fit <- conditional_lowrank_fixture(TRUE, TRUE, stream = TRUE)
  ssn_create_bigdist(fit$ssn.object, predpts = "pred1km", among_predpts = TRUE,
                    overwrite = TRUE, no_cores = 1, verbose = FALSE)
  original_cross <- get_dist_pred_object
  original_block <- factor_conditional_glm_block
  original_base <- get_conditional_glm_base
  original_D <- get_D
  selected <- list()
  blocks <- bases <- curvature <- 0L
  local_mocked_bindings(
    covmatrix = function(...) stop("Unexpected full covariance"),
    get_conditional_glm_joint = function(...) stop("Unexpected full latent preparation"),
    get_block_obs_covariance = function(...) stop("Unexpected all-observation cross covariance"),
    get_dist_pred_object = function(object, ...) {
      expect_lte(NROW(object$ssn.object$obs), 6)
      selected[[length(selected) + 1L]] <<- get_decorrelate_rows(object$ssn.object$obs)$pid
      original_cross(object, ...)
    },
    get_conditional_glm_base = function(...) { bases <<- bases + 1L; original_base(...) },
    get_D = function(family, w, y, ...) {
      curvature <<- curvature + 1L
      expect_equal(length(w), 6)
      expect_identical(y, fit$y[1:6])
      original_D(family, w, y, ...)
    },
    factor_conditional_glm_block = function(...) { blocks <<- blocks + 1L; original_block(...) })
  set.seed(2)
  value <- conditional(fit, "pred1km", samples = 7, local = conditional_lowrank_settings(chunk = 2))
  expect_equal(dim(value), c(4, 7))
  expect_true(length(selected) > 0)
  expect_true(all(vapply(selected, identical, logical(1), get_decorrelate_rows(fit$ssn.object$obs)$pid[1:6])))
  expect_equal(c(bases, curvature, blocks), c(1, 1, 2))
})

test_that("low-rank missing-set predictions retain rows and beta-only calls need no distances", {
  set.seed(2)
  fit <- conditional_lowrank_fixture(TRUE, TRUE, missing = TRUE)
  settings <- conditional_lowrank_settings(size = 5, block = 1, chunk = 2)
  set.seed(2)
  missing <- conditional(fit, samples = 1, local = settings)
  expect_equal(dim(missing), c(2, 1))
  expect_identical(rownames(missing), as.character(fit$missing_index))
  local_mocked_bindings(get_conditional_lowrank_covariance = function(...) stop("Unexpected covariance"),
    get_conditional_glm_base = function(...) stop("Unexpected latent preparation"))
  fit$ssn.object$preds <- list()
  expect_equal(dim(conditional(fit, output = "beta", samples = 1, local = settings)), c(2, 1))
  expect_equal(dim(conditional(fit, output = "object", samples = 1, local = settings)), c(10, 1))
})

test_that("one-row bases and prediction blocks retain matrix dimensions", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) {
    fit <- conditional_lowrank_fixture(glm, TRUE)
    settings <- conditional_lowrank_settings(size = 1, block = 1, chunk = 1)
    set.seed(2)
    value <- conditional(fit, "pred1km", output = "all", samples = 1, local = settings)
    expect_equal(dim(value$beta), c(2, 1))
    expect_equal(dim(value$newdata), c(4, 1))
    expect_true(all(is.finite(value$newdata)))
  }
})

test_that("low-rank triangular solves do not accept singular covariance roots", {
  fit <- conditional_lowrank_fixture(TRUE)
  local_mocked_bindings(conditional_vecchia_covariance = function(context, rows, columns = rows) {
    matrix(1, length(rows), length(columns))
  })
  expect_error(get_conditional_glm_base(fit, 1:6), "not positive")
  expect_error(factor_conditional_glm_block(list(cross = matrix(1, 1, 2),
    new = matrix(1, 2, 2)), matrix(1)), "not positive")
})
