skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_vecchia_fixture <- function(glm = FALSE, effects = FALSE, stream = FALSE, missing = FALSE) {
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
    root <- tempfile("conditional-vecchia-")
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

conditional_vecchia_oracle <- function(fit, size, observed_order = seq_len(fit$n), prediction_order = 1:4) {
  X <- model.matrix(fit)
  n <- NROW(X)
  Xp <- cbind(1, fit$ssn.object$preds$pred1km$x)
  Sigma <- as.matrix(covmatrix(fit))
  spatial <- conditional_spatial_object(fit)
  S <- as.matrix(covmatrix(spatial))
  cross <- as.matrix(covmatrix(spatial, "pred1km", cov_type = "obs.pred"))
  pred <- as.matrix(covmatrix(spatial, "pred1km", cov_type = "pred.pred"))
  full <- rbind(cbind(S, cross), cbind(t(cross), pred))
  part <- if (is.null(fit$partition_factor)) rep(1, n + 4) else
    c(as.character(fit$ssn.object$obs$part), as.character(fit$ssn.object$preds$pred1km$part))
  select <- function(i, pool, count, covariance, tie = seq_len(NROW(covariance))) {
    pool <- pool[pool != i & part[pool] == part[i]]
    head(pool[order(-abs(covariance[i, pool]), tie[pool])], count)
  }
  glm <- inherits(fit, "ssn_glm")
  eta <- if (glm) fitted(fit, type = "link") else model.response(model.frame(fit))
  latent_mean <- if (glm) eta + fit$laplace$delta[names(eta)] else eta
  w <- as.vector(latent_mean - model.offset(model.frame(fit)))
  T <- matrix(0, n, NCOL(X))
  F <- matrix(0, n, n)
  operators <- list()
  for (k in seq_len(n)) {
    i <- observed_order[k]
    N <- select(i, observed_order[seq_len(k - 1L)], size, Sigma)
    K <- c(i, select(i, seq_len(n), size - 1L, Sigma))
    A <- unique(c(i, N, K))
    if (glm) {
      P <- solve(Sigma[A, A, drop = FALSE])
      V <- solve(P + diag(exp(eta[A]), length(A)))
      M <- V %*% P %*% X[A, , drop = FALSE]
      j <- match(N, A)
      g <- if (length(N)) as.vector(solve(V[j, j, drop = FALSE], V[j, 1])) else numeric()
      variance <- V[1, 1] - sum(g * V[j, 1])
      coefficient <- as.vector(M[1, ] - g %*% M[j, , drop = FALSE])
      T[i, ] <- coefficient + g %*% T[N, , drop = FALSE]
      F[i, ] <- g %*% F[N, , drop = FALSE]
      F[i, i] <- sqrt(variance)
      operators[[k]] <- list(N = N, K = K, A = A, weights = g, variance = variance,
        coefficient = coefficient, intercept = w[i] - sum(g * w[N]))
    }
  }
  W <- rbind(diag(n), matrix(0, 4, n))
  E <- matrix(0, n + 4, 4)
  for (k in 1:4) {
    i <- n + prediction_order[k]
    N <- select(i, c(seq_len(n), n + prediction_order[seq_len(k - 1L)]), size, full,
      c(seq_len(n), n + match(1:4, prediction_order)))
    g <- if (length(N)) as.vector(solve(full[N, N, drop = FALSE], full[N, i])) else numeric()
    W[i, ] <- g %*% W[N, , drop = FALSE]
    E[i, ] <- g %*% E[N, , drop = FALSE]
    E[i, prediction_order[k]] <- sqrt(full[i, i] - sum(g * full[N, i]))
  }
  list(T = T, V = tcrossprod(F), W = W[n + 1:4, , drop = FALSE],
    R = tcrossprod(E[n + 1:4, , drop = FALSE]), center = as.vector(w), X = X, Xp = Xp,
    S = S, Sigma = Sigma, operators = operators)
}

conditional_vecchia_moments <- function(fit, oracle, blocks = list(seq_len(fit$n))) {
  W <- oracle$W
  R <- oracle$R
  if (!is.null(fit$random)) {
    design <- get_conditional_random_design(fit, fit$ssn.object$preds$pred1km)
    Z <- as.matrix(design$observed)
    Zp <- as.matrix(design$newdata)
    P <- matrix(0, fit$n, fit$n)
    for (b in blocks) P[b, b] <- solve(oracle$S[b, b, drop = FALSE])
    H <- solve(diag(1 / design$variance) + t(Z) %*% P %*% Z)
    U <- H %*% t(Z) %*% P
    A <- Zp - W %*% Z
    R <- R + A %*% H %*% t(A)
    W <- W + A %*% U
  }
  B <- oracle$Xp + W %*% (oracle$T - oracle$X)
  C <- as.matrix(vcov(fit))
  beta <- if (inherits(fit, "ssn_glm")) fit$laplace$betahat else coef(fit)
  mean <- as.vector(oracle$Xp %*% beta + W %*% (oracle$center - oracle$X %*% beta) +
    fit$ssn.object$preds$pred1km$off)
  covariance <- R + W %*% oracle$V %*% t(W) + B %*% C %*% t(B)
  list(mean = c(as.vector(beta), mean), covariance = rbind(cbind(C, C %*% t(B)),
    cbind(B %*% C, covariance)))
}

expect_vecchia_moments <- function(draws, mean, covariance) {
  count <- NCOL(draws)
  expect_lt(max(abs(rowMeans(draws) - mean) / sqrt(diag(covariance) / count)), 6)
  mcse <- sqrt((outer(diag(covariance), diag(covariance)) + covariance^2) / (count - 1))
  expect_lt(max(abs(cov(t(draws)) - covariance) / mcse), 6)
}

test_that("observed Vecchia operators match independent two-neighborhood conditionals", {
  fit <- conditional_vecchia_fixture(TRUE, TRUE, stream = TRUE)
  settings <- get_local_list_conditional(list(size = 3, ordering = "none"), fit, fit$ssn.object$preds$pred1km)
  reference <- conditional_vecchia_oracle(fit, 3)
  context <- get_conditional_vecchia_glm_context(fit, "pred1km", settings)
  expect_identical(context$order, seq_len(12))
  for (k in 1:12) {
    op <- context$operators[[k]]
    target <- reference$operators[[k]]
    expect_equal(sort(op$neighbors), sort(target$N))
    expect_equal(op$weights[match(target$N, op$neighbors)], target$weights, tolerance = 1e-10)
    expect_equal(op$variance, target$variance, tolerance = 1e-10)
    expect_equal(op$intercept, target$intercept, tolerance = 1e-10)
    expect_equal(op$coefficient, target$coefficient, tolerance = 1e-10)
  }
  expect_true(any(vapply(seq_along(reference$operators), function(i) any(reference$operators[[i]]$K > i), logical(1))))
})

test_that("truncated Vecchia draws recover joint prediction and shared-effect moments", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) for (effects in c(FALSE, TRUE)) {
    fit <- conditional_vecchia_fixture(glm, effects, stream = TRUE)
    settings <- list(size = 3, ordering = "none", chunk_size = 333)
    if (effects) settings <- c(settings, list(groups_random = 3, method_random = "random"))
    set.seed(2)
    resolved <- get_local_list_conditional(settings, fit, fit$ssn.object$preds$pred1km)
    oracle <- conditional_vecchia_oracle(fit, 3)
    moments <- conditional_vecchia_moments(fit, oracle, resolved$random_blocks)
    set.seed(2)
    value <- conditional(fit, "pred1km", samples = 6000, local = settings,
      output = c("all", if (glm) "latent", if (effects) "random"))
    expect_vecchia_moments(rbind(value$beta, value$newdata), moments$mean, moments$covariance)
    h <- c(0, 0, rep(.25, 4))
    expect_vecchia_moments(matrix(colMeans(value$newdata), 1), sum(h * moments$mean),
      matrix(as.numeric(t(h) %*% moments$covariance %*% h), 1))
    if (glm) {
      C <- as.matrix(vcov(fit))
      V <- oracle$V + oracle$T %*% C %*% t(oracle$T)
      expect_vecchia_moments(rbind(value$beta, value$latent),
        c(fit$laplace$betahat, oracle$center + model.offset(model.frame(fit))),
        rbind(cbind(C, C %*% t(oracle$T)), cbind(oracle$T %*% C, V)))
    }
    if (effects) expect_true(!is.null(attr(value$random, "effect_info")))
  }
})

test_that("full-neighborhood operators recover exact beta-first covariance", {
  set.seed(2)
  fit <- conditional_vecchia_fixture(TRUE, TRUE)
  oracle <- conditional_vecchia_oracle(fit, 30)
  P <- solve(oracle$Sigma)
  V <- solve(P + diag(exp(fitted(fit, type = "link"))))
  M <- V %*% P %*% model.matrix(fit)
  expect_equal(oracle$V, V, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(oracle$T, M, ignore_attr = TRUE, tolerance = 1e-10)
  moments <- conditional_vecchia_moments(fit, oracle)
  set.seed(2)
  value <- conditional(fit, "pred1km", output = c("beta", "newdata"), samples = 6000,
    local = list(size = 30, ordering = "none", groups_random = 1))
  expect_vecchia_moments(rbind(value$beta, value$newdata), moments$mean, moments$covariance)
  set.seed(2)
  exact <- conditional(fit, "pred1km", output = c("all", "latent", "random"), samples = 7, local = FALSE)
  set.seed(2)
  expect_identical(conditional(fit, "pred1km", output = c("all", "latent", "random"), samples = 7,
    local = list(method = "all")), exact)
  expect_error(conditional(fit, "pred1km", local = list(method = "all", size_random = 4)), "only used for Vecchia")
})

test_that("Vecchia caches contexts and operators across chunks with indexed bmat reads", {
  set.seed(2)
  fit <- conditional_vecchia_fixture(TRUE, TRUE, stream = TRUE)
  ssn_create_bigdist(fit$ssn.object, predpts = "pred1km", among_predpts = TRUE,
    overwrite = TRUE, no_cores = 1, verbose = FALSE)
  original_context <- get_conditional_vecchia_covariance
  original_design <- get_conditional_random_design
  original_scores <- conditional_vecchia_scores
  original_gaussian <- prepare_conditional_local_gaussian
  original_D <- get_D
  original_cross <- get_dist_pred_object
  original_square <- get_dist_object_bigdata_cross
  contexts <- designs <- curvature <- gaussian <- 0L
  scores <- integer()
  stage <- "effects"
  read_shapes <- list()
  local_mocked_bindings(
    covmatrix = function(...) stop("Unexpected full covariance"),
    get_conditional_glm_joint = function(...) stop("Unexpected full latent factor"),
    get_conditional_vecchia_covariance = function(...) {
      contexts <<- contexts + 1L
      value <- original_context(...)
      expect_true(all(unlist(value$backends) == "bigdata"))
      value
    },
    get_conditional_random_design = function(...) { designs <<- designs + 1L; original_design(...) },
    get_D = function(...) { curvature <<- curvature + 1L; original_D(...) },
    conditional_vecchia_scores = function(context, i, method, columns = seq_len(context$n + NROW(context$newdata))) {
      scores <<- c(scores, i)
      stage <<- if (i <= context$n) "observed" else "prediction"
      if (i <= context$n) expect_identical(columns, seq_len(fit$n))
      else expect_true(all(columns <= i - 1L))
      original_scores(context, i, method, columns)
    },
    prepare_conditional_local_gaussian = function(Sigma, ...) {
      gaussian <<- gaussian + 1L
      expect_lte(NROW(Sigma), 6)
      original_gaussian(Sigma, ...)
    },
    get_dist_pred_object = function(...) {
      expect_false(stage == "observed")
      original_cross(...)
    },
    get_dist_object_bigdata_cross = function(d1, d2, ..., backend = "bigdata") {
      expect_identical(backend, "bigdata")
      read_shapes[[length(read_shapes) + 1L]] <<- c(NROW(d1), NROW(d2))
      original_square(d1, d2, ..., backend = backend)
    })
  set.seed(2)
  value <- conditional(fit, "pred1km", samples = 7,
    local = list(size = 3, ordering = "none", size_random = 4, method_random = "random", chunk_size = 2))
  expect_equal(dim(value), c(4, 7))
  expect_equal(c(contexts, designs, curvature, gaussian), c(1, 1, 1, 12))
  expect_equal(scores, seq_len(16))
  expect_true(all(vapply(read_shapes, function(x) max(x) <= 6, logical(1))))
})

test_that("Vecchia outputs retain identities and avoid unrequested prediction reads", {
  set.seed(2)
  fit <- conditional_vecchia_fixture(TRUE, TRUE, missing = TRUE)
  settings <- list(size = 3, ordering = "pid", chunk_size = 2)
  set.seed(2)
  value <- conditional(fit, output = c("all", "latent", "random"), samples = 5, local = settings)
  expect_identical(rownames(value$newdata), as.character(fit$missing_index))
  expect_identical(rownames(value$latent), rownames(model.matrix(fit)))
  expect_equal(value$object, matrix(fitted(fit, type = "link"), 10, 5), ignore_attr = TRUE)
  set.seed(2)
  response <- conditional(fit, output = c("all", "latent", "random"), samples = 5, local = settings, type = "response")
  expect_equal(response$newdata, exp(value$newdata))
  expect_identical(response$latent, value$latent)
  expect_identical(response$random, value$random)
  fit$ssn.object$preds <- list()
  local_mocked_bindings(get_dist_pred_object = function(...) stop("Unexpected prediction read"))
  expect_equal(dim(conditional(fit, output = "latent", samples = 1, local = settings)), c(10, 1))
  expect_true(all(is.finite(conditional(fit, output = "random", samples = 1, local = settings))))
  local_mocked_bindings(get_conditional_vecchia_covariance = function(...) stop("Unexpected covariance context"))
  expect_equal(dim(conditional(fit, output = "beta", samples = 1, local = settings)), c(2, 1))
})

test_that("observed anisotropy ordering differs from raw prediction ordering as specified", {
  fit <- conditional_vecchia_fixture(TRUE)
  fit$anisotropy <- TRUE
  fit$coefficients$params_object$euclid[["rotate"]] <- .4
  fit$coefficients$params_object$euclid[["scale"]] <- .3
  seen <- list()
  local_mocked_bindings(get_decorrelate_order = function(rows, ordering, data) {
    seen[[length(seen) + 1L]] <<- sf::st_coordinates(data)
    rev(seq_len(NROW(rows)))
  })
  settings <- get_local_list_conditional(list(size = 2, ordering = "coordinate"), fit, fit$ssn.object$preds$pred1km)
  context <- get_conditional_vecchia_glm_context(fit, "pred1km", settings)
  xy <- sf::st_coordinates(fit$ssn.object$obs)
  transformed <- transform_anis(xy[, 1], xy[, 2], .4, .3)
  expect_equal(seen[[1]], sf::st_coordinates(fit$ssn.object$preds$pred1km))
  expect_equal(unname(seen[[2]]), unname(cbind(transformed$xcoord_val, transformed$ycoord_val)))
  expect_identical(context$order, 12:1)
  expect_identical(context$prediction$order, 4:1)
  oracle <- conditional_vecchia_oracle(fit, 2, 12:1, 4:1)
  for (i in 1:12) expect_equal(context$operators[[i]]$variance, oracle$operators[[i]]$variance, tolerance = 1e-10)
})

test_that("spatial Vecchia operators allow empty eligible neighborhoods and stable cancellation", {
  fit <- conditional_vecchia_fixture(FALSE, TRUE)
  fit$ssn.object$preds$pred1km$part <- factor(c("new", "two", "one", "one"))
  covariance <- get_conditional_vecchia_covariance(fit, "pred1km", spatial = TRUE)
  settings <- get_local_list_conditional(list(size = 2, ordering = "none"), fit, fit$ssn.object$preds$pred1km)
  prediction <- prepare_conditional_vecchia(covariance, settings)
  expect_length(prediction$operators[[1]]$neighbors, 0)
  expect_equal(prediction$operators[[1]]$variance, as.numeric(conditional_vecchia_covariance(covariance, 13)))
  V <- matrix(c(1, 1 - 1e-12, 1 - 1e-12, 1), 2)
  op <- prepare_conditional_local_weights(V, 1, 2, "test")
  expect_gt(op$variance, 0)
  expect_equal(op$variance, 1 - (1 - 1e-12)^2, tolerance = 1e-15)
  expect_error(prepare_conditional_local_weights(matrix(c(1, 2, 2, 1), 2), 1, 2, "test"), "negative or nonfinite")
  gaussian <- prepare_conditional_local_gaussian(diag(2), matrix(1, 2, 1), c(0, 0), c(.1, -.2))
  expect_equal(gaussian$V, diag(1 / c(.9, 1.2)))
})

test_that("prediction covariance ties use positions in the reordered sequence", {
  Sigma <- matrix(.5, 4, 4)
  diag(Sigma) <- 1
  Sigma[1, 2:4] <- Sigma[2:4, 1] <- .05
  local_mocked_bindings(conditional_vecchia_covariance = function(context, rows, columns = rows) {
    Sigma[rows, columns, drop = FALSE]
  })
  context <- list(n = 1, newdata = data.frame(x = 1:3), partition = NULL)
  prediction <- prepare_conditional_vecchia(context, list(order = c(2L, 1L, 3L), size = 1, method = "covariance"))
  expect_equal(prediction$operators[[3]]$neighbors, 2)
  expect_equal(prediction$operators[[3]]$weights, .5)
  expect_equal(prediction$operators[[3]]$variance, .75)
})
