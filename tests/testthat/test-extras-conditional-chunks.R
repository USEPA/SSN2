skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_chunks_fixture <- function(glm = FALSE, effects = FALSE) {
  network <- mf04p
  network$obs <- network$obs[round(seq(1, NROW(network$obs), length.out = 8)), ]
  network$preds <- list(pred1km = network$preds$pred1km[c(5, 1, 8, 3), ])
  network$obs$x <- seq(-1, 1, length.out = 8)
  network$obs$off <- seq(-.2, .2, length.out = 8)
  network$obs$y <- c(1, 3, 2, 4, 2, 5, 4, 3)
  network$obs$g <- factor(c("a", "a", "b", "b", "c", "c", "a", "b"))
  network$preds$pred1km$x <- c(-.5, .2, .8, 0)
  network$preds$pred1km$off <- c(.1, -.3, .2, .4)
  network$preds$pred1km$g <- factor(c("c", "new", "c", "new"))
  args <- list(formula = y ~ x + offset(off), ssn.object = network, additive = "afvArea",
    tailup_initial = tailup_initial("none"), taildown_initial = taildown_initial("none"),
    euclid_initial = euclid_initial("exponential", .4, 12000, known = "given"),
    nugget_initial = nugget_initial("nugget", .25, known = "given"))
  if (effects) {
    args$random <- ~g
    args$randcov_initial <- randcov_initial(g = .6, known = "given")
  }
  if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
  do.call(if (glm) ssn_glm else ssn_lm, args)
}

conditional_chunks_settings <- function(approximation, chunk) {
  if (approximation == "low-rank") {
    list(approximation = approximation, size_base = 4, reorder_base = "none",
      size_new = 2, reorder_new = "none", kmeans_new = FALSE, chunk_size = chunk)
  } else list(approximation = approximation, size = 3, ordering = "none", chunk_size = chunk)
}

conditional_chunks_oracle <- function(fit, approximation) {
  glm <- inherits(fit, "ssn_glm")
  base <- if (approximation == "low-rank") 1:4 else 1:8
  n <- length(base)
  X <- model.matrix(fit)[base, , drop = FALSE]
  Xp <- cbind(1, fit$ssn.object$preds$pred1km$x)
  Sigma <- as.matrix(covmatrix(fit))[base, base]
  spatial <- conditional_spatial_object(fit)
  S <- as.matrix(covmatrix(spatial))[base, base]
  cross <- as.matrix(covmatrix(spatial, "pred1km", cov_type = "obs.pred"))[base, , drop = FALSE]
  pred <- as.matrix(covmatrix(spatial, "pred1km", cov_type = "pred.pred"))
  eta <- if (glm) fitted(fit, type = "link") else model.response(model.frame(fit))
  latent_mean <- if (glm) eta + fit$laplace$delta[names(eta)] else eta
  center <- (latent_mean - model.offset(model.frame(fit)))[base]
  beta <- if (glm) fit$laplace$betahat else coef(fit)
  V <- matrix(0, n, n)
  M <- matrix(0, n, NCOL(X))
  if (approximation == "low-rank") {
    if (glm) {
      V <- solve(solve(Sigma) + diag(exp(eta[base])))
      M <- V %*% solve(Sigma, X)
    }
    W <- t(cross) %*% solve(S)
    R <- matrix(0, 4, 4)
    for (b in list(1:2, 3:4)) R[b, b] <- pred[b, b] - W[b, ] %*% cross[, b]
  } else {
    select <- function(i, pool, size, covariance) {
      pool <- pool[pool != i]
      head(pool[order(-abs(covariance[i, pool]), pool)], size)
    }
    F <- matrix(0, n, n)
    if (glm) for (i in seq_len(n)) {
      N <- select(i, seq_len(i - 1L), 3, Sigma)
      A <- unique(c(i, N, select(i, seq_len(n), 2, Sigma)))
      P <- solve(Sigma[A, A, drop = FALSE])
      VA <- solve(P + diag(exp(eta[A]), length(A)))
      MA <- VA %*% P %*% X[A, , drop = FALSE]
      j <- match(N, A)
      g <- if (length(N)) as.vector(solve(VA[j, j, drop = FALSE], VA[j, 1])) else numeric()
      M[i, ] <- MA[1, ] - g %*% MA[j, , drop = FALSE] + g %*% M[N, , drop = FALSE]
      F[i, ] <- g %*% F[N, , drop = FALSE]
      F[i, i] <- sqrt(VA[1, 1] - sum(g * VA[j, 1]))
    }
    V <- tcrossprod(F)
    full <- rbind(cbind(S, cross), cbind(t(cross), pred))
    W <- rbind(diag(n), matrix(0, 4, n))
    E <- matrix(0, n + 4, 4)
    for (k in 1:4) {
      i <- n + k
      N <- select(i, seq_len(i - 1L), 3, full)
      g <- as.vector(solve(full[N, N, drop = FALSE], full[N, i]))
      W[i, ] <- g %*% W[N, , drop = FALSE]
      E[i, ] <- g %*% E[N, , drop = FALSE]
      E[i, k] <- sqrt(full[i, i] - sum(g * full[N, i]))
    }
    W <- W[n + 1:4, , drop = FALSE]
    R <- tcrossprod(E[n + 1:4, , drop = FALSE])
  }
  if (!is.null(fit$random)) {
    groups <- c("a", "b", "c", "new")
    Z <- outer(as.character(fit$ssn.object$obs$g[base]), groups, "==") * 1
    Zp <- outer(as.character(fit$ssn.object$preds$pred1km$g), groups, "==") * 1
    H <- solve(diag(1 / .6, 4) + t(Z) %*% solve(S, Z))
    A <- Zp - W %*% Z
    R <- R + A %*% H %*% t(A)
    W <- W + A %*% H %*% t(Z) %*% solve(S)
  }
  B <- Xp + W %*% (M - X)
  C <- as.matrix(vcov(fit))
  mu <- as.vector(Xp %*% beta + W %*% (center - X %*% beta) +
    fit$ssn.object$preds$pred1km$off)
  covariance <- R + W %*% V %*% t(W) + B %*% C %*% t(B)
  list(mean = c(as.vector(beta), mu), covariance = rbind(cbind(C, C %*% t(B)),
    cbind(B %*% C, covariance)))
}

expect_chunks_moments <- function(draws, reference, label) {
  count <- NCOL(draws)
  V <- reference$covariance
  expect_lt(max(abs(rowMeans(draws) - reference$mean) / sqrt(diag(V) / count)), 6, label = label)
  mcse <- sqrt((outer(diag(V), diag(V)) + V^2) / (count - 1))
  expect_lt(max(abs(cov(t(draws)) - V) / mcse), 6, label = label)
  h <- c(0, 0, rep(.25, 4))
  variance <- as.numeric(t(h) %*% V %*% h)
  expect_lt(abs(var(as.vector(h %*% draws)) - variance) / (variance * sqrt(2 / (count - 1))),
            6, label = paste(label, "prediction average"))
}

conditional_chunks_probe <- function(fit, settings, samples = 7) {
  state <- new.env(parent = emptyenv())
  state$counts <- integer()
  state$arrays <- list()
  namespace <- asNamespace("SSN2")
  counted <- function(name) {
    original <- get(name, namespace)
    state$counts[name] <- 0L
    function(...) {
      state$counts[name] <- state$counts[name] + 1L
      original(...)
    }
  }
  allocated <- function(name) {
    original <- get(name, namespace)
    env <- new.env(parent = environment(original))
    factor_name <- paste0("chol:", name)
    state$counts[factor_name] <- 0L
    env$chol <- function(...) {
      state$counts[factor_name] <- state$counts[factor_name] + 1L
      base::chol(...)
    }
    env$matrix <- function(...) {
      value <- base::matrix(...)
      args <- list(...)
      state$arrays[[length(state$arrays) + 1L]] <- list(name = name, dim = dim(value),
        initialized = !all(is.na(args[[1L]])))
      value
    }
    environment(original) <- env
    original
  }
  preparations <- c("get_conditional_vecchia_covariance", "get_conditional_lowrank_covariance",
    "get_conditional_glm_base", "get_conditional_vecchia_glm_context", "get_D",
    "prepare_conditional_random", "prepare_conditional_vecchia", "prepare_conditional_local_gaussian",
    "prepare_conditional_local_weights", "conditional_local_factor", "factor_conditional_glm_block")
  draws <- c("get_conditional_lm_lowrank", "get_conditional_glm_lowrank", "draw_conditional_lowrank_blocks",
    "get_conditional_vecchia_lm", "get_conditional_vecchia_glm", "draw_conditional_vecchia",
    "draw_conditional_random")
  bindings <- c(setNames(lapply(preparations, counted), preparations),
                setNames(lapply(draws, allocated), draws))
  do.call(local_mocked_bindings, c(bindings, list(.package = "SSN2", .env = environment())))
  set.seed(2)
  value <- conditional(fit, "pred1km", samples = samples, output = c("beta", "newdata"), local = settings)
  list(value = value, counts = state$counts, arrays = state$arrays)
}

test_that("chunk boundaries preserve joint moments across model and approximation paths", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) for (effects in c(FALSE, TRUE)) {
    fit <- conditional_chunks_fixture(glm, effects)
    for (approximation in c("low-rank", "vecchia")) {
      reference <- conditional_chunks_oracle(fit, approximation)
      for (chunk in c(1, 17, 1201)) {
        set.seed(2)
        value <- conditional(fit, "pred1km", samples = 1200, output = c("beta", "newdata"),
          local = conditional_chunks_settings(approximation, chunk))
        label <- paste(glm, effects, approximation, chunk)
        expect_equal(dim(value$beta), c(2, 1200), label = label)
        expect_equal(dim(value$newdata), c(4, 1200), label = label)
        expect_true(all(is.finite(value$newdata)), label = label)
        expect_chunks_moments(rbind(value$beta, value$newdata), reference, label)
      }
    }
  }
})

test_that("chunk loops reuse preparations and bound non-output draw arrays", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) for (effects in c(FALSE, TRUE)) {
    fit <- conditional_chunks_fixture(glm, effects)
    for (approximation in c("low-rank", "vecchia")) {
      baseline <- NULL
      for (chunk in c(1, 3, 20)) {
        probe <- conditional_chunks_probe(fit, conditional_chunks_settings(approximation, chunk))
        label <- paste(glm, effects, approximation, chunk)
        if (is.null(baseline)) baseline <- probe$counts
        expect_identical(probe$counts, baseline, label = label)
        expect_equal(unname(probe$counts["get_conditional_vecchia_covariance"]), 1, label = label)
        expect_equal(unname(probe$counts["get_D"]), as.integer(glm), label = label)
        expect_equal(unname(probe$counts["prepare_conditional_random"]), as.integer(effects), label = label)
        lowrank <- approximation == "low-rank"
        expect_equal(unname(probe$counts["get_conditional_lowrank_covariance"]), as.integer(lowrank), label = label)
        expect_equal(unname(probe$counts["factor_conditional_glm_block"]), if (lowrank) 2 else 0, label = label)
        expect_equal(unname(probe$counts["get_conditional_glm_base"]), as.integer(glm && lowrank), label = label)
        expect_equal(unname(probe$counts["get_conditional_vecchia_glm_context"]), as.integer(glm && !lowrank), label = label)
        expect_equal(unname(probe$counts["prepare_conditional_vecchia"]), as.integer(!lowrank), label = label)
        expect_equal(unname(probe$counts["prepare_conditional_local_gaussian"]),
          if (glm && !lowrank) 8 else 0, label = label)
        initialized <- Filter(function(a) a$initialized, probe$arrays)
        expect_gt(length(initialized), 0, label = label)
        widths <- vapply(initialized, function(a) a$dim[2], integer(1))
        expect_lte(max(widths), min(chunk, 7), label = label)
        expect_setequal(unique(widths), if (chunk > 7) 7L else unique(c(chunk, 7L %% chunk)[c(chunk, 7L %% chunk) > 0]))
        if (!lowrank) {
          internal <- Filter(function(a) a$dim[1] %in% c(8, 12), probe$arrays)
          expect_gt(length(internal), 0, label = label)
          expect_true(all(vapply(internal, function(a) a$dim[2] <= min(chunk, 7), logical(1))), label = label)
        } else {
          # Base residuals and shared effects may retain all simulation columns.
          retained <- Filter(function(a) a$name %in% c("get_conditional_lm_lowrank", "get_conditional_glm_lowrank") &&
            !a$initialized && identical(a$dim, c(4L, 7L)), probe$arrays)
          expect_gte(length(retained), 1, label = label)
        }
      }
    }
  }
})

test_that("parallel low-rank preparation preserves the parent RNG stream", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) for (effects in c(FALSE, TRUE)) {
    fit <- conditional_chunks_fixture(glm, effects)
    settings <- conditional_chunks_settings("low-rank", 3)
    set.seed(2)
    serial <- conditional(fit, "pred1km", samples = 7, output = c("beta", "newdata"), local = settings)
    serial_seed <- .Random.seed
    settings$parallel <- TRUE
    settings$ncores <- 2
    set.seed(2)
    parallel <- conditional(fit, "pred1km", samples = 7, output = c("beta", "newdata"), local = settings)
    expect_identical(parallel, serial)
    expect_identical(.Random.seed, serial_seed)
    settings <- conditional_chunks_settings("vecchia", 3)
    set.seed(2)
    serial <- conditional(fit, "pred1km", samples = 7, output = c("beta", "newdata"), local = settings)
    settings$parallel <- TRUE
    set.seed(2)
    expect_warning(parallel <- conditional(fit, "pred1km", samples = 7,
      output = c("beta", "newdata"), local = settings), "local\\$parallel is not used")
    expect_identical(parallel, serial)
  }
})

test_that("singleton chunks retain output dimensions and requested latent fields", {
  set.seed(2)
  for (glm in c(FALSE, TRUE)) for (effects in c(FALSE, TRUE)) {
    fit <- conditional_chunks_fixture(glm, effects)
    fit$ssn.object$preds$pred1km <- fit$ssn.object$preds$pred1km[1, ]
    for (approximation in c("low-rank", "vecchia")) {
      set.seed(2)
      value <- conditional(fit, "pred1km", samples = 1, output = c("beta", "newdata"),
        local = conditional_chunks_settings(approximation, 1))
      expect_equal(dim(value$beta), c(2, 1))
      expect_equal(dim(value$newdata), c(1, 1))
      expect_true(all(is.finite(value$newdata)))
    }
    if (glm) for (chunk in c(1, 3, 20)) {
      set.seed(2)
      value <- conditional(fit, samples = 7, output = c("beta", "latent", if (effects) "random"),
        local = conditional_chunks_settings("vecchia", chunk))
      expect_equal(dim(value$latent), c(8, 7))
      expect_true(all(is.finite(value$latent)))
      if (effects) {
        expect_equal(dim(value$random), c(3, 7))
        expect_false(is.null(attr(value$random, "effect_info")))
      }
    }
  }
})
