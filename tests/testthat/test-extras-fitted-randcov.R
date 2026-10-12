skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

fitted_randcov_data <- function(obs, params, partition = NULL,
                                block = rep(1L, NROW(obs)), permutation = NULL) {
  ordered <- if (is.null(permutation)) obs else obs[permutation, , drop = FALSE]
  index <- if (is.null(permutation)) block else block[permutation]
  Z <- lapply(names(params), function(term) get_randcov_Z(term, ordered, ZZt = FALSE))
  names(Z) <- names(params)
  list(ssn.object = list(obs = obs), partition_factor_initial = partition,
       partition_factor = partition, local_index = index,
       order_bigdata = permutation,
       randcov_list = lapply(split(seq_len(NROW(obs)), index), function(i) {
         lapply(Z, function(z) list(Z = z$Z[i, , drop = FALSE]))
       }))
}

fitted_randcov_fixture <- function(glm, local, partitioned) {
  set.seed(2)
  network <- mf04p
  network$obs <- network$obs[round(seq(1, NROW(network$obs), length.out = 24)), ]
  d <- expand.grid(g = c("b", "a"), part = c("two", "one"),
                   block = c(7, 2, 10), replicate = 1:2)
  network$obs$g <- factor(d$g, levels = c("b", "a", "unused"))
  network$obs$part <- factor(d$part, levels = c("two", "one"))
  network$obs$block <- d$block
  network$obs$x <- seq(-1.2, 1.1, length.out = 24)
  network$obs$off <- rep(c(-.3, .2, .1), 8)
  network$obs$y <- c(1, 4, 2, 6, 3, 2, 5, 1, 7, 4, 2, 5,
                     3, 6, 1, 4, 2, 8, 5, 3, 4, 2, 6, 1)
  network$obs$y[c(3, 20)] <- NA
  network$obs <- network$obs[sample.int(24), ]
  args <- list(formula = y ~ x + offset(off), ssn.object = network,
    random = ~(x | g), partition_factor = if (partitioned) ~part else NULL,
    randcov_initial = randcov_initial(`1 | g` = .3, `x | g` = .15, known = "given"),
    tailup_initial = tailup_initial("none"), taildown_initial = taildown_initial("none"),
    euclid_initial = euclid_initial("none"),
    nugget_initial = nugget_initial("nugget", .25, known = "given"),
    local = if (local) list(index = network$obs$block, var_adjust = "none") else FALSE)
  if (glm) args$family <- "poisson" else args$ddf <- "asymptotic"
  do.call(if (glm) ssn_glm else ssn_lm, args)
}

for (glm in c(FALSE, TRUE)) for (local in c(FALSE, TRUE)) {
  test_that(paste(if (glm) "Poisson" else "LM", if (local) "local" else "global",
                  "fitted effects agree with independent analytical BLUPs"), {
    for (partitioned in c(FALSE, TRUE)) {
      fit <- fitted_randcov_fixture(glm, local, partitioned)
      obs <- fit$ssn.object$obs
      n <- NROW(obs)
      expect_equal(n, 22)
      expect_length(fit$missing_index, 2)
      expect_true(is.unsorted(as.numeric(ssn_get_netgeom(obs)$pid)))
      expect_gt(length(unique(ssn_get_netgeom(obs)$NetworkID)), 1)
      expect_equal(as.numeric(fit$coefficients$params_object$randcov), c(.3, .15))
      groups <- if (partitioned) c("b", "a") else c("a", "b")
      partitions <- if (partitioned) c("two", "one") else "all"
      masks <- list()
      for (g in groups) for (part in partitions) {
        label <- if (partitioned) paste0("g", g, ":part", part) else g
        masks[[label]] <- as.numeric(obs$g == g &
          if (partitioned) obs$part == part else rep(TRUE, n))
      }
      Z <- do.call(cbind, masks)
      Zs <- obs$x * Z
      Sigma <- diag(.25, n) + .3 * tcrossprod(Z) + .15 * tcrossprod(Zs)
      expect_equal(as.matrix(covmatrix(fit)), Sigma, ignore_attr = TRUE, tolerance = 1e-10)
      X <- cbind(1, obs$x)
      residual <- (if (glm) as.numeric(fitted(fit, type = "link")) else obs$y) -
        obs$off - as.vector(X %*% coef(fit))
      blocks <- if (local) split(seq_len(n), fit$local_index) else list(seq_len(n))
      if (local) {
        expect_equal(length(blocks), 3)
        expect_true(all(vapply(blocks, function(i) all(colSums(Z[i, , drop = FALSE]) > 0), logical(1))))
        expect_equal(outer(fit$local_index, fit$local_index, `==`), outer(obs$block, obs$block, `==`))
      }
      contributions <- lapply(blocks, function(i) {
        precision_residual <- solve(Sigma[i, i, drop = FALSE], residual[i])
        cbind(.3 * as.vector(crossprod(Z[i, , drop = FALSE], precision_residual)),
              .15 * as.vector(crossprod(Zs[i, , drop = FALSE], precision_residual)))
      })
      actual <- fitted(fit, type = "randcov")
      expect_named(actual, c("1 | g", "x | g"))
      for (j in 1:2) {
        values <- vapply(contributions, function(x) x[, j], numeric(NCOL(Z)))
        expected <- apply(values, 1, function(x) if (any(x != 0)) mean(x[x != 0]) else 0)
        expect_identical(names(actual[[j]]), colnames(Z))
        expect_equal(actual[[j]], setNames(expected, colnames(Z)), tolerance = 1e-8)
        expect_true(all(is.finite(actual[[j]])))
        if (local) expect_gt(max(abs(rowSums(values) - expected)), 1e-5)
      }
    }
  })
}

test_that("fitted effects average nonzero block contributions including cancellations", {
  obs <- data.frame(g = factor(rep(c("b", "a", "c"), 3), levels = c("c", "b", "a", "unused")),
                    part = "one", netID = rep(1:3, each = 3))
  block <- rep(c(2, 7, 10), each = 3)
  r <- list(c(2, 2, 0), c(-2, 4, 0), c(0, 0, 0))
  for (partition in list(NULL, ~part)) {
    object <- fitted_randcov_data(obs, c("1 | g" = 1), partition, block)
    actual <- get_fitted_randcov(object, c("1 | g" = 1), r)[[1]]
    expected <- if (is.null(partition)) c(a = 3, b = 0, c = 0) else
      c("gc:partone" = 0, "gb:partone" = 0, "ga:partone" = 3)
    expect_identical(actual, expected)
    expect_true(all(is.finite(actual)))
    zero <- get_fitted_randcov(object, c("1 | g" = 0), r)[[1]]
    expect_identical(zero, setNames(numeric(length(expected)), names(expected)))
    zero_residual <- get_fitted_randcov(object, c("1 | g" = 1), rep(list(numeric(3)), 3))[[1]]
    expect_identical(zero_residual, zero)
  }
})

test_that("unpartitioned fitted names and support retain the stored Z structure", {
  obs <- data.frame(g = factor(c("z", "a", "z", "b", "a", "b"),
                              levels = c("z", "b", "a", "unused")),
                    x = c(0, 2, 0, -1, 3, 1), netID = c(1, 1, 2, 2, 3, 3))
  params <- c("1 | g" = 2, "x | g" = .5)
  object <- fitted_randcov_data(obs, params)
  r <- c(1, 2, 3, 4, -2, 4)
  expected <- list("1 | g" = c(a = 0, b = 16, z = 8), "x | g" = c(a = -1, b = 0))
  expect_identical(get_fitted_randcov(object, params, list(r)), expected)
  object$partition_factor <- ~netID
  expect_identical(get_fitted_randcov(object, params, list(r)), expected)
  expect_null(get_fitted_randcov(object, NULL, list(r)))
  expect_null(get_fitted_randcov(object, numeric(), list(r)))
  # A stored supported zero column must survive even when absent from the data.
  for (term in names(params)) {
    Z <- object$randcov_list[[1]][[term]]$Z
    object$randcov_list[[1]][[term]]$Z <- cbind(Z, retained = numeric(NROW(obs)))
    expected[[term]] <- c(expected[[term]], retained = 0)
    expected[[term]] <- expected[[term]][order(names(expected[[term]]))]
  }
  actual <- get_fitted_randcov(object, params, list(r))
  expect_identical(actual, expected)
})

test_that("partitioned helper aligns local permutations and ignores global reporting order", {
  obs <- data.frame(g = c("a", "b", "a", "b", "a", "b", "a", "b"),
                    part = c("one", "two", "two", "one", "one", "two", "two", "one"),
                    netID = c(1, 2, 1, 2, 2, 1, 2, 1))
  block <- c(7, 2, 7, 2, 2, 7, 2, 7)
  permutation <- c(7, 2, 5, 4, 8, 1, 6, 3)
  r <- c(1, 2, 4, 8, 16, 32, 64, 128)
  params <- c("1 | g" = 2)
  object <- fitted_randcov_data(obs, params, ~part, block, permutation)
  expect_null(object$obdata)
  actual <- get_fitted_randcov(object, params, split(r[permutation], block[permutation]))[[1]]
  expect_equal(actual, c("ga:partone" = 17, "ga:parttwo" = 68,
                         "gb:partone" = 136, "gb:parttwo" = 34))
  object <- fitted_randcov_data(obs, params, ~part)
  object$order <- rev(seq_len(NROW(obs)))
  expected <- c("ga:partone" = 34, "ga:parttwo" = 136,
                "gb:partone" = 272, "gb:parttwo" = 68)
  expect_equal(get_fitted_randcov(object, params, list(r))[[1]], expected)
  object$partition_factor_initial <- NULL
  expect_equal(get_fitted_randcov(object, params, list(r))[[1]], expected)
  object$local_index <- NULL
  expect_equal(get_fitted_randcov(object, params, list(r))[[1]], expected)
})

test_that("nested slopes and compound partitions use independent grouping tuples", {
  obs <- expand.grid(g = c("a", "b"), h = c("u", "v"),
                     p = c("one", "two"), q = c("left", "right"), block = c(2, 7))
  obs$x <- rep(c(-2, 0, 1, 3), 8)
  params <- c("1 | g" = .5, "x | g" = .25, "1 | g:h" = 2, "x | g:h" = 1)
  expect_identical(get_randcov_names(~(x | g/h)), names(params))
  object <- fitted_randcov_data(obs, params, ~p:q, obs$block, rev(seq_len(NROW(obs))))
  r <- seq_len(NROW(obs)) - 9
  actual <- get_fitted_randcov(object, params,
    split(rev(r), rev(obs$block)))
  expect_named(actual, names(params))
  for (term in names(params)) {
    compound <- grepl("g:h", term, fixed = TRUE)
    slope <- startsWith(term, "x")
    expected <- numeric()
    for (g in c("a", "b")) for (h in if (compound) c("u", "v") else "all") {
      for (p in c("one", "two")) for (q in c("left", "right")) {
        mask <- obs$g == g & obs$p == p & obs$q == q &
          if (compound) obs$h == h else rep(TRUE, NROW(obs))
        z <- mask * if (slope) obs$x else 1
        values <- vapply(c(2, 7), function(b) params[[term]] * sum(z[obs$block == b] * r[obs$block == b]), 0)
        label <- paste0("g", g, if (compound) paste0(":h", h), ":p", p, ":q", q)
        expected[label] <- if (any(values != 0)) mean(values[values != 0]) else 0
      }
    }
    expect_equal(actual[[term]][names(expected)], expected)
    expect_setequal(names(actual[[term]]), names(expected))
    expect_true(all(is.finite(actual[[term]])))
  }
})

test_that("colliding partition display labels never merge distinct effect tuples", {
  obs <- data.frame(g = c("a:partb", "a", "a:partb", "a"),
                    part = c("c", "b:partc", "c", "b:partc"),
                    netID = c(1, 1, 2, 2))
  params <- c("1 | g" = 1)
  object <- fitted_randcov_data(obs, params, ~part, c(2, 2, 7, 7))
  actual <- get_fitted_randcov(object, params, list(c(2, 10), c(4, 14)))[[1]]
  expect_identical(actual, c("ga:partb:partc" = 12, "ga:partb:partc.1" = 3))
  expect_identical(anyDuplicated(names(actual)), 0L)
})
