skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"))

conditional_laplace_fixture <- local({
  network <- NULL
  fits <- list()
  function(iid = FALSE, local = FALSE, effects = FALSE) {
    key <- paste(iid, local, effects)
    if (!is.null(fits[[key]])) return(fits[[key]])
    if (is.null(network)) {
      s <- mf04p
      root <- tempfile("conditional-laplace-")
      dir.create(root)
      stopifnot(file.copy(s$path, root, recursive = TRUE, copy.mode = FALSE))
      s$path <- file.path(root, basename(s$path))
      s$obs <- s$obs[round(seq(1, NROW(s$obs), length.out = 16)), ]
      s$preds <- list(correction = s$preds$pred1km[c(5, 1, 8, 3), ])
      set.seed(2)
      s$obs$x <- rnorm(16)
      s$obs$off <- seq(-.2, .3, length.out = 16)
      s$obs$y <- c(1, 3, 2, 4, 2, 5, 4, 2, 5, 3, 6, 4, 2, 5, 1, 4)
      s$obs$g <- factor(rep(letters[1:4], 4))
      s$obs$part <- factor(rep(c("one", "two"), 8))
      s$preds$correction$x <- c(-.5, .2, .8, 0)
      s$preds$correction$off <- c(.1, -.3, .2, .4)
      s$preds$correction$g <- factor(c("a", "b", "c", "d"))
      s$preds$correction$part <- factor(c("one", "two", "one", "two"))
      ssn_create_distmat(s, predpts = "correction", among_predpts = TRUE, overwrite = TRUE)
      network <<- s
    }
    args <- list(formula = y ~ x + offset(off), ssn.object = network, family = "poisson", additive = "afvArea",
      tailup_initial = if (iid) tailup_initial("none") else tailup_initial("exponential", .1, 10000, known = "given"),
      taildown_initial = if (iid) taildown_initial("none") else taildown_initial("exponential", .1, 10000, known = "given"),
      euclid_initial = if (iid) euclid_initial("none") else euclid_initial("exponential", .2, 12000, known = "given"),
      nugget_initial = nugget_initial("nugget", .25, known = "given"),
      local = if (local) list(index = rep(c(2, 1), 8), var_adjust = "none") else FALSE)
    if (effects) {
      args$random <- ~g
      args$randcov_initial <- randcov_initial(g = .1, known = "given")
      args$partition_factor <- ~part
    }
    fit <- do.call(ssn_glm, args)
    fits[[key]] <<- fit
    fit
  }
})

test_that("conditional corrections shift draws without changing their random errors", {
  paths <- list(exact = FALSE,
    lowrank = list(approximation = "low-rank", method_base = "all", method_new = "all"),
    vecchia = list(approximation = "vecchia", size = 5, ordering = "none"))
  for (fitting_local in c(FALSE, TRUE)) for (effects in c(FALSE, TRUE)) {
    fit <- conditional_laplace_fixture(local = fitting_local, effects = effects)
    original <- fit
    old <- fit
    old$laplace <- NULL
    for (path in paths) {
      args <- list(object = fit, newdata = "correction", samples = 8, local = path,
        output = c("newdata", "beta", "latent", "object"))
      set.seed(2)
      on <- do.call(conditional, args)
      set.seed(2)
      off <- do.call(conditional, c(args, list(fit_correct = FALSE)))
      expect_identical(on$object, off$object)
      expect_equal(on$object, matrix(fitted(fit, "link"), fit$n, 8), ignore_attr = TRUE)
      for (component in c("beta", "latent", "newdata")) {
        difference <- on[[component]] - off[[component]]
        expect_equal(difference, matrix(difference[, 1], nrow(difference), ncol(difference)),
          ignore_attr = TRUE, tolerance = 1e-10)
        expect_equal(cov(t(on[[component]])), cov(t(off[[component]])), tolerance = 1e-10)
      }
      expect_equal(as.numeric(on$latent[, 1] - off$latent[, 1]), unname(fit$laplace$delta), tolerance = 1e-10)
      expect_equal(as.numeric(on$beta[, 1] - off$beta[, 1]), unname(fit$laplace$betahat - coef(fit)), tolerance = 1e-10)
      if (identical(path, FALSE)) {
        shift <- predict(fit, "correction") - predict(fit, "correction", fit_correct = FALSE)
        expect_equal(as.numeric(on$newdata[, 1] - off$newdata[, 1]), as.numeric(shift), tolerance = 1e-9)
      }
      args$object <- old
      set.seed(2)
      expect_identical(do.call(conditional, args), off)
    }
    expect_identical(fit, original)
  }
})

test_that("independent conditional draws retain mode curvature and old-object fallback", {
  for (local in c(FALSE, TRUE)) {
    fit <- conditional_laplace_fixture(iid = TRUE, local = local)
    old <- fit
    old$laplace <- NULL
    output <- c("beta", "latent", "newdata", "object")
    set.seed(2)
    on <- conditional(fit, "correction", samples = 8, output = output)
    set.seed(2)
    off <- conditional(fit, "correction", samples = 8, output = output, fit_correct = FALSE)
    expect_identical(on$object, off$object)
    expect_equal(on$latent - off$latent, matrix(fit$laplace$delta, fit$n, 8), ignore_attr = TRUE)
    expect_equal(on$beta - off$beta, matrix(fit$laplace$betahat - coef(fit), fit$p, 8), ignore_attr = TRUE)
    for (component in c("beta", "latent", "newdata")) {
      expect_equal(cov(t(on[[component]])), cov(t(off[[component]])), tolerance = 1e-10)
    }
    set.seed(2)
    expect_identical(conditional(old, "correction", samples = 8, output = output), off)
    for (bad in list(NA, 1, c(TRUE, FALSE), NULL, "TRUE")) {
      expect_error(conditional(fit, "correction", samples = 1, fit_correct = bad), "fit_correct")
    }
  }
})

test_that("a reduced simulation base subsets the corrected latent center by pid", {
  fit <- conditional_laplace_fixture()
  base <- seq_len(8)
  settings <- list(approximation = "low-rank", size_base = 8, reorder_base = "none",
    size_new = 2, reorder_new = "none", kmeans_new = FALSE)
  output <- c("beta", "newdata", "object")
  set.seed(2)
  on <- conditional(fit, "correction", samples = 8, output = output, local = settings)
  set.seed(2)
  off <- conditional(fit, "correction", samples = 8, output = output, local = settings, fit_correct = FALSE)
  X <- model.matrix(fit)[base, , drop = FALSE]
  Xnew <- get_newdata_model_matrix(fit, fit$ssn.object$preds$correction)$newdata_model
  S <- as.matrix(covmatrix(fit))[base, base, drop = FALSE]
  C <- t(as.matrix(covmatrix(fit, "correction", cov_type = "obs.pred"))[base, , drop = FALSE])
  beta_shift <- fit$laplace$betahat - coef(fit)
  expected <- as.numeric(Xnew %*% beta_shift + C %*% solve(S,
    fit$laplace$delta[base] - X %*% beta_shift))
  expect_equal(on$newdata - off$newdata, matrix(expected, length(expected), 8),
    ignore_attr = TRUE, tolerance = 1e-9)
  expect_equal(cov(t(on$newdata)), cov(t(off$newdata)), tolerance = 1e-10)
  expect_identical(on$object, off$object)
  old <- fit
  old$laplace <- NULL
  set.seed(2)
  expect_identical(conditional(old, "correction", samples = 8, output = output, local = settings), off)
})

test_that("log-link conditional response means agree with averaged predictions", {
  for (iid in c(FALSE, TRUE)) {
    fit <- conditional_laplace_fixture(iid = iid)
    expected <- predict(fit, "correction", type = "response", local = FALSE)
    set.seed(2)
    link <- conditional(fit, "correction", samples = 6000, local = FALSE)
    set.seed(2)
    response <- conditional(fit, "correction", samples = 6000, local = FALSE, type = "response")
    expect_equal(response, exp(link))
    mcse <- apply(response, 1, sd) / sqrt(ncol(response))
    expect_lt(max(abs(rowMeans(response) - expected) / mcse), 6)
    prediction <- predict(fit, "correction", se.fit = TRUE, local = FALSE)
    se_mcse <- prediction$se.fit / sqrt(2 * (ncol(link) - 1))
    expect_lt(max(abs(apply(link, 1, sd) - prediction$se.fit) / se_mcse), 6)
  }
})
