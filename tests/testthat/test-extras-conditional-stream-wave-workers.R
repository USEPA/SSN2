skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")
skip_if_not_installed("pracma")

test_that("installed wave simulation workers preserve conditional draws and unconditional moments", {
  skip_if(requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("SSN2"),
          "installed-package worker coverage")
  network <- mf04p
  network$preds$pred1km <- network$preds$pred1km[1:6, ]
  network$obs$count <- rep(1:5, length.out = nrow(network$obs))
  parameters <- list(
    tailup_params = tailup_params("wave", .2, 3000),
    taildown_params = taildown_params("wave", .3, 5000),
    euclid_params = euclid_params("none"), nugget_params = nugget_params("nugget", .25))
  local <- list(approximation = "low-rank", size_base = 10, size_new = 2,
    reorder_base = "none", reorder_new = "none", kmeans_new = FALSE, ncores = 2)
  for (glm in c(FALSE, TRUE)) {
    arguments <- list(formula = count ~ 1, ssn.object = network, additive = "afvArea",
      tailup_initial = tailup_initial("wave", .2, 3000, known = "given"),
      taildown_initial = taildown_initial("wave", .3, 5000, known = "given"),
      nugget_initial = nugget_initial("nugget", .25, known = "given"))
    if (glm) arguments$family <- "poisson"
    fit <- do.call(if (glm) ssn_glm else ssn_lm, arguments)
    set.seed(2)
    serial <- conditional(fit, "pred1km", samples = 7, output = "all", local = local)
    set.seed(2)
    parallel <- conditional(fit, "pred1km", samples = 7, output = "all",
      local = c(local, list(parallel = TRUE)))
    expect_identical(parallel, serial)
    expect_true(all(is.finite(parallel$newdata)))
  }
  cluster <- make_ssn_cluster
  local_mocked_bindings(make_ssn_cluster = function(ncores) {
    cl <- cluster(ncores)
    parallel::clusterSetRNGStream(cl, iseed = 2)
    cl
  })
  S <- as.matrix(covmatrix(fit))
  groups <- get_ssn_simulate_local_lowrank(local, nrow(S), network$obs)$index
  base <- groups$base
  expected <- S[, base] %*% solve(S[base, base], S[base, ])
  for (group in c(list(base), groups$new)) expected[group, group] <- S[group, group]
  samples <- 4000
  arguments <- c(list(ssn.object = network, additive = "afvArea", samples = samples), parameters)
  for (parallel in c(FALSE, TRUE)) {
    set.seed(2)
    draws <- do.call(ssn_rnorm, c(arguments, list(local = modifyList(local, list(parallel = parallel)))))
    expect_true(all(is.finite(draws)))
    expect_lt(max(abs(rowMeans(draws)) / sqrt(diag(expected) / samples)), 6)
    mcse <- sqrt((outer(diag(expected), diag(expected)) + expected^2) / (samples - 1))
    expect_lt(max(abs(cov(t(draws)) - expected) / mcse), 6)
  }
})
