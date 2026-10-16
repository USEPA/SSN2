skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"))

nugget_search_initial <- function(de = 1) {
  list(tailup_initial = tailup_initial("none"),
    taildown_initial = taildown_initial("none"),
    euclid_initial = euclid_initial("exponential", de = de, range = 1000, rotate = 0, scale = 1, known = "given"),
    nugget_initial = nugget_initial("nugget", .3))
}

test_that("sole-nugget search bounds use relevant covariance floors", {
  initial <- nugget_search_initial()
  for (scale in c(1e-6, 1, 10)) for (diagtol in c(0, 1e-4, .01)) {
    initial$euclid_initial$initial[["de"]] <- scale
    lower <- nugget_only_lower(c(nugget_log = log(.3)), initial, diagtol)
    expect_equal(lower, log(max(1e-4 * scale, diagtol)))
    params <- get_params_object_known(initial)
    params$nugget[["nugget"]] <- 0
    actual <- cov_matrix(params$nugget, list(network_index = c(1, 1)), de_scale = scale, diagtol = diagtol)
    expect_equal(as.numeric(diag(actual)), rep(exp(lower), 2))
  }
  initial <- nugget_search_initial(0)
  expect_null(nugget_only_lower(c(nugget_log = 0), initial, 0))
  expect_null(nugget_only_lower(c(euclid_range_log = 0), initial, 1e-4))
  expect_null(nugget_only_lower(c(nugget_log = 0, rand_log = 0), initial, 1e-4))
  initial$nugget_initial <- nugget_initial("none")
  expect_null(nugget_only_lower(c(nugget_log = 0), initial, 1e-4))
})

test_that("bounded nugget searches retain optimizer controls and diagnostics", {
  args <- list(par = c(nugget_log = -20), fn = function(theta) (exp(theta) - .7)^2)
  controls <- check_optim_method(args$par, get_optim_dotlist(hessian = TRUE, control = list(reltol = 1e-9)))
  bounded <- controls
  bounded$lower <- log(1e-4)
  expect_identical(optim_with_nugget_bound(args, controls, log(1e-4)), do.call(optim, c(args, bounded)))
  expect_identical(optim_with_nugget_bound(args, controls), do.call(optim, c(args, controls)))
  for (lower in list(NA_real_, Inf, c(-2, -1), 50, 51)) {
    expect_error(optim_with_nugget_bound(args, controls, lower), "Invalid bounds")
  }
  for (upper in list(NULL, NA_real_, Inf, c(0, 1), -20)) {
    bad <- controls
    bad$upper <- upper
    expect_error(optim_with_nugget_bound(args, bad, log(1e-4)), "Invalid bounds")
  }
  args$par <- c(range_log = 0)
  expect_error(optim_with_nugget_bound(args, controls, log(1e-4)), "Invalid bounds")
  args$par <- c(nugget_log = 0)
  args$fn <- function(theta) stop("objective failed")
  expect_error(optim_with_nugget_bound(args, controls, log(1e-4)), "objective failed")
})

gaussian_nugget_data <- function() {
  set.seed(2)
  s <- mf04p
  s$obs$.keep <- seq_len(nrow(s$obs)) <= 36
  s <- ssn_subset(s, tempfile(fileext = ".ssn"), .keep)
  ssn_create_distmat(s)
  s$obs$x <- rnorm(36)
  s$obs$y <- rnorm(36)
  s$obs$off <- seq(-.2, .2, length.out = 36)
  s$obs$g <- factor(rep(1:6, each = 6))
  s$obs$part <- factor(rep(1:2, each = 18))
  fit <- ssn_lm(y ~ x, s, tailup_type = "none", ddf = "asymptotic",
    euclid_initial = euclid_initial("exponential", de = 1, range = 1000, known = "given"),
    nugget_initial = nugget_initial("nugget", .5, known = "given"))
  s$obs$y <- as.numeric(1 + .5 * s$obs$x + t(chol(as.matrix(covmatrix(fit)))) %*% rnorm(36))
  s
}

test_that("Gaussian nugget searches recover reference optimums in different response units", {
  for (scale in c(1, .001)) {
    s <- gaussian_nugget_data()
    s$obs$y <- scale * s$obs$y
    for (method in c("ml", "reml")) {
      args <- list(formula = y ~ x, ssn.object = s, estmethod = method, ddf = "asymptotic", tailup_type = "none",
        euclid_initial = euclid_initial("exponential", de = scale^2, range = 1000, known = "given"))
      fixed <- function(theta) {
        args$nugget_initial <- nugget_initial("nugget", exp(theta), known = "given")
        do.call(ssn_lm, args)
      }
      reference <- optimize(function(theta) -2 * as.numeric(logLik(fixed(theta))), log(scale^2 * c(.01, 10)), tol = 1e-7)
      for (start in c(1e-8, .3, 5)) {
        args$nugget_initial <- nugget_initial("nugget", start * scale^2)
        fit <- do.call(ssn_lm, args)
        nugget <- coef(fit, type = "nugget")[["nugget"]]
        ref <- fixed(log(nugget))
        expect_equal(log(nugget), reference$minimum, tolerance = 1e-5)
        expect_equal(fit$optim$value, reference$objective, tolerance = 1e-9)
        expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-10)
        expect_equal(coef(fit), coef(ref), tolerance = 1e-10)
        expect_equal(vcov(fit), vcov(ref), tolerance = 1e-10)
        expect_equal(fitted(fit), fitted(ref), tolerance = 1e-10)
        expect_false(fit$is_known$nugget[["nugget"]])
        if (scale < 1) expect_lt(nugget, 1e-4)
      }
    }
  }
})

test_that("Gaussian nugget estimates remain coherent with network and local covariance", {
  s <- gaussian_nugget_data()
  for (path in c("network", "blocks", "partition", "random", "anisotropic")) for (method in c("ml", "reml")) {
    args <- list(formula = y ~ x + offset(off), ssn.object = s, estmethod = method, ddf = "asymptotic",
      tailup_initial = tailup_initial("exponential", de = .2, range = 1000, known = "given"),
      taildown_initial = taildown_initial("exponential", de = .2, range = 1000, known = "given"),
      additive = "afvArea", nugget_initial = nugget_initial("nugget", .3))
    if (path == "blocks") args$local <- list(index = rep(1:3, each = 12))
    if (path == "partition") args$partition_factor <- ~part
    if (path == "random") {
      args$random <- ~g
      args$randcov_initial <- randcov_initial(g = .1, known = "given")
    }
    if (path == "anisotropic") {
      args$anisotropy <- TRUE
      args$euclid_initial <- euclid_initial("exponential", de = .2, range = 1000, rotate = .3, scale = .7, known = "given")
    }
    fit <- do.call(ssn_lm, args)
    nugget <- coef(fit, type = "nugget")[["nugget"]]
    args$nugget_initial <- nugget_initial("nugget", nugget, known = "given")
    ref <- do.call(ssn_lm, args)
    expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-9)
    expect_equal(coef(fit), coef(ref), tolerance = 1e-9)
    expect_equal(vcov(fit), vcov(ref), tolerance = 1e-9)
  }
})

test_that("Gaussian floor optimums preserve fixed values and unrelated searches", {
  s <- gaussian_nugget_data()
  s$obs$y <- 1 + s$obs$x + .001 * sin(seq_len(nrow(s$obs)))
  for (method in c("ml", "reml")) {
    args <- list(formula = y ~ x, ssn.object = s, estmethod = method, ddf = "asymptotic", tailup_type = "none",
      euclid_initial = euclid_initial("exponential", de = 1, range = 1000, known = "given"),
      nugget_initial = nugget_initial("nugget", .3))
    fit <- do.call(ssn_lm, args)
    nugget <- coef(fit, type = "nugget")[["nugget"]]
    expect_equal(nugget / 1e-4, 1, tolerance = 1e-6)
    args$nugget_initial <- nugget_initial("nugget", nugget, known = "given")
    ref <- do.call(ssn_lm, args)
    expect_equal(fit$optim$value, -2 * as.numeric(logLik(ref)), tolerance = 1e-10)
  }
  s <- gaussian_nugget_data()
  for (path in c("iid", "fixed", "range", "multiple", "zero_floor")) {
    args <- list(formula = y ~ x, ssn.object = s, ddf = "asymptotic", tailup_type = "none")
    if (path != "iid") args$euclid_initial <- euclid_initial("exponential", de = if (path == "zero_floor") 0 else 1,
      range = 1000, known = if (path == "multiple") character(0) else if (path == "range") "de" else "given")
    args$nugget_initial <- nugget_initial("nugget", .3, known = if (path %in% c("fixed", "range")) "given" else character(0))
    fit <- do.call(ssn_lm, args)
    ref <- with_mocked_bindings(do.call(ssn_lm, args),
      optim_with_nugget_bound = function(optim_args, optim_dotlist, lower) do.call("optim", c(optim_args, optim_dotlist)))
    for (field in c("coefficients", "optim", "is_known", "fitted", "vcov")) expect_identical(fit[[field]], ref[[field]])
  }
})
