skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

simulation_split_network <- function() {
  network <- mf04p
  index <- rev(round(seq(1, NROW(network$obs), length.out = 16)))
  network$obs <- network$obs[index, ]
  network$obs$g <- factor(rep(c("a", "b"), 8))
  network$obs$h <- rep(c("u", "v", "w", "u"), 4)
  network$obs$x <- rep(c(-2, .5, 1, -2), 4)
  network$obs$part <- factor(rep(c("a|b", "a", "a|b", "a"), 4))
  network$obs$other <- rep(c("c", "b|c", "c", "b|c"), 4)
  network$obs$other[16] <- "singleton"
  network$path <- withr::local_tempdir(.local_envir = parent.frame())
  stopifnot(file.copy(file.path(mf04p$path, "binaryID.db"), network$path))
  ssn_create_distmat(network, overwrite = TRUE)
  network
}

simulation_split_params <- function(component = "mixed", effects = TRUE) {
  params <- list(tailup = tailup_params("none"),
    taildown = taildown_params("none"), euclid = euclid_params("none"),
    nugget = nugget_params("nugget", nugget = .2),
    randcov = if (effects) c("1 | g" = .4, "x | g" = .3, "1 | h" = .1) else NULL)
  if (component %in% c("up", "mixed", "wave"))
    params$tailup <- tailup_params("exponential", de = .4, range = 12000)
  if (component %in% c("down", "mixed", "wave"))
    params$taildown <- taildown_params("exponential", de = .3, range = 9000)
  if (component %in% c("euclid", "mixed", "anisotropic", "wave"))
    params$euclid <- euclid_params("exponential", de = .2, range = 16000)
  if (component == "anisotropic")
    params$euclid <- euclid_params("exponential", de = .7, range = 16000,
      rotate = .7, scale = .25)
  if (component %in% c("wave-up", "wave"))
    params$tailup <- tailup_params("wave", de = .4, range = 2500)
  if (component %in% c("wave-down", "wave"))
    params$taildown <- taildown_params("wave", de = .3, range = 3500)
  params
}

simulation_split_fit <- function(network, params, partition = NULL) {
  anisotropy <- !inherits(params$euclid, "euclid_none") &&
    (params$euclid[["rotate"]] != 0 || params$euclid[["scale"]] != 1)
  list(ssn.object = network, coefficients = list(params_object = params),
       additive = "afvArea", anisotropy = anisotropy, partition_factor = partition)
}

simulation_split_native <- function(fit, backend = "dense") {
  params <- fit$coefficients$params_object
  types <- vapply(params[1:4], function(x) SSN2:::remove_covtype(class(x)), character(1))
  initial <- SSN2:::get_initial_object(types[1], types[2], types[3], types[4],
    NULL, NULL, NULL, NULL)
  distance <- SSN2:::get_dist_object(fit$ssn.object, initial, fit$additive,
    fit$anisotropy, backend = backend)
  random <- SSN2:::get_randcov_Zs(fit$ssn.object$obs, names(params$randcov))
  partition <- SSN2:::partition_matrix(fit$partition_factor, fit$ssn.object$obs)
  de <- sum(vapply(params[1:3], function(x) x[["de"]], numeric(1)))
  unname(as.matrix(SSN2:::get_cov_matrix(params, distance, random, partition,
    fit$anisotropy, de, diagtol = 0)))
}

simulation_split_args <- function(fit, samples, mean = 0) {
  params <- fit$coefficients$params_object
  random <- params$randcov
  if (length(random)) names(random) <- sub("^x ", "0 + x ", names(random))
  list(ssn.object = fit$ssn.object, tailup_params = params$tailup,
    taildown_params = params$taildown, euclid_params = params$euclid,
    nugget_params = params$nugget, randcov_params = random,
    partition_factor = fit$partition_factor, additive = fit$additive,
    mean = mean, samples = samples)
}

test_that("dense and filematrix distances preserve simulation identity and row order", {
  skip_if_not_installed("filematrix")
  network <- simulation_split_network()
  fit <- simulation_split_fit(network, simulation_split_params("mixed"), ~part)
  reference <- simulation_split_native(fit)
  set.seed(2)
  dense <- do.call(ssn_rnorm, simulation_split_args(fit, 3L))
  ssn_create_bigdist(network, overwrite = TRUE, no_cores = 1, verbose = FALSE)
  set.seed(2)
  expect_equal(do.call(ssn_rnorm, simulation_split_args(fit, 3L)), dense, tolerance = 1e-12)
  files <- list.files(file.path(network$path, "distance", "obs"),
    pattern = "[.]RData$", full.names = TRUE)
  expect_gt(length(files), 0L)
  expect_true(all(file.remove(files)))
  expect_identical(SSN2:::get_decorrelate_observed_backend(fit), "bigdata")
  expect_equal(simulation_split_native(fit, "bigdata"), reference, tolerance = 1e-12)
  set.seed(2)
  expect_equal(do.call(ssn_rnorm, simulation_split_args(fit, 3L)), dense, tolerance = 1e-12)
})

simulation_split_iid <- function(n = 16L) {
  list(ssn.object = list(obs = data.frame(g = factor(seq_len(n)),
    h = rep(c("a", "b"), length.out = n), x = rep(c(-1, 2), length.out = n),
    part = rep(c("p", "q"), length.out = n), afvArea = rep(1, n))),
    tailup_params = tailup_params("none"), taildown_params = taildown_params("none"),
    euclid_params = euclid_params("none"), nugget_params = nugget_params("nugget", nugget = .4),
    randcov_params = c("1 | g" = .3, "0 + x | h" = .2), partition_factor = ~part,
    additive = "afvArea")
}

test_that("iid effects use sparse designs without covariance or distance setup", {
  unexpected <- function(...) stop("unexpected dense or spatial setup")
  validate <- SSN2:::get_randcov_Zs
  draw <- SSN2:::draw_ssn_simulate_random
  local_mocked_bindings(get_decorrelate_order = unexpected,
    get_dist_object = unexpected, get_dist_object_bigdata_cross = unexpected,
    select_square_dist_backend = unexpected, get_cov_matrix = unexpected,
    partition_matrix = unexpected, chol = unexpected,
    get_randcov_Zs = function(data, randcov_names, ZZt = TRUE, ...) {
      expect_false(ZZt)
      validate(data, randcov_names, ZZt = ZZt, ...)
    }, draw_ssn_simulate_random = function(design, samples) {
      expect_s4_class(design$observed, "sparseMatrix")
      expect_lt(as.numeric(object.size(design$observed)), 1e6)
      draw(design, samples)
    }, .package = "SSN2")
  for (n in c(16L, 6001L)) for (samples in c(1L, 3L)) {
    args <- simulation_split_iid(n)
    design <- SSN2:::get_randcov_design(args$ssn.object$obs, args$randcov_params, ~part)
    q <- NCOL(design$observed)
    for (mean in list(2, seq_len(n) / n)) {
      set.seed(2)
      reference <- matrix(rnorm(n * samples, sd = sqrt(.4)), n, samples) + mean
      coefficients <- matrix(rnorm(q * samples), q, samples) * sqrt(design$variance)
      reference <- reference + as.matrix(design$observed %*% coefficients)
      if (samples == 1) reference <- as.vector(reference)
      seed <- .Random.seed
      set.seed(2)
      actual <- do.call(ssn_rnorm, c(args, list(samples = samples, mean = mean)))
      expect_equal(actual, reference, tolerance = 1e-12)
      expect_identical(.Random.seed, seed)
    }
  }
  for (component in c("tailup", "taildown", "euclid")) {
    args <- simulation_split_iid()
    constructor <- get(paste0(component, "_params"), asNamespace("SSN2"))
    args[[paste0(component, "_params")]] <- constructor("exponential", de = 0, range = 1)
    set.seed(2)
    expect_length(do.call(ssn_rnorm, args), 16L)
  }
})

test_that("empty and zero-variance random designs preserve their draw contracts", {
  observed <- simulation_split_iid()$ssn.object$obs
  empty <- SSN2:::get_randcov_design(observed, NULL)
  set.seed(2)
  seed <- .Random.seed
  expect_null(SSN2:::draw_ssn_simulate_random(empty, 3L))
  expect_identical(.Random.seed, seed)
  design <- SSN2:::get_randcov_design(observed, c("1 | g" = 0, "x | h" = .2), ~part)
  q <- NCOL(design$observed)
  for (samples in c(1L, 3L)) {
    set.seed(2)
    expected <- unname(as.matrix(design$observed %*%
      (matrix(rnorm(q * samples), q, samples) * sqrt(design$variance))))
    seed <- .Random.seed
    set.seed(2)
    actual <- SSN2:::draw_ssn_simulate_random(design, samples)
    expect_identical(actual, expected)
    expect_identical(.Random.seed, seed)
    expect_equal(dim(actual), c(16L, samples))
  }
})

test_that("iid random effects retain group slope partition and network validation", {
  args <- simulation_split_iid()
  expect_error(do.call(ssn_rnorm, c(args, list(network = "pred1km"))), "network must be")
  for (column in c("g", "x", "part")) {
    invalid <- args
    invalid$ssn.object$obs[[column]][1] <- NA
    expect_error(suppressWarnings(do.call(ssn_rnorm, invalid)))
  }
  invalid <- args
  invalid$ssn.object$obs$g <- seq_len(16)
  expect_error(do.call(ssn_rnorm, invalid), "categorical|factor")
  invalid <- args
  invalid$ssn.object$obs$x <- factor(rep(letters[1:3], length.out = 16))
  expect_error(do.call(ssn_rnorm, invalid), "numeric|slope")
  invalid <- args
  invalid$ssn.object$obs$x[1] <- Inf
  expect_error(do.call(ssn_rnorm, invalid), "slope")
  invalid <- args
  invalid$partition_factor <- ~absent
  expect_error(do.call(ssn_rnorm, invalid), "absent")
  invalid <- args
  invalid$tailup_params <- tailup_params("exponential", de = 0, range = 1)
  invalid$additive <- NULL
  expect_error(do.call(ssn_rnorm, invalid), "additive")
  invalid$additive <- "absent"
  expect_error(do.call(ssn_rnorm, invalid), "additive|absent")
  invalid$additive <- "afvArea"
  invalid$ssn.object$obs$afvArea[1] <- Inf
  expect_error(do.call(ssn_rnorm, invalid), "additive|finite")
  for (component in c("tailup", "taildown", "euclid")) {
    constructor <- get(paste0(component, "_params"), asNamespace("SSN2"))
    expect_error(constructor("exponential", de = 0, range = 0), "range")
    invalid <- args
    invalid[[paste0(component, "_params")]] <- constructor("exponential", de = 0, range = 1)
    invalid[[paste0(component, "_params")]][["range"]] <- 0
    expect_error(do.call(ssn_rnorm, invalid), "range")
  }
})

test_that("family wrappers preserve the latent interface and response support", {
  skip_if_not_installed("statmod")
  network <- simulation_split_network()
  fit <- simulation_split_fit(network, simulation_split_params("mixed"), ~part)
  families <- c("gaussian", "poisson", "nbinomial", "binomial", "beta", "Gamma", "inverse.gaussian")
  for (samples in c(1L, 3L)) {
    args <- simulation_split_args(fit, samples, mean = -.5)
    set.seed(2)
    latent <- do.call(ssn_rnorm, args)
    for (family in families) {
      reference <- local({
        set.seed(2)
        do.call(ssn_rnorm, args)
        local_mocked_bindings(ssn_rnorm = function(ssn.object, mean, samples,
          randcov_params, partition_factor, ...) {
          expect_identical(ssn.object, network)
          expect_identical(mean, args$mean)
          expect_identical(samples, args$samples)
          expect_identical(randcov_params, args$randcov_params)
          expect_identical(partition_factor, args$partition_factor)
          latent
        }, .package = "SSN2")
        do.call(ssn_simulate, c(list(family = family, size = 4, dispersion = 2), args))
      })
      set.seed(2)
      actual <- do.call(ssn_simulate, c(list(family = family, size = 4, dispersion = 2), args))
      expect_identical(actual, reference, info = family)
      expect_true(all(is.finite(actual)))
      if (family %in% c("poisson", "nbinomial", "binomial"))
        expect_true(all(actual >= 0 & actual == floor(actual)))
      if (family == "binomial") expect_true(all(actual <= 4))
      if (family == "beta") expect_true(all(actual > 0 & actual < 1))
      if (family %in% c("Gamma", "inverse.gaussian")) expect_true(all(actual > 0))
    }
  }
})
