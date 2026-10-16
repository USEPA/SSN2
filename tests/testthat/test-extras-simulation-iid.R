skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

simulation_iid_args <- function(n = 12L, samples = 3L) {
  list(ssn.object = list(obs = data.frame(group = rep(c("a", "b"), length.out = n))),
    tailup_params = tailup_params("none"), taildown_params = taildown_params("none"),
    euclid_params = euclid_params("none"), nugget_params = nugget_params("nugget", nugget = 0.4),
    mean = seq_len(n) / n, samples = samples)
}

test_that("independent simulation avoids spatial setup and matches direct normal draws", {
  unexpected <- function(...) stop("unexpected spatial setup")
  local_mocked_bindings(get_ssn_simulate_local = unexpected,
    get_ssn_simulate_covariance_fit = unexpected, get_dist_object = unexpected,
    select_square_dist_backend = unexpected, get_cov_matrix = unexpected,
    partition_matrix = unexpected, .package = "SSN2")
  for (n in c(1L, 12L, 4999L, 5000L, 5001L)) {
    for (samples in c(1L, 3L)) {
      args <- simulation_iid_args(n, samples)
      args$partition_factor <- ~group
      set.seed(2)
      reference <- vapply(seq_len(samples), function(i) {
        args$mean + rnorm(n, sd = sqrt(0.4))
      }, numeric(n))
      if (samples == 1) reference <- as.vector(reference)
      reference_seed <- .Random.seed
      for (local in list(NULL, FALSE, TRUE, list(), list(approximation = "low-rank"),
                        list(approximation = "vecchia", parallel = TRUE),
                        list(size = -1, ordering = "invalid"), "unused")) {
        set.seed(2)
        expect_message(value <- do.call(ssn_rnorm, c(args, list(local = local))), NA)
        expect_identical(value, reference)
        expect_identical(.Random.seed, reference_seed)
      }
      args$mean <- 2
      set.seed(2)
      scalar <- do.call(ssn_rnorm, args)
      expect_equal(scalar, reference - seq_len(n) / n + 2)
    }
  }
})

test_that("independent simulation retains partition and network validation", {
  args <- simulation_iid_args()
  expect_error(do.call(ssn_rnorm, c(args, list(partition_factor = ~unknown_group))), "unknown_group")
  args$ssn.object$obs$group[1] <- NA
  expect_error(do.call(ssn_rnorm, c(args, list(partition_factor = ~group))), "Missing values")
  expect_error(do.call(ssn_rnorm, c(args, list(network = "pred1km"))), "network must be")
})

test_that("empty random effects and multiple partition factors preserve independent draws", {
  args <- simulation_iid_args()
  args$ssn.object$obs$other <- rep(c("x", "y", "z"), 4)
  set.seed(2)
  reference <- do.call(ssn_rnorm, args)
  for (random in list(NULL, numeric())) {
    set.seed(2)
    value <- do.call(ssn_rnorm, c(args,
      list(randcov_params = random, partition_factor = ~group + other)))
    expect_identical(value, reference)
  }
})

test_that("zero nugget and all-none unconditional simulation retain their existing behavior", {
  args <- simulation_iid_args(nrow(mf04p$obs))
  args$ssn.object <- mf04p
  args$local <- FALSE
  args$nugget_params <- nugget_params("nugget", nugget = 0)
  set.seed(2)
  value <- do.call(ssn_rnorm, args)
  expect_identical(value, matrix(args$mean, nrow(mf04p$obs), args$samples))
  args$nugget_params <- nugget_params("none")
  expect_error(do.call(ssn_rnorm, args), "At least one covariance type")
})

test_that("independent simulation excludes spatial and random effects and degenerate variance", {
  local_mocked_bindings(get_ssn_simulate_local = function(...) stop("general simulation"),
                        .package = "SSN2")
  for (component in c("tailup", "taildown", "euclid")) {
    args <- simulation_iid_args()
    args$additive <- "group"
    constructor <- get(paste0(component, "_params"), asNamespace("SSN2"))
    args[[paste0(component, "_params")]] <- constructor("exponential", de = 0, range = 1)
    expect_error(do.call(ssn_rnorm, args), "general simulation")
  }
  for (variance in c(0, 1)) {
    args <- simulation_iid_args()
    args$randcov_params <- spmodel::randcov_params(group = variance)
    expect_error(do.call(ssn_rnorm, args), "general simulation")
  }
  for (nugget in list(nugget_params("none"), nugget_params("nugget", nugget = 0))) {
    args <- simulation_iid_args()
    args$nugget_params <- nugget
    expect_error(do.call(ssn_rnorm, args), "general simulation")
  }
})

test_that("independent family generators retain the latent normal and response transformations", {
  skip_if_not_installed("statmod")
  families <- c("gaussian", "poisson", "nbinomial", "binomial", "beta", "Gamma", "inverse.gaussian")
  for (samples in c(1L, 3L)) {
    args <- simulation_iid_args(samples = samples)
    reference <- local({
      local_mocked_bindings(ssn_rnorm = function(ssn.object, mean, samples, ...) {
        n <- NROW(ssn.object$obs)
        value <- vapply(seq_len(samples), function(i) {
          mean + as.numeric(diag(sqrt(0.4), n) %*% rnorm(n))
        }, numeric(n))
        if (samples == 1) value <- as.vector(value)
        value
      }, .package = "SSN2")
      lapply(families, function(family) {
        set.seed(2)
        do.call(ssn_simulate, c(list(family = family, size = 4, dispersion = 2), args))
      })
    })
    for (i in seq_along(families)) {
      for (local in list(FALSE, TRUE, list(approximation = "low-rank", size_base = 2))) {
        set.seed(2)
        value <- do.call(ssn_simulate,
          c(list(family = families[i], size = 4, dispersion = 2, local = local), args))
        expect_identical(value, reference[[i]], info = families[i])
      }
    }
  }
})

test_that("independent Gaussian simulation matches spmodel", {
  args <- simulation_iid_args()
  set.seed(2)
  reference <- spmodel::sprnorm(spmodel::spcov_params("none", ie = 0.4),
    mean = args$mean, samples = args$samples, data = args$ssn.object$obs,
    partition_factor = ~group)
  set.seed(2)
  value <- do.call(ssn_rnorm, c(args, list(partition_factor = ~group)))
  expect_identical(value, reference)
})
