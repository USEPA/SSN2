test_that("omitting family simulates a Gaussian response", {
  args <- list(
    ssn.object = mf04p, tailup_params = tailup_params("none"),
    taildown_params = taildown_params("none"), euclid_params = euclid_params("none"),
    nugget_params = nugget_params("nugget", nugget = 1), mean = 2, samples = 2
  )
  set.seed(2)
  default <- do.call(ssn_simulate, args)
  set.seed(2)
  explicit <- do.call(ssn_simulate, c(list(family = "Gaussian"), args))
  set.seed(2)
  direct <- do.call(ssn_rnorm, args)
  expect_identical(default, explicit)
  expect_identical(default, direct)
  expect_equal(dim(default), c(nrow(mf04p$obs), 2))
})

test_that("simulating works", {
  tu <- tailup_params("exponential", de = 1, range = 1)
  td <- taildown_params("exponential", de = 1, range = 1)
  eu <- euclid_params("exponential", de = 1, range = 1, rotate = 0, scale = 1)
  eu2 <- euclid_params("exponential", de = 1, range = 1, rotate = pi / 2, scale = 0.5)
  nug <- nugget_params("nugget", nugget = 1)
  rand <- spmodel::randcov_params("netID" = 1)

  # mean seq
  set.seed(2)
  n_obs <- NROW(mf04p$obs)
  mean_seq <- rnorm(n = n_obs, 0, sd = 0.25)

  # set netID as factor
  mf04p$obs$netID <- as.factor(mf04p$obs$netID)

  # partition factor
  pf <- ~netID

  # ssn_rnorm
  set.seed(2)
  sim1 <- ssn_simulate(
    family = "gaussian", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(-1.794, 0.370), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = Gaussian, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c(4.228, 3.354), tolerance = 0.01)


  # ssn_rpois
  set.seed(2)
  sim1 <- ssn_simulate(
    family = "poisson", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(1, 1), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = poisson, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 1, "2" = 13), tolerance = 0.01)

  # ssn_rnbinom
  set.seed(2)
  sim1 <- ssn_simulate(
    family = "binomial", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(1, 1), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = binomial, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 1, "2" = 1), tolerance = 0.01)

  # ssn_rbeta
  set.seed(2)
  sim1 <- ssn_simulate(
    family = "beta", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(0.3150, 0.4132), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = beta, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 0.5265, "2" = 0.1137), tolerance = 0.01)

  # ssn_rgamma
  set.seed(2)
  sim1 <- ssn_simulate(
    family = "Gamma", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(0.4821, 0.4448), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = Gamma, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 0.5403, "2" = 0.02175), tolerance = 0.01)

  # ssn_rinvgauss
  set.seed(2)
  sim1 <- ssn_simulate(
    family = "inverse.gaussian", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(0.9632, 1.0676), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = inverse.gaussian, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 3.2906, "2" = 0.002061), tolerance = 0.01)
})
