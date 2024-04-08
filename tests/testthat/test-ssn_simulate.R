test_that("simulating works", {

  tu <- tailup_params("exponential", de = 1, range = 1)
  td <- taildown_params("exponential", de = 1, range = 1)
  eu <- euclid_params("exponential", de = 1, range = 1, rotate = 0, scale = 1)
  eu2 <- euclid_params("exponential", de = 1, range = 1, rotate = pi / 2, scale = 0.5)
  nug <- nugget_params("nugget", nugget = 1)
  rand <- spmodel::randcov_params("netID" = 1)

  # mean seq
  set.seed(0)
  mean_seq <- rnorm(n = NROW(mf04p$obs), 0, sd = 0.25)

  # set netID as factor
  mf04p$obs$netID <- as.factor(mf04p$obs$netID)

  # partition factor
  pf <- ~netID

  # ssn_rnorm
  expect_vector(ssn_simulate(
    family = "gaussian", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  ))
  expect_true(inherits(ssn_simulate(
    family = Gaussian, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  ), "matrix"))

  # ssn_rpois
  expect_vector(ssn_simulate(
    family = "poisson", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  ))
  expect_true(inherits(ssn_simulate(
    family = poisson, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  ), "matrix"))

  # ssn_rnbinom
  expect_vector(ssn_simulate(
    family = "nbinomial", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1,
    dispersion = 1
  ))
  expect_true(inherits(ssn_simulate(
    family = nbinomial, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    dispersion = 1, randcov_params = rand, partition_factor = pf
  ), "matrix"))

  # ssn_rbinom
  expect_vector(ssn_simulate(
    family = "binomial", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  ))
  expect_true(inherits(ssn_simulate(
    family = binomial, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  ), "matrix"))

  # ssn_rbeta
  expect_vector(ssn_simulate(
    family = "beta", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1,
    dispersion = 1
  ))
  expect_true(inherits(ssn_simulate(
    family = beta, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    dispersion = 1, randcov_params = rand, partition_factor = pf
  ), "matrix"))

  # ssn_rgamma
  expect_vector(ssn_simulate(
    family = "Gamma", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1,
    dispersion = 1
  ))
  expect_true(inherits(ssn_simulate(
    family = Gamma, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    dispersion = 1, randcov_params = rand, partition_factor = pf
  ), "matrix"))

  # ssn_rinvgauss
  expect_vector(ssn_simulate(
    family = "inverse.gaussian", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1,
    dispersion = 1
  ))
  expect_true(inherits(ssn_simulate(
    family = inverse.gaussian, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    dispersion = 1, randcov_params = rand, partition_factor = pf
  ), "matrix"))
})
