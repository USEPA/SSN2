test_that("simulating works", {

  tu <- tailup_params("exponential", de = 1, range = 1)
  td <- taildown_params("exponential", de = 1, range = 1)
  eu <- euclid_params("exponential", de = 1, range = 1, rotate = 0, scale = 1)
  eu2 <- euclid_params("exponential", de = 1, range = 1, rotate = pi / 2, scale = 0.5)
  nug <- nugget_params("nugget", nugget = 1)
  rand <- spmodel::randcov_params("netID" = 1)

  # mean seq
  n_obs <- NROW(mf04p$obs)
  mean_seq <- rnorm(n = n_obs, 0, sd = 0.25)

  # set netID as factor
  mf04p$obs$netID <- as.factor(mf04p$obs$netID)

  # partition factor
  pf <- ~netID

  # ssn_rnorm
  set.seed(0)
  sim1 <- ssn_simulate(
    family = "gaussian", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(-2.131, -3.127), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = Gaussian, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c(3.221, -1.821), tolerance = 0.01)


  # ssn_rpois
  set.seed(0)
  sim1 <- ssn_simulate(
    family = "poisson", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(19, 1), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = poisson, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 3, "2" = 6), tolerance = 0.01)

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
