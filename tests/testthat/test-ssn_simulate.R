test_that("simulating works", {
  tu <- tailup_params("exponential", de = 1, range = 1)
  td <- taildown_params("exponential", de = 1, range = 1)
  eu <- euclid_params("exponential", de = 1, range = 1, rotate = 0, scale = 1)
  eu2 <- euclid_params("exponential", de = 1, range = 1, rotate = pi / 2, scale = 0.5)
  nug <- nugget_params("nugget", nugget = 1)
  rand <- spmodel::randcov_params("netID" = 1)

  # mean seq
  set.seed(0)
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
  expect_equal(sim1[1:2], c(2.525, -0.652), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = Gaussian, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c(-2.066, 3.221), tolerance = 0.01)


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
  set.seed(0)
  sim1 <- ssn_simulate(
    family = "binomial", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(1, 0), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = binomial, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 1, "2" = 1), tolerance = 0.01)

  # ssn_rbeta
  set.seed(0)
  sim1 <- ssn_simulate(
    family = "beta", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(0.9997, 0.7529), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = beta, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 0.0001, "2" = 0.0141), tolerance = 0.01)

  # ssn_rgamma
  set.seed(0)
  sim1 <- ssn_simulate(
    family = "Gamma", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(47.063, 0.427), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = Gamma, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 0.275, "2" = 0.087), tolerance = 0.01)

  # ssn_rinvgauss
  set.seed(0)
  sim1 <- ssn_simulate(
    family = "inverse.gaussian", ssn.object = mf04p, network = "obs",
    tu, td, eu, nug, additive = afvArea, mean = 0, samples = 1
  )

  expect_equal(length(sim1), n_obs)
  expect_equal(sim1[1:2], c(34.695, 0.123), tolerance = 0.01)

  sim2 <- ssn_simulate(
    family = inverse.gaussian, ssn.object = mf04p, network = "obs",
    tu, td, eu2, nug, additive = "afvArea", mean = mean_seq, samples = 2,
    randcov_params = rand, partition_factor = pf
  )

  expect_equal(dim(sim2), c(n_obs, 2))
  expect_equal(sim2[1, ], c("1" = 7.874, "2" = 0.465), tolerance = 0.01)
})
