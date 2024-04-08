# test some utility functions

# get an initial value object
initial_object_val <- get_initial_object(
  tailup_type = "exponential",
  taildown_type = "exponential",
  euclid_type = "exponential",
  nugget_type = "nugget",
  tailup_initial = NULL,
  taildown_initial = NULL,
  euclid_initial = NULL,
  nugget_initial = NULL
)

# fit an example model
ssn_mod <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")



################################################################################
############################ check_optim_method
################################################################################

test_that("check optim method works", {
  optim_dotlist <- get_optim_dotlist()

  optim_par <- c(a = 1, b = 2)
  optim_dotlist_val <- check_optim_method(optim_par, optim_dotlist)
  expect_equal(optim_dotlist_val$method, optim_dotlist$method)
  expect_equal(optim_dotlist_val$lower, optim_dotlist$lower)
  expect_equal(optim_dotlist_val$upper, optim_dotlist$upper)

  optim_par <- c(a = 1)
  optim_dotlist_val <- check_optim_method(optim_par, optim_dotlist)
  expect_equal(optim_dotlist_val$method, "Brent")
  expect_equal(optim_dotlist_val$lower, -50)
  expect_equal(optim_dotlist_val$upper, 50)
})

################################################################################
############################ params, cov_matrix, cov_vector work
################################################################################

test_that("params, cov_matrix, cov_vector work", {
  tailup_par <- tailup_params("exponential", 1, 1)
  taildown_par <- taildown_params("exponential", 1, 1)
  euclid_par <- euclid_params("exponential", 1, 1, 0, 1)
  nugget_par <- nugget_params("nugget", 0.1)

  # create dist object
  dist_object <- get_dist_object(
    ssn.object = mf04p,
    initial_object = initial_object_val,
    additive = "afvArea",
    anisotropy = FALSE
  )

  n_obs <- NROW(mf04p$obs)
  n_obs_dim <- c(n_obs, n_obs)

  expect_equal(dim(cov_matrix(tailup_par, dist_object)), n_obs_dim)
  expect_equal(dim(cov_matrix(taildown_par, dist_object)), n_obs_dim)
  expect_equal(dim(cov_matrix(euclid_par, dist_object, anisotropy = FALSE)), n_obs_dim)
  expect_equal(dim(cov_matrix(nugget_par, dist_object, de_scale = 0)), n_obs_dim)


  # create distance object
  dist_pred_object <- get_dist_pred_object(
    object = ssn_mod,
    newdata_name = "pred1km",
    initial_object = initial_object_val
  )

  n_obs <- NROW(ssn_mod$ssn.object$obs)
  n_pred <- NROW(ssn_mod$ssn.object$preds[["pred1km"]])
  n_dim <- c(n_pred, n_obs)

  expect_equal(dim(cov_vector(tailup_par, dist_pred_object)), n_dim)
  expect_equal(dim(cov_vector(taildown_par, dist_pred_object)), n_dim)
  expect_equal(dim(cov_vector(euclid_par, dist_pred_object, anisotropy = FALSE)), n_dim)
})
