ssn_create_distmat(
  ssn.object = mf04p,
  predpts = c("pred1km"),
  overwrite = TRUE,
  among_predpts = TRUE
)

# fit an example model
ssn_mod <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")

# create an example initial object
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

test_that("dist object output appropriate", {
  # create distance object
  object <- get_dist_object(
    ssn.object = mf04p,
    initial_object = initial_object_val,
    additive = "afvArea",
    anisotropy = FALSE
  )

  # store names and dimensions
  names_vec <- c(
    "distjunc_mat", "mask_mat", "a_mat", "b_mat",
    "hydro_mat", "w_mat", "euclid_mat", "network_index", "pid",
    "dist_order", "inv_dist_order"
  )
  n_obs <- NROW(mf04p$obs)
  n_obs_dim <- c(n_obs, n_obs)

  # run test on object structure
  expect_identical(names(object), names_vec)
  expect_equal(dim(object$distjunc_mat), n_obs_dim)
  expect_equal(dim(object$mask_mat), n_obs_dim)
  expect_equal(dim(object$a_mat), n_obs_dim)
  expect_equal(dim(object$b_mat), n_obs_dim)
  expect_equal(dim(object$hydro_mat), n_obs_dim)
  expect_equal(dim(object$w_mat), n_obs_dim)
  expect_equal(dim(object$euclid_mat), n_obs_dim)
  expect_equal(length(object$network_index), n_obs)
  expect_equal(length(object$pid), n_obs)
  expect_equal(length(object$dist_order), n_obs)
  expect_equal(length(object$inv_dist_order), n_obs)
})

test_that("dist pred object output appropriate", {
  # create distance object
  dist_pred_object <- get_dist_pred_object(
    object = ssn_mod,
    newdata_name = "pred1km",
    initial_object = initial_object_val
  )

  # store names and dimensions
  names_vec <- c(
    "distjunca_pred_mat", "distjuncb_pred_mat", "mask_pred_mat",
    "a_pred_mat", "b_pred_mat",
    "hydro_pred_mat", "w_pred_mat", "euclid_pred_mat",
    "network_index", "pid", "dist_order", "inv_dist_order",
    "network_index_pred", "pid_pred", "dist_order_pred", "inv_dist_order_pred"
  )

  n_obs <- NROW(ssn_mod$ssn.object$obs)
  n_pred <- NROW(ssn_mod$ssn.object$preds[["pred1km"]])
  n_dim <- c(n_pred, n_obs)

  # run test on object structure
  expect_identical(names(dist_pred_object), names_vec)
  expect_equal(dim(dist_pred_object$distjunca_pred_mat), n_dim)
  expect_equal(dim(t(dist_pred_object$distjuncb_pred_mat)), n_dim)
  expect_equal(dim(dist_pred_object$mask_pred_mat), n_dim)
  expect_equal(dim(dist_pred_object$a_pred_mat), n_dim)
  expect_equal(dim(dist_pred_object$b_pred_mat), n_dim)
  expect_equal(dim(dist_pred_object$hydro_pred_mat), n_dim)
  expect_equal(dim(dist_pred_object$w_pred_mat), n_dim)
  expect_equal(dim(dist_pred_object$euclid_pred_mat), n_dim)
  expect_equal(length(dist_pred_object$network_index), n_obs)
  expect_equal(length(dist_pred_object$pid), n_obs)
  expect_equal(length(dist_pred_object$dist_order), n_obs)
  expect_equal(length(dist_pred_object$inv_dist_order), n_obs)
  expect_equal(length(dist_pred_object$network_index_pred), n_pred)
  expect_equal(length(dist_pred_object$pid_pred), n_pred)
  expect_equal(length(dist_pred_object$dist_order_pred), n_pred)
  expect_equal(length(dist_pred_object$inv_dist_order_pred), n_pred)
})

test_that("dist pred bk object output appropriate", {
  # create distance object
  object <- get_dist_predbk_object(
    object = ssn_mod,
    newdata_name = "pred1km",
    initial_object = initial_object_val
  )

  # store names and dimensions
  names_vec <- c(
    "distjunc_mat", "mask_mat", "a_mat", "b_mat",
    "hydro_mat", "w_mat", "euclid_mat", "network_index", "pid",
    "dist_order", "inv_dist_order"
  )
  n_pred <- NROW(mf04p$preds[["pred1km"]])
  n_pred_dim <- c(n_pred, n_pred)

  # run test on object structure
  expect_identical(names(object), names_vec)
  expect_equal(dim(object$distjunc_mat), n_pred_dim)
  expect_equal(dim(object$mask_mat), n_pred_dim)
  expect_equal(dim(object$a_mat), n_pred_dim)
  expect_equal(dim(object$b_mat), n_pred_dim)
  expect_equal(dim(object$hydro_mat), n_pred_dim)
  expect_equal(dim(object$w_mat), n_pred_dim)
  expect_equal(dim(object$euclid_mat), n_pred_dim)
  expect_equal(length(object$network_index), n_pred)
  expect_equal(length(object$pid), n_pred)
  expect_equal(length(object$dist_order), n_pred)
  expect_equal(length(object$inv_dist_order), n_pred)
})
