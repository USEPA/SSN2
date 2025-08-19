test_that("covariance matrix functions run", {
  ################
  ##### create distance object
  ################

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

  object <- get_dist_object(
    ssn.object = mf04p,
    initial_object = initial_object_val,
    additive = "afvArea",
    anisotropy = FALSE
  )

  object_anis <- get_dist_object(
    ssn.object = mf04p,
    initial_object = initial_object_val,
    additive = "afvArea",
    anisotropy = TRUE
  )

  n_obs <- NROW(mf04p$obs)
  n_obs_dim <- c(n_obs, n_obs)

  ################
  ##### tailup
  ################

  tailup_covs <- c(
    "linear", "spherical", "exponential",
    "mariah", "epa", "gaussian", "none"
  )
  lapply(tailup_covs, function(x) {
    if (x == "none") {
      params <- tailup_params(x)
      covmx <- cov_matrix(params, object)
      expect_equal(covmx, 0)
    } else {
      de <- 1
      params <- tailup_params(x, de = de, range = 1)
      covmx <- cov_matrix(params, object)
      expect_equal(dim(covmx), n_obs_dim)
      expect_equal(unique(diag(covmx)), de) # de value is overall variance
    }
  })

  ################
  ##### taildown
  ################

  taildown_covs <- c(
    "linear", "spherical", "exponential",
    "mariah", "epa", "gaussian", "none"
  )
  lapply(taildown_covs, function(x) {
    if (x == "none") {
      params <- taildown_params(x)
      covmx <- cov_matrix(params, object)
      expect_equal(covmx, 0)
    } else {
      de <- 1
      params <- taildown_params(x, de = de, range = 1)
      covmx <- cov_matrix(params, object)
      expect_equal(dim(covmx), n_obs_dim)
      expect_equal(unique(diag(covmx)), de)
    }
  })

  ################
  ##### euclid
  ################

  euclid_covs <- c(
    "exponential", "spherical", "gaussian", "cubic",
    "pentaspherical", "cosine", "wave", "jbessel", "gravity",
    "rquad", "magnetic", "none"
  )
  lapply(euclid_covs, function(x) {
    if (x == "none") {
      params <- euclid_params(x)
      covmx <- cov_matrix(params, object)
      expect_equal(covmx, 0)
    } else if (x == "spherical") {
      de <- 1
      params <- euclid_params(x, de = de, range = 1, rotate = 0.5, scale = 0.5)
      covmx <- cov_matrix(params, object_anis, anisotropy = TRUE)
      expect_equal(dim(covmx), n_obs_dim)
      expect_equal(unique(diag(covmx)), de)
    } else {
      de <- 1
      params <- euclid_params(x, de = 1, range = 1, rotate = 0, scale = 1)
      covmx <- cov_matrix(params, object, anisotropy = FALSE)
      expect_equal(dim(covmx), n_obs_dim)
      expect_equal(unique(diag(covmx)), de)
    }
  })

  ################
  ##### nugget
  ################

  nugget_covs <- c("nugget", "none")
  lapply(nugget_covs, function(x) {
    if (x == "none") {
      params <- nugget_params(x)
      covmx <- cov_matrix(params, object, de_scale = 1)
      expect_equal(dim(covmx), n_obs_dim)
    } else {
      params <- nugget_params(x, nugget = 1)
      covmx <- cov_matrix(params, object, de_scale = 1)
      expect_equal(dim(covmx), n_obs_dim)
    }
  })
})

test_that("covariance vector functions run", {
  ssn_mod <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", taildown_type = "exponential",
    euclid_type = "exponential", anisotropy = FALSE,
    random = ~ as.factor(netID), additive = "afvArea",
    estmethod = "ml", partition_factor = ~ as.factor(netID)
  )

  ssn_mod_anis <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    euclid_type = "exponential", anisotropy = TRUE
  )


  ################
  ##### create distance object
  ################

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

  object <- get_dist_pred_object(
    object = ssn_mod,
    newdata_name = "pred1km",
    initial_object = initial_object_val
  )

  initial_object_val <- get_initial_object(
    tailup_type = "none",
    taildown_type = "none",
    euclid_type = "exponential",
    nugget_type = "nugget",
    tailup_initial = NULL,
    taildown_initial = NULL,
    euclid_initial = NULL,
    nugget_initial = NULL
  )

  object_anis <- get_dist_pred_object(
    object = ssn_mod_anis,
    newdata_name = "pred1km",
    initial_object = initial_object_val
  )

  n_obs <- NROW(ssn_mod$ssn.object$obs)
  n_pred <- NROW(ssn_mod$ssn.object$preds[["pred1km"]])
  n_dim <- c(n_pred, n_obs)


  ################
  ##### tailup
  ################

  tailup_covs <- c(
    "linear", "spherical", "exponential",
    "mariah", "epa",  "gaussian", "none"
  )
  lapply(tailup_covs, function(x) {
    if (x == "none") {
      params <- tailup_params(x)
      covv <- cov_vector(params, object)
      expect_equal(covv, 0)
    } else {
      de <- 1
      params <- tailup_params(x, de = de, range = 1)
      covv <- cov_vector(params, object)
      expect_equal(dim(covv), n_dim)
    }
  })

  ################
  ##### taildown
  ################

  taildown_covs <- c(
    "linear", "spherical", "exponential",
    "mariah", "epa",  "gaussian", "none"
  )
  lapply(taildown_covs, function(x) {
    if (x == "none") {
      params <- taildown_params(x)
      covv <- cov_vector(params, object)
      expect_equal(covv, 0)
    } else {
      de <- 1
      params <- taildown_params(x, de = de, range = 1)
      covv <- cov_vector(params, object)
      expect_equal(dim(covv), n_dim)
    }
  })

  ################
  ##### euclid
  ################

  euclid_covs <- c(
    "exponential", "spherical", "gaussian", "cubic",
    "pentaspherical", "cosine", "wave", "jbessel", "gravity",
    "rquad", "magnetic", "none"
  )
  lapply(euclid_covs, function(x) {
    if (x == "none") {
      params <- euclid_params(x)
      covv <- cov_vector(params, object)
      expect_equal(covv, 0)
    } else if (x == "spherical") {
      de <- 1
      params <- euclid_params(x, de = 1, range = 1, rotate = 0.5, scale = 0.5)
      covv <- cov_vector(params, object_anis, anisotropy = TRUE)
      expect_equal(dim(covv), n_dim)
    } else {
      de <- 1
      params <- euclid_params(x, de = 1, range = 1, rotate = 0, scale = 1)
      covv <- cov_vector(params, object, anisotropy = FALSE)
      expect_equal(dim(covv), n_dim)
    }
  })

  ################
  ##### prediction
  ################

  preds <- predict(ssn_mod, "pred1km", interval = "prediction", level = 0.9)
  expect_equal(dim(preds), c(n_pred, 3))
  expect_identical(colnames(preds), c("fit", "lwr", "upr"))
  preds_anis <- predict(ssn_mod_anis, "pred1km")
  expect_equal(length(preds_anis), n_pred)
})


test_that("initial objects", {
  de <- 1
  range <- 1
  nugget <- 1
  dispersion <- 1
  tu <- tailup_initial("exponential", de = de, range = range, known = "given")
  td <- taildown_initial("exponential", de = de, range = range, known = c("de", "range"))
  eu <- euclid_initial("exponential",
    de = de, range = range,
    rotate = 0, scale = 1, known = "given"
  )
  nu <- nugget_initial("nugget", nugget, known = "nugget")
  disp <- dispersion_initial("Gamma", dispersion = dispersion, known = "given")
  ssn_mod <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_initial = tu, taildown_initial = td,
    euclid_initial = eu, nugget_initial = nu,
    additive = "afvArea"
  )
  expect_equal(as.vector(coef(ssn_mod, "tailup")), c(de, range))
  expect_equal(as.vector(coef(ssn_mod, "taildown")), c(de, range))
  expect_equal(as.vector(coef(ssn_mod, "euclid")), c(de, range, 0, 1))
  expect_equal(as.vector(coef(ssn_mod, "nugget")), c(de))

  ssn_mod <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_initial = tu, taildown_initial = td,
    euclid_initial = eu, nugget_initial = nu,
    dispersion_initial = disp,
    additive = "afvArea"
  )
  expect_equal(as.vector(coef(ssn_mod, "tailup")), c(de, range))
  expect_equal(as.vector(coef(ssn_mod, "taildown")), c(de, range))
  expect_equal(as.vector(coef(ssn_mod, "euclid")), c(de, range, 0, 1))
  expect_equal(as.vector(coef(ssn_mod, "nugget")), c(de))
  expect_equal(as.vector(coef(ssn_mod, "dispersion")), c(dispersion))
})


test_that("print an ssn", {
  expect_output(print(mf04p))
  expect_output(print(summary(mf04p)))
  expect_output(print(names(mf04p)))
  expect_output(print(ssn_names(mf04p)))
})
