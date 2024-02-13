copy_lsn_to_temp()
temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(
  temp_path,
  predpts = c("pred1km"),
  overwrite = TRUE
)

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

  ################
  ##### tailup
  ################

  tailup_covs <- c("linear", "spherical", "exponential",
                   "mariah", "epa", "none")
  lapply(tailup_covs, function(x) {
    if (x == "none") {
      params <- tailup_params(x)
      expect_equal(cov_matrix(params, object), 0)
    } else {
      params <- tailup_params(x, de = 1, range = 1)
      expect_equal(dim(cov_matrix(params, object)), c(45, 45))
    }
  })

  ################
  ##### taildown
  ################

  taildown_covs <- c("linear", "spherical", "exponential",
                   "mariah", "epa", "none")
  lapply(taildown_covs, function(x) {
    if (x == "none") {
      params <- taildown_params(x)
      expect_equal(cov_matrix(params, object), 0)
    } else {
      params <- taildown_params(x, de = 1, range = 1)
      expect_equal(dim(cov_matrix(params, object)), c(45, 45))
    }
  })

  ################
  ##### euclid
  ################

  euclid_covs <- c("exponential", "spherical", "gaussian", "cubic",
                    "pentaspherical", "cosine", "wave", "jbessel", "gravity",
                   "rquad", "magnetic", "none")
  lapply(euclid_covs, function(x) {
    if (x == "none") {
      params <- euclid_params(x)
      expect_equal(cov_matrix(params, object), 0)
    } else if (x == "spherical") {
      params <- euclid_params(x, de = 1, range = 1, rotate = 0.5, scale = 0.5)
      expect_equal(dim(cov_matrix(params, object_anis, anisotropy = TRUE)), c(45, 45))
    } else {
      params <- euclid_params(x, de = 1, range = 1, rotate = 0, scale = 1)
      expect_equal(dim(cov_matrix(params, object, anisotropy = FALSE)), c(45, 45))
    }
  })

  ################
  ##### nugget
  ################

  nugget_covs <- c("nugget", "none")
  lapply(nugget_covs, function(x) {
    if (x == "none") {
      params <- nugget_params(x)
      expect_equal(dim(cov_matrix(params, object, de_scale = 1)), c(45, 45))
    } else {
      params <- nugget_params(x, nugget = 1)
      expect_equal(dim(cov_matrix(params, object, de_scale = 1)), c(45, 45))
    }
  })

})

test_that("covariance vector functions run", {



  ssn_mod <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
                    tailup_type = "exponential", taildown_type = "exponential",
                    euclid_type = "exponential", anisotropy = FALSE,
                    random = ~ as.factor(netID), additive = "afvArea",
                    estmethod = "ml", partition_factor = ~ as.factor(netID))

  ssn_mod_anis <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
                         euclid_type = "exponential", anisotropy = TRUE)


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


  ################
  ##### tailup
  ################

  tailup_covs <- c("linear", "spherical", "exponential",
                   "mariah", "epa", "none")
  lapply(tailup_covs, function(x) {
    if (x == "none") {
      params <- tailup_params(x)
      expect_equal(cov_vector(params, object), 0)
    } else {
      params <- tailup_params(x, de = 1, range = 1)
      expect_equal(dim(cov_vector(params, object)), c(175, 45))
    }
  })

  ################
  ##### taildown
  ################

  taildown_covs <- c("linear", "spherical", "exponential",
                     "mariah", "epa", "none")
  lapply(taildown_covs, function(x) {
    if (x == "none") {
      params <- taildown_params(x)
      expect_equal(cov_vector(params, object), 0)
    } else {
      params <- taildown_params(x, de = 1, range = 1)
      expect_equal(dim(cov_vector(params, object)), c(175, 45))
    }
  })

  ################
  ##### euclid
  ################

  euclid_covs <- c("exponential", "spherical", "gaussian", "cubic",
                   "pentaspherical", "cosine", "wave", "jbessel", "gravity",
                   "rquad", "magnetic", "none")
  lapply(euclid_covs, function(x) {
    if (x == "none") {
      params <- euclid_params(x)
      expect_equal(cov_vector(params, object), 0)
    } else if (x == "spherical") {
      params <- euclid_params(x, de = 1, range = 1, rotate = 0.5, scale = 0.5)
      expect_equal(dim(cov_vector(params, object_anis, anisotropy = TRUE)), c(175, 45))
    } else {
      params <- euclid_params(x, de = 1, range = 1, rotate = 0, scale = 1)
      expect_equal(dim(cov_vector(params, object, anisotropy = FALSE)), c(175, 45))
    }
  })

  ################
  ##### prediction
  ################

  expect_no_error(predict(ssn_mod, "pred1km", interval = "prediction", level = 0.9))
  expect_no_error(predict(ssn_mod_anis, "pred1km"))



})


test_that("initial objects", {
  tu <- tailup_initial("exponential", de = 1, range = 1, known = "given")
  td <- taildown_initial("exponential", de = 1, range = 1, known = c("de", "range"))
  eu <- euclid_initial("exponential", de = 1, range = 1,
                       rotate = 0, scale = 1, known = "given")
  nu <- nugget_initial("nugget", 1, known = "nugget")
  disp <- dispersion_initial("Gamma", dispersion = 1, known = "given")
  ssn_mod <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
                    tailup_initial = tu, taildown_initial = td,
                    euclid_initial = eu, nugget_initial = nu,
                    additive = "afvArea")
  expect_s3_class(ssn_mod, "ssn_lm")

  ssn_mod <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p,
                    tailup_initial = tu, taildown_initial = td,
                    euclid_initial = eu, nugget_initial = nu,
                    dispersion_initial = disp,
                    additive = "afvArea")
  expect_s3_class(ssn_mod, "ssn_glm")
})


test_that("print an ssn", {
  expect_output(print(mf04p))
  expect_output(print(summary(mf04p)))
  expect_output(print(names(mf04p)))
  expect_output(print(ssn_names(mf04p)))
})



test_that("extra test fits", {
  ssn_mod <- ssn_glm(Summer_mn > 11 ~ ELEV_DEM, mf04p, family = "binomial",
                     tailup_type = "exponential", additive = "afvArea",
                     estmethod = "ml")
  expect_s3_class(ssn_mod, "ssn_glm")

  ssn_mod <- ssn_glm(round(Summer_mn) ~ ELEV_DEM, mf04p, family = "poisson",
                     taildown_type = "exponential")
  expect_s3_class(ssn_mod, "ssn_glm")

  ssn_mod <- ssn_glm(round(Summer_mn) ~ ELEV_DEM, mf04p, family = "nbinomial",
                     euclid_type = "exponential")
  expect_s3_class(ssn_mod, "ssn_glm")

  ssn_mod <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "inverse.gaussian",
                     tailup_type = "exponential", additive = "afvArea")
  expect_s3_class(ssn_mod, "ssn_glm")

  ssn_mod <- ssn_glm(ratio ~ ELEV_DEM, mf04p, family = "beta",
                     taildown_type = "exponential")
  expect_s3_class(ssn_mod, "ssn_glm")

  ssn_mod <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, family = "binomial",
                     tailup_type = "exponential", additive = "afvArea",
                     random = ~ as.factor(netID))
  expect_output(print(ssn_mod))
  expect_output(print(summary(ssn_mod)))
  expect_type(fitted(ssn_mod, type = "randcov"), "list")

  ssn_mod <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma",
                    random = ~ as.factor(netID))
  expect_output(print(ssn_mod))
  expect_output(print(summary(ssn_mod)))
  expect_type(fitted(ssn_mod, type = "randcov"), "list")

})
