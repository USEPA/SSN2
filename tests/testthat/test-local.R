test_that("generics work local big data", {

  ssn_create_bigdist(
    mf04p,
    predpts = c("CapeHorn"),
    overwrite = TRUE,
    among_predpts = TRUE
  )


  set.seed(2)

  form <- Summer_mn ~ ELEV_DEM

  ################## Linear Models

  ######### Fit
  ssn_mod1_bd <- ssn_lm(form, mf04p,
                     tailup_type = "exponential",
                     taildown_type = "exponential", euclid_type = "exponential",
                     nugget_type = "nugget", additive = "afvArea", estmethod = "ml",
                     local = TRUE
  )
  ssn_mod2_bd <- ssn_lm(form, mf04p,
                     tailup_type = "exponential",
                     taildown_type = "none", euclid_type = "none",
                     nugget_type = "nugget", additive = "afvArea",
                     local = list(parallel = TRUE, ncores = 2, size = 20)
  )


  ######### Fit Tests
  expect_equal(coef(ssn_mod1_bd), c("(Intercept)" = 74.66, "ELEV_DEM" = -0.0310), tolerance = 0.01)
  expect_equal(length(unique(ssn_mod1_bd$local_index)), ceiling(NROW(mf04p$obs) / 100))
  expect_equal(coef(ssn_mod2_bd), c("(Intercept)" = 85.64, "ELEV_DEM" = -0.0365), tolerance = 0.01)
  expect_equal(length(unique(ssn_mod2_bd$local_index)), ceiling(NROW(mf04p$obs) / 20))


  ######### Predict/Predict Tests
  n_CH <- nrow(mf04p$preds$CapeHorn)

  preds_p1_pred <- predict(ssn_mod1_bd, "CapeHorn", interval = "prediction", level = 0.9, local = TRUE)
  expect_equal(dim(preds_p1_pred), c(n_CH, 3))
  expect_equal(preds_p1_pred[1, ], c("fit" = 9.97, "lwr" = 9.59, "upr" = 10.34), tolerance = 0.01)

  preds_block1_pred <- predict(ssn_mod1_bd, "CapeHorn", interval = "prediction", level = 0.9, block = TRUE, local = TRUE)
  expect_equal(dim(preds_block1_pred), c(1, 3))
  expect_equal(preds_block1_pred[1, ], c("fit" = 10.02, "lwr" = 9.92, "upr" = 10.12), tolerance = 0.01)

  preds_p2_pred <- predict(ssn_mod2_bd, "CapeHorn", interval = "prediction", level = 0.9, local = list(size = 40, parallel = TRUE, ncores = 2))
  expect_equal(dim(preds_p2_pred), c(n_CH, 3))
  expect_equal(preds_p2_pred[1, ], c("fit" = 9.93, "lwr" = 9.47, "upr" = 10.39), tolerance = 0.01)

  preds_block2_pred <- predict(ssn_mod2_bd, "CapeHorn", interval = "prediction", level = 0.9, block = TRUE, local = list(size = 40, parallel = TRUE, ncores = 2))
  expect_equal(dim(preds_block2_pred), c(1, 3))
  expect_equal(preds_block2_pred[1, ], c("fit" = 10.02, "lwr" = 9.88, "upr" = 10.16), tolerance = 0.01)

  ################## Generalized Linear Models

  ######### Fit
  ssn_mod1_bd <- ssn_glm(form, mf04p,
                     tailup_type = "exponential",
                     taildown_type = "exponential", euclid_type = "exponential",
                     nugget_type = "nugget", additive = "afvArea", estmethod = "ml",
                     local = TRUE, family = "Gamma"
  )
  ssn_mod2_bd <- ssn_glm(form, mf04p,
                     tailup_type = "exponential",
                     taildown_type = "none", euclid_type = "none",
                     nugget_type = "nugget", additive = "afvArea",
                     local = list(parallel = TRUE, ncores = 2, size = 20), family = Gamma
  )


  ######### Fit Tests
  expect_equal(coef(ssn_mod1_bd), c("(Intercept)" = 7.72, "ELEV_DEM" = -0.0026), tolerance = 0.01)
  expect_equal(length(unique(ssn_mod1_bd$local_index)), ceiling(NROW(mf04p$obs) / 100))
  expect_equal(coef(ssn_mod2_bd), c("(Intercept)" = 8.61, "ELEV_DEM" = -0.0030), tolerance = 0.01)
  expect_equal(length(unique(ssn_mod2_bd$local_index)), ceiling(NROW(mf04p$obs) / 20))


  ######### Predict/Predict Tests
  n_CH <- nrow(mf04p$preds$CapeHorn)

  preds_p1_pred <- predict(ssn_mod1_bd, "CapeHorn", interval = "prediction", level = 0.9, local = TRUE)
  expect_equal(dim(preds_p1_pred), c(n_CH, 3))
  expect_equal(preds_p1_pred[1, ], c("fit" = 2.30, "lwr" = 2.28, "upr" = 2.33), tolerance = 0.01)

  # no block prediction (yet)

  preds_p2_pred <- predict(ssn_mod2_bd, "CapeHorn", interval = "prediction", level = 0.9, local = list(size = 40, parallel = TRUE, ncores = 2), type = "response")
  expect_equal(dim(preds_p2_pred), c(n_CH, 3))
  expect_equal(preds_p2_pred[1, ], c("fit" = 9.99, "lwr" = 9.63, "upr" = 10.35), tolerance = 0.01)

  # no block prediction (yet)

})
