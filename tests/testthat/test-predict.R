test_that("prediction works", {
  # Pin legacy starts to keep the historical prediction expectations stable.
  expect_warning(ssn_mod1 <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 500),
    tailup_initial = tailup_initial("exponential", 2.2492, 22085.989, known = "none"),
    taildown_initial = taildown_initial("exponential", 0.0833, 7361.9963, known = "none"),
    euclid_initial = euclid_initial("exponential", 0.0833, 6988.4664, known = "none"),
    nugget_initial = nugget_initial("nugget", 0.0833, known = "none")
  ), "convergence code 1", fixed = TRUE)

  n_p1 <- nrow(mf04p$preds$pred1km)
  n_CH <- nrow(mf04p$preds$CapeHorn)

  preds_CH <- predict(ssn_mod1, "CapeHorn")
  expect_equal(length(preds_CH), n_CH)
  expect_equal(preds_CH[1:2], c("1" = 9.972, "2" = 9.973), tolerance = 0.01)

  preds_p1_block <- predict(ssn_mod1, "pred1km", block = TRUE)
  expect_equal(length(preds_p1_block), 1)
  expect_equal(preds_p1_block, c("1" = 10.295), tolerance = 0.01)

  preds_CH_conf <- predict(ssn_mod1, "CapeHorn", interval = "confidence")
  expect_equal(dim(preds_CH_conf), c(n_CH, 3))
  expect_equal(preds_CH_conf[1, ], c("fit" = 12.411, "lwr" = 9.951, "upr" = 14.871), tolerance = 0.01)

  preds_p1_pred <- predict(ssn_mod1, "pred1km", interval = "prediction", level = 0.9)
  expect_equal(dim(preds_p1_pred), c(n_p1, 3))
  expect_equal(preds_p1_pred[1, ], c("fit" = 14.690, "lwr" = 14.416, "upr" = 14.964), tolerance = 0.01)

  preds_list <- predict(ssn_mod1)
  expect_identical(names(preds_list), c("pred1km", "CapeHorn"))

  ssn_mod2 <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, "Gamma",
    tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea",
    tailup_initial = tailup_initial("exponential", 0.05, 22085.989, known = "none"),
    taildown_initial = taildown_initial("exponential", 0.05, 22085.989, known = "none"),
    euclid_initial = euclid_initial("exponential", 0.05, 20965.3991, known = "none"),
    nugget_initial = nugget_initial("nugget", 0.05, known = "none"),
    dispersion_initial = dispersion_initial("Gamma", 100, known = "none")
  )

  preds_CH <- predict(ssn_mod2, "CapeHorn")
  expect_equal(length(preds_CH), n_CH)
  expect_equal(preds_CH[1:2], c("1" = 2.301, "2" = 2.301), tolerance = 0.01)

  preds_CH <- predict(ssn_mod2, "CapeHorn", type = "response")
  expect_equal(length(preds_CH), n_CH)
  expect_equal(preds_CH[1:2], c("1" = 9.992, "2" = 9.993), tolerance = 0.01)

  preds_CH_conf <- predict(ssn_mod2, "CapeHorn", interval = "confidence")
  expect_equal(dim(preds_CH_conf), c(n_CH, 3))
  expect_equal(preds_CH_conf[1, ], c("fit" = 2.495, "lwr" = 1.597, "upr" = 3.392), tolerance = 0.01)

  preds_p1_pred <- predict(ssn_mod2, "pred1km", interval = "prediction", level = 0.9)
  expect_equal(dim(preds_p1_pred), c(n_p1, 3))
  expect_equal(preds_p1_pred[1, ], c("fit" = 2.687, "lwr" = 2.672, "upr" = 2.701), tolerance = 0.01)

  preds_list <- predict(ssn_mod2)
  expect_identical(names(preds_list), c("pred1km", "CapeHorn"))
})
