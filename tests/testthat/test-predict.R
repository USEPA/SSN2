test_that("prediction works", {

  ssn_create_distmat(
    ssn.object = mf04p,
    predpts = c("pred1km"),
    overwrite = TRUE,
    among_predpts = TRUE
  )

  ssn_create_distmat(
    ssn.object = mf04p,
    predpts = c("CapeHorn"),
    overwrite = TRUE,
    only_predpts = TRUE
  )

  ssn_mod1 <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )

  expect_vector(predict(ssn_mod1, "CapeHorn"))
  expect_vector(predict(ssn_mod1, "pred1km", block = TRUE))
  expect_true(inherits(predict(ssn_mod1, "CapeHorn", interval = "confidence"), "matrix"))
  expect_true(inherits(predict(ssn_mod1, "pred1km", interval = "prediction", level = 0.9), "matrix"))
  expect_type(predict(ssn_mod1), type = "list")

  ssn_mod2 <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, "Gamma",
    tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )

  expect_vector(predict(ssn_mod2, "CapeHorn", type = "link"))
  expect_true(inherits(predict(ssn_mod2, "pred1km", interval = "confidence", level = 0.9), "matrix"))
  expect_true(inherits(predict(ssn_mod2, "CapeHorn", interval = "prediction"), "matrix"))
  expect_type(predict(ssn_mod2), type = "list")
})
