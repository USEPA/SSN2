test_that("prediction works", {
  # Copy the mf04p .ssn data to a local directory and read it into R
  # When modeling with your .ssn object, you will load it using the relevant
  # path to the .ssn data on your machine
  copy_lsn_to_temp()
  temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
  mf04p <- ssn_import(
    temp_path,
    predpts = c("pred1km", "CapeHorn", "Knapp"),
    overwrite = TRUE
  )

  ssn_mod1 <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )

  expect_vector(predict(ssn_mod1, "pred1km"))
  expect_vector(predict(ssn_mod1, "pred1km", block = TRUE))
  expect_true(inherits(predict(ssn_mod1, "pred1km", interval = "confidence"), "matrix"))
  expect_true(inherits(predict(ssn_mod1, "pred1km", interval = "prediction", level = 0.9), "matrix"))
  expect_type(predict(ssn_mod1), type = "list")
  expect_type(predict(ssn_mod1, block = TRUE), type = "list")


  ssn_mod2 <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, "Gamma",
    tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )

  expect_vector(predict(ssn_mod2, "pred1km", type = "link"))
  expect_true(inherits(predict(ssn_mod2, "pred1km", interval = "confidence", level = 0.9), "matrix"))
  expect_true(inherits(predict(ssn_mod2, "pred1km", interval = "prediction"), "matrix"))
  expect_type(predict(ssn_mod2), type = "list")
})
