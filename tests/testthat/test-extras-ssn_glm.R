#
#
# test_local <- FALSE # FALSE for CRAN
#
# if (test_local) {
#
#   ssn_create_distmat(
#     ssn.object = mf04p,
#     predpts = c("pred1km", "CapeHorn"),
#     overwrite = TRUE
#   )
#
#   # set a seed
#   set.seed(2)
#
#   test_that("ssn_glm models fit Gaussian", {
#     ssn_mod <- ssn_glm(Summer_mn ~ ELEV_DEM,
#       family = "gaussian", mf04p, tailup_type = "exponential",
#       taildown_type = "exponential", euclid_type = "exponential",
#       nugget_type = "nugget", additive = "afvArea"
#     )
#     expect_s3_class(ssn_mod, "ssn_lm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("ssn_glm models fit poisson", {
#     ssn_mod <- ssn_glm(round(Summer_mn) ~ ELEV_DEM,
#       family = "poisson", mf04p, tailup_type = "exponential",
#       taildown_type = "exponential", euclid_type = "exponential",
#       nugget_type = "nugget", additive = "afvArea"
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("ssn_glm models fit negative binomial", {
#     ssn_mod <- ssn_glm(round(Summer_mn) ~ ELEV_DEM,
#       family = "nbinomial", mf04p, tailup_type = "exponential",
#       nugget_type = "nugget", additive = "afvArea"
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("ssn_glm models fit binomial", {
#     ssn_mod <- ssn_glm(Summer_mn > 11 ~ ELEV_DEM,
#       family = "binomial", mf04p, tailup_type = "exponential",
#       taildown_type = "exponential", euclid_type = "exponential",
#       nugget_type = "nugget", additive = "afvArea", estmethod = "ml"
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("ssn_glm models fit beta", {
#     mf04p$obs$betavar <- runif(NROW(mf04p$obs), min = 0.25, max = 0.75)
#     ssn_mod <- ssn_glm(betavar ~ ELEV_DEM,
#       family = "beta", mf04p,
#       taildown_type = "exponential", euclid_type = "exponential",
#       nugget_type = "nugget"
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("ssn_glm models fit gamma", {
#     ssn_mod <- ssn_glm(Summer_mn ~ ELEV_DEM,
#       family = "Gamma", mf04p, tailup_type = "exponential",
#       nugget_type = "none", additive = "afvArea"
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("ssn_glm models fit inverse gaussian", {
#     ssn_mod <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, inverse.gaussian,
#      euclid_type = "exponential",
#       nugget_type = "nugget",
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#
#   test_that("random effects work", {
#     ssn_mod <- ssn_glm(Summer_mn > 11 ~ ELEV_DEM, mf04p,
#       family = "binomial", tailup_type = "exponential",
#       taildown_type = "exponential",
#       nugget_type = "nugget", additive = "afvArea",
#       random = ~ as.factor(netID)
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("partition factors work", {
#     ssn_mod <- ssn_glm(Summer_mn > 11 ~ ELEV_DEM, mf04p,
#       family = "binomial",
#       taildown_type = "exponential",
#       partition_factor = ~ as.factor(netID)
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("anisotropy works", {
#     ssn_mod <- ssn_glm(Summer_mn > 11 ~ ELEV_DEM, mf04p,
#       family = "binomial", tailup_type = "exponential",
#       euclid_type = "exponential",
#       nugget_type = "nugget", additive = "afvArea",
#       anisotropy = TRUE
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#   })
#
#   test_that("fixing parameters works", {
#     tu <- tailup_initial("exponential", de = 1, known = "de")
#     ssn_mod <- ssn_glm(Summer_mn > 11 ~ ELEV_DEM, mf04p,
#       family = "binomial", tailup_initial = tu,
#       taildown_type = "exponential",
#       nugget_type = "nugget", additive = "afvArea"
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#     expect_equal(coef(ssn_mod, type = "tailup")[["de"]], 1)
#   })
#
#   test_that("missing data works", {
#     mf04p$obs$Summer_mn[1] <- NA
#     ssn_mod <- ssn_glm(Summer_mn > 11 ~ ELEV_DEM, mf04p,
#       family = "binomial",
#       taildown_type = "exponential",
#       nugget_type = "nugget"
#     )
#     expect_s3_class(ssn_mod, "ssn_glm")
#     expect_vector(predict(ssn_mod, "pred1km"))
#     expect_vector(predict(ssn_mod, ".missing"))
#   })
#
#   # previously from test-extras
#   test_that("extra test fits", {
#     ssn_mod <- ssn_glm(Summer_mn > 11 ~ ELEV_DEM, mf04p, family = "binomial",
#                        tailup_type = "exponential", additive = "afvArea",
#                        estmethod = "ml")
#     expect_s3_class(ssn_mod, "ssn_glm")
#
#     ssn_mod <- ssn_glm(round(Summer_mn) ~ ELEV_DEM, mf04p, family = "poisson",
#                        taildown_type = "exponential")
#     expect_s3_class(ssn_mod, "ssn_glm")
#
#     ssn_mod <- ssn_glm(round(Summer_mn) ~ ELEV_DEM, mf04p, family = "nbinomial",
#                        euclid_type = "exponential")
#     expect_s3_class(ssn_mod, "ssn_glm")
#
#     ssn_mod <- ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "inverse.gaussian",
#                        tailup_type = "exponential", additive = "afvArea")
#     expect_s3_class(ssn_mod, "ssn_glm")
#
#     ssn_mod <- ssn_glm(ratio ~ ELEV_DEM, mf04p, family = "beta",
#                        taildown_type = "exponential")
#     expect_s3_class(ssn_mod, "ssn_glm")
#   })
# }
