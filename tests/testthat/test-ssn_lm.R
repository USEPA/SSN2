test_that("generics work ssn_lm point data", {

  form <- Summer_mn ~ ELEV_DEM
  ssn_mod1 <- ssn_lm(form, mf04p,
    tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )
  ssn_mod2 <- ssn_lm(form, mf04p,
    tailup_type = "exponential",
    taildown_type = "none", euclid_type = "none",
    nugget_type = "nugget", additive = "afvArea"
  )

  # AIC
  expect_vector(AIC(ssn_mod1))
  expect_equal(AIC(ssn_mod1), 84.934, tolerance = 0.01)
  expect_s3_class(AIC(ssn_mod1, ssn_mod2), "data.frame")
  expect_equal(AIC(ssn_mod1, ssn_mod2)$AIC, c(84.934, 82.893), tolerance = 0.01)

  # anova
  anova1 <- anova(ssn_mod1)
  expect_s3_class(anova1, "data.frame")
  expect_s3_class(anova1, "anova.ssn_lm")
  expect_s3_class(tidy(anova1), "data.frame")
  expect_equal(tidy(anova1)$statistic, c(30.6, 21.1), tolerance = 0.01)
  anova12 <- anova(ssn_mod1, ssn_mod2)
  expect_s3_class(anova12, "data.frame")
  expect_s3_class(anova12, "anova.ssn_lm")
  expect_s3_class(tidy(anova12), "data.frame")
  expect_equal(tidy(anova12)$statistic, 5.96, tolerance = 0.01)

  # augment
  aug_ssn_mod1 <- augment(ssn_mod1)
  expect_s3_class(aug_ssn_mod1, "sf")
  expect_equal(aug_ssn_mod1$.fitted[1], 14.261, tolerance = 0.01)
  aug_pred_ssn_mod1 <- augment(ssn_mod1, newdata = "pred1km")
  expect_s3_class(aug_pred_ssn_mod1, "sf")
  expect_equal(aug_pred_ssn_mod1$.fitted[1], 14.690, tolerance = 0.01)

  # coef
  expect_vector(coef(ssn_mod1))
  expect_equal(coef(ssn_mod1), c("(Intercept)" = 70.54, "ELEV_DEM" = -0.0289), tolerance = 0.01)
  expect_s3_class(coef(ssn_mod1, type = "tailup"), "tailup_exponential")
  expect_equal(unclass(coef(ssn_mod1, type = "tailup")), c("de" = 1.387, "range" = 4.0744e06), tolerance = 0.01)
  expect_s3_class(coef(ssn_mod1, type = "taildown"), "taildown_exponential")
  expect_equal(unclass(coef(ssn_mod1, type = "taildown")), c("de" = 2.419, "range" = 80694), tolerance = 0.01)
  expect_s3_class(coef(ssn_mod1, type = "euclid"), "euclid_exponential")
  expect_equal(unclass(coef(ssn_mod1, type = "euclid")), c("de" = 0.445, "range" = 7.97e05, "rotate" = 0, "scale" = 1), tolerance = 0.01)
  expect_type(coef(ssn_mod1, type = "ssn"), "list")
  expect_null(coef(ssn_mod1, type = "randcov"))
 # coefficients alias
  expect_vector(coefficients(ssn_mod1))
  expect_equal(coefficients(ssn_mod1), c("(Intercept)" = 70.54, "ELEV_DEM" = -0.0289), tolerance = 0.01)
  expect_s3_class(coefficients(ssn_mod1, type = "tailup"), "tailup_exponential")
  expect_equal(unclass(coefficients(ssn_mod1, type = "tailup")), c("de" = 1.387, "range" = 4.0744e06), tolerance = 0.01)
  expect_s3_class(coefficients(ssn_mod1, type = "taildown"), "taildown_exponential")
  expect_equal(unclass(coefficients(ssn_mod1, type = "taildown")), c("de" = 2.419, "range" = 80694), tolerance = 0.01)
  expect_s3_class(coefficients(ssn_mod1, type = "euclid"), "euclid_exponential")
  expect_equal(unclass(coefficients(ssn_mod1, type = "euclid")), c("de" = 0.445, "range" = 7.97e05, "rotate" = 0, "scale" = 1), tolerance = 0.01)
  expect_type(coefficients(ssn_mod1, type = "ssn"), "list")
  expect_null(coefficients(ssn_mod1, type = "randcov"))

  # confint
  expect_equal(dim(confint(ssn_mod1)), c(2, 2))
  expect_equal(dim(confint(ssn_mod1, parm = c("ELEV_DEM"), level = 0.9)), c(1, 2))

  # cooks.distance
  expect_vector(cooks.distance(ssn_mod1))
  expect_equal(cooks.distance(ssn_mod1)[1], c("1" = 0.004), tolerance = 0.01)

  # covmatrix
  expect_equal(dim(covmatrix(ssn_mod1)), c(45, 45))
  expect_equal(dim(covmatrix(ssn_mod1, "pred1km")), c(175, 45))
  expect_equal(dim(covmatrix(ssn_mod1, "pred1km", cov_type = "obs.pred")), c(45, 175))
  expect_equal(dim(covmatrix(ssn_mod1, "pred1km", cov_type = "pred.pred")), c(175, 175))

  # deviance
  expect_vector(deviance(ssn_mod1))
  expect_equal(deviance(ssn_mod1), 41.63, tolerance = 0.01)

  # fitted
  expect_vector(fitted(ssn_mod1))
  expect_equal(fitted(ssn_mod1)[1], c("1" = 14.26), tolerance = 0.01)
  expect_vector(fitted(ssn_mod1, type = "tailup"))
  expect_equal(fitted(ssn_mod1, type = "tailup")[1], c("1" = -0.130), tolerance = 0.01)
  expect_vector(fitted(ssn_mod1, type = "taildown"))
  expect_equal(fitted(ssn_mod1, type = "taildown")[1], c("1" = 0.772), tolerance = 0.01)
  expect_null(fitted(ssn_mod1, type = "randcov"))
  # fitted.values alias
  expect_vector(fitted.values(ssn_mod1))
  expect_equal(fitted.values(ssn_mod1)[1], c("1" = 14.26), tolerance = 0.01)
  expect_vector(fitted.values(ssn_mod1, type = "tailup"))
  expect_equal(fitted.values(ssn_mod1, type = "tailup")[1], c("1" = -0.130), tolerance = 0.01)
  expect_vector(fitted.values(ssn_mod1, type = "taildown"))
  expect_equal(fitted.values(ssn_mod1, type = "taildown")[1], c("1" = 0.772), tolerance = 0.01)
  expect_null(fitted.values(ssn_mod1, type = "randcov"))

  # formula
  expect_type(formula(ssn_mod1), "language")
  expect_equal(formula(ssn_mod1), form)

  # getCall
  expect_type(getCall(ssn_mod1), "language")
  expect_equal(getCall(ssn_mod1), ssn_mod1$call)

  # glance
  glance_ssn_mod1 <- glance(ssn_mod1)
  names_glance <- c("n", "p", "npar", "value", "AIC", "AICc", "logLik", "deviance", "pseudo.r.squared")
  expect_s3_class(glance_ssn_mod1, "data.frame")
  expect_equal(dim(glance_ssn_mod1), c(1, 9))
  expect_identical(names_glance, names(glance_ssn_mod1))

  # glances
  expect_identical(glance_ssn_mod1, glances(ssn_mod1)[, -1])
  glance_ssn_mod12 <- glances(ssn_mod1, ssn_mod2)
  expect_s3_class(glances(ssn_mod1, ssn_mod2), "data.frame")
  expect_equal(dim(glance_ssn_mod12), c(2, 10))
  expect_identical(names_glance, names(glance_ssn_mod12))

  # hatvalues
  expect_vector(hatvalues(ssn_mod1))

  # influence
  expect_s3_class(influence(ssn_mod1), "data.frame")

  # labels
  expect_type(labels(ssn_mod1), "character")

  # logLik
  expect_vector(logLik(ssn_mod1))

  # loocv
  expect_s3_class(loocv(ssn_mod1), "data.frame")
  expect_type(loocv(ssn_mod1, cv_predict = TRUE, se.fit = TRUE), "list")

  # model.frame
  expect_s3_class(model.frame(ssn_mod1), "data.frame")

  # model.matrix
  expect_true(inherits(model.matrix(ssn_mod1), "matrix"))

  # model.offset
  expect_null(model.offset(model.frame(ssn_mod1)))

  # model.response
  expect_vector(model.response(model.frame(ssn_mod1)))

  # plot
  expect_invisible(plot(ssn_mod1, which = 1))
  expect_invisible(plot(ssn_mod1, which = 2))
  expect_invisible(plot(ssn_mod1, which = 3))
  expect_invisible(plot(ssn_mod1, which = 4))
  expect_invisible(plot(ssn_mod1, which = 5))
  expect_invisible(plot(ssn_mod1, which = 6))

  # predict
  expect_vector(predict(ssn_mod1, newdata = "pred1km"))
  expect_type(predict(ssn_mod1, newdata = "pred1km", se.fit = TRUE), "list")
  expect_type(predict(ssn_mod1, newdata = "pred1km", interval = "prediction", se.fit = TRUE), "list")
  expect_true(inherits(predict(ssn_mod1, newdata = "pred1km", interval = "confidence", level = 0.9), "matrix"))
  expect_vector(predict(ssn_mod1, newdata = "pred1km", block = TRUE))
  expect_type(predict(ssn_mod1, newdata = "pred1km", block = TRUE, se.fit = TRUE), "list")
  expect_true(inherits(predict(ssn_mod1, newdata = "pred1km", block = TRUE, interval = "confidence"), "matrix"))
  expect_true(inherits(predict(ssn_mod1, newdata = "pred1km", block = TRUE, interval = "prediction"), "matrix"))


  # print
  expect_output(print(ssn_mod1))
  expect_output(print(summary(ssn_mod1)))
  expect_output(print(anova(ssn_mod1)))

  # pseudoR2
  expect_vector(pseudoR2(ssn_mod1))

  # residuals
  expect_vector(residuals(ssn_mod1))
  expect_vector(residuals(ssn_mod1, type = "pearson"))
  expect_vector(residuals(ssn_mod1, type = "standardized"))
  expect_vector(resid(ssn_mod1))
  expect_vector(resid(ssn_mod1, type = "pearson"))
  expect_vector(resid(ssn_mod1, type = "standardized"))
  expect_vector(rstandard(ssn_mod1))

  # summary
  expect_type(summary(ssn_mod1), "list")

  # terms
  expect_type(terms(ssn_mod1), "language")

  # tidy
  expect_s3_class(tidy(ssn_mod1), "data.frame")
  expect_s3_class(tidy(ssn_mod1, conf.int = TRUE, level = 0.9), "data.frame")
  expect_s3_class(tidy(ssn_mod1, effects = "ssn"), "data.frame")
  expect_s3_class(tidy(ssn_mod1, effects = "tailup"), "data.frame")
  expect_s3_class(tidy(ssn_mod1, effects = "taildown"), "data.frame")
  expect_s3_class(tidy(ssn_mod1, effects = "euclid"), "data.frame")
  expect_s3_class(tidy(ssn_mod1, effects = "nugget"), "data.frame")


  # update
  expect_s3_class(update(ssn_mod2), "ssn_lm")

  # varcomp
  expect_s3_class(varcomp(ssn_mod1), "data.frame")

  # vcov
  expect_true(inherits(vcov(ssn_mod1), "matrix"))
})
