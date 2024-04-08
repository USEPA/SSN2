test_that("generics work ssn_lm point data", {

  ssn_mod1 <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )
  ssn_mod2 <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential",
    taildown_type = "none", euclid_type = "none",
    nugget_type = "nugget", additive = "afvArea"
  )

  # AIC
  expect_vector(AIC(ssn_mod1))
  expect_s3_class(AIC(ssn_mod1, ssn_mod2), "data.frame") # turn reml fixed effects warning off

  # anova
  expect_s3_class(anova(ssn_mod1), "data.frame")
  expect_s3_class(anova(ssn_mod1), "anova.ssn_lm")
  expect_s3_class(tidy(anova(ssn_mod1)), "data.frame")
  expect_s3_class(anova(ssn_mod1, ssn_mod2), "data.frame")
  expect_s3_class(anova(ssn_mod1, ssn_mod2), "anova.ssn_lm")
  expect_s3_class(tidy(anova(ssn_mod1, ssn_mod2)), "data.frame")

  # augment
  expect_s3_class(augment(ssn_mod1), "data.frame")
  expect_s3_class(augment(ssn_mod1, newdata = "pred1km"), "data.frame")

  # coef
  expect_vector(coef(ssn_mod1))
  expect_s3_class(coef(ssn_mod1, type = "tailup"), "tailup_exponential")
  expect_s3_class(coef(ssn_mod1, type = "euclid"), "euclid_exponential")
  expect_type(coef(ssn_mod1, type = "ssn"), "list")
  expect_null(coef(ssn_mod1, type = "randcov"))
  expect_vector(coefficients(ssn_mod1))
  expect_s3_class(coefficients(ssn_mod1, type = "taildown"), "taildown_exponential")
  expect_s3_class(coef(ssn_mod1, type = "nugget"), "nugget_nugget")
  expect_null(coefficients(ssn_mod1, type = "randcov"))

  # confint
  expect_true(inherits(confint(ssn_mod1), "matrix"))
  expect_true(inherits(confint(ssn_mod1, parm = c("x"), level = 0.9), "matrix"))

  # cooks.distance
  expect_vector(cooks.distance(ssn_mod1))

  # covmatrix
  expect_true(inherits(covmatrix(ssn_mod1), "matrix"))
  expect_true(inherits(covmatrix(ssn_mod1, "pred1km"), "matrix"))
  expect_true(inherits(covmatrix(ssn_mod1, "pred1km", type = "obs.pred"), "matrix"))
  expect_true(inherits(covmatrix(ssn_mod1, "pred1km", cov_type = "pred.pred"), "matrix"))

  # deviance
  expect_vector(deviance(ssn_mod1))

  # fitted
  expect_vector(fitted(ssn_mod1))
  expect_vector(fitted(ssn_mod1, type = "tailup"))
  expect_vector(fitted(ssn_mod1, type = "taildown"))
  expect_null(fitted(ssn_mod1, type = "randcov"))
  expect_vector(fitted.values(ssn_mod1))
  expect_vector(fitted.values(ssn_mod1, type = "euclid"))
  expect_vector(fitted.values(ssn_mod1, type = "nugget"))
  expect_null(fitted.values(ssn_mod1, type = "randcov"))

  # formula
  expect_type(formula(ssn_mod1), "language")

  # getCall
  expect_type(getCall(ssn_mod1), "language")

  # glance
  expect_s3_class(glance(ssn_mod1), "data.frame")

  # glances
  expect_s3_class(glances(ssn_mod1), "data.frame")
  expect_s3_class(glances(ssn_mod1, ssn_mod2), "data.frame")

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
