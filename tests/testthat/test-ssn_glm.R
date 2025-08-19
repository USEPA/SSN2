test_that("generics work ssn_glm point data", {
  set.seed(2)

  form <- Summer_mn ~ ELEV_DEM
  ssn_mod1 <- ssn_glm(form,
    family = "Gamma", mf04p, tailup_type = "exponential",
    taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea"
  )
  ssn_mod2 <- ssn_glm(form,
    family = Gamma, mf04p, tailup_type = "exponential",
    taildown_type = "none", euclid_type = "none",
    nugget_type = "nugget", additive = "afvArea"
  )

  # AIC
  expect_vector(AIC(ssn_mod1))
  expect_equal(AIC(ssn_mod1), 181.112, tolerance = 0.01)
  expect_s3_class(AIC(ssn_mod1, ssn_mod2), "data.frame")
  expect_equal(AIC(ssn_mod1, ssn_mod2)$AIC, c(181.112, 181.888), tolerance = 0.01)

  # AICc
  expect_vector(AICc(ssn_mod1))
  expect_equal(AICc(ssn_mod1), 185.112, tolerance = 0.01)
  expect_s3_class(AICc(ssn_mod1, ssn_mod2), "data.frame")
  expect_equal(AICc(ssn_mod1, ssn_mod2)$AIC, c(185.112, 182.888), tolerance = 0.01)

  # anova
  anova1 <- anova(ssn_mod1)
  expect_s3_class(anova1, "data.frame")
  expect_s3_class(anova1, "anova.ssn_glm")
  expect_s3_class(tidy(anova1), "data.frame")
  expect_equal(tidy(anova1)$statistic, c(44.9, 19.9), tolerance = 0.01)
  anova12 <- anova(ssn_mod1, ssn_mod2)
  expect_s3_class(anova12, "data.frame")
  expect_s3_class(anova12, "anova.ssn_glm")
  expect_s3_class(tidy(anova12), "data.frame")
  expect_equal(tidy(anova12)$statistic, 8.78, tolerance = 0.01)

  # augment
  aug_ssn_mod1 <- augment(ssn_mod1)
  expect_s3_class(aug_ssn_mod1, "sf")
  expect_true(all(c(".fitted", ".resid", ".hat", ".cooksd", ".std.resid") %in% names(aug_ssn_mod1)))
  expect_equal(aug_ssn_mod1$.fitted[1], log(14.918), tolerance = 0.01)
  aug_pred_ssn_mod1 <- augment(ssn_mod1, newdata = "pred1km")
  expect_s3_class(aug_pred_ssn_mod1, "sf")
  expect_true(all(c(".fitted") %in% names(aug_ssn_mod1)))
  expect_equal(aug_pred_ssn_mod1$.fitted[1], 2.687, tolerance = 0.01)
  aug_ssn_mod1 <- augment(ssn_mod1, type.predict = "response")
  expect_s3_class(aug_ssn_mod1, "sf")
  expect_equal(aug_ssn_mod1$.fitted[1], 14.918, tolerance = 0.01)
  aug_pred_ssn_mod1 <- augment(ssn_mod1, newdata = "pred1km", type.predict = "response")
  expect_equal(aug_pred_ssn_mod1$.fitted[1], exp(2.687), tolerance = 0.01)

  # BIC
  expect_vector(BIC(ssn_mod1))
  expect_equal(BIC(ssn_mod1), 195.57, tolerance = 0.01)
  expect_s3_class(BIC(ssn_mod1, ssn_mod2), "data.frame")
  expect_equal(BIC(ssn_mod1, ssn_mod2)$BIC, c(195.57, 189.11), tolerance = 0.01)

  # coef
  expect_vector(coef(ssn_mod1))
  expect_equal(coef(ssn_mod1), c("(Intercept)" = 7.379, "ELEV_DEM" = -0.002), tolerance = 0.01)
  expect_s3_class(coef(ssn_mod1, type = "tailup"), "tailup_exponential")
  expect_equal(unclass(coef(ssn_mod1, type = "tailup")), c("de" = 1.113e-02, "range" = 1.439e07), tolerance = 0.01)
  expect_s3_class(coef(ssn_mod1, type = "taildown"), "taildown_exponential")
  expect_equal(unclass(coef(ssn_mod1, type = "taildown")), c("de" = 1.826e-02, "range" = 7.101e04), tolerance = 0.01)
  expect_s3_class(coef(ssn_mod1, type = "euclid"), "euclid_exponential")
  expect_equal(unclass(coef(ssn_mod1, type = "euclid")), c("de" = 1.158e-07, "range" = 4.94e03, "rotate" = 0, "scale" = 1), tolerance = 0.01)
  expect_type(coef(ssn_mod1, type = "ssn"), "list")
  expect_null(coef(ssn_mod1, type = "randcov"))
  expect_vector(coef(ssn_mod1, type = "dispersion"))
  expect_s3_class(coef(ssn_mod1, type = "dispersion"), "Gamma")
  expect_equal(unclass(coef(ssn_mod1, type = "dispersion")), c("dispersion" = 82810.24), tolerance = 0.01)
  expect_error(coef(ssn_mod1, type = "error"), 'Invalid type argument. The type argument must be "fixed", "ssn", "tailup",  "taildown",  "euclid",  "nugget", "dispersion", or "randcov".')


  # coefficients alias
  expect_vector(coefficients(ssn_mod1))
  expect_equal(coefficients(ssn_mod1), c("(Intercept)" = 7.379, "ELEV_DEM" = -0.002), tolerance = 0.01)
  expect_s3_class(coefficients(ssn_mod1, type = "tailup"), "tailup_exponential")
  expect_equal(unclass(coefficients(ssn_mod1, type = "tailup")), c("de" = 1.113e-02, "range" = 1.439e07), tolerance = 0.01)
  expect_s3_class(coefficients(ssn_mod1, type = "taildown"), "taildown_exponential")
  expect_equal(unclass(coefficients(ssn_mod1, type = "taildown")), c("de" = 1.826e-02, "range" = 7.101e04), tolerance = 0.01)
  expect_s3_class(coefficients(ssn_mod1, type = "euclid"), "euclid_exponential")
  expect_equal(unclass(coefficients(ssn_mod1, type = "euclid")), c("de" = 1.158e-07, "range" = 4.94e03, "rotate" = 0, "scale" = 1), tolerance = 0.01)
  expect_type(coefficients(ssn_mod1, type = "ssn"), "list")
  expect_null(coefficients(ssn_mod1, type = "randcov"))
  expect_vector(coefficients(ssn_mod1, type = "dispersion"))
  expect_s3_class(coefficients(ssn_mod1, type = "dispersion"), "Gamma")
  expect_equal(unclass(coefficients(ssn_mod1, type = "dispersion")), c("dispersion" = 82810.24), tolerance = 0.01)
  expect_error(coefficients(ssn_mod1, type = "error"), 'Invalid type argument. The type argument must be "fixed", "ssn", "tailup",  "taildown",  "euclid",  "nugget", "dispersion", or "randcov".')

  # confint
  expect_equal(dim(confint(ssn_mod1)), c(2, 2))
  expect_equal(dim(confint(ssn_mod1, parm = c("ELEV_DEM"), level = 0.9)), c(1, 2))

  # cooks.distance
  expect_vector(cooks.distance(ssn_mod1))
  expect_equal(cooks.distance(ssn_mod1)[1], c("1" = 4.303e-10), tolerance = 0.01)

  # covmatrix
  expect_equal(dim(covmatrix(ssn_mod1)), c(45, 45))
  expect_equal(dim(covmatrix(ssn_mod1, "pred1km")), c(175, 45))
  expect_equal(dim(covmatrix(ssn_mod1, "pred1km", cov_type = "obs.pred")), c(45, 175))
  expect_equal(dim(covmatrix(ssn_mod1, "pred1km", cov_type = "pred.pred")), c(175, 175))
  expect_error(covmatrix(ssn_mod1, "pred1km", cov_type = "error"), 'Invalid "cov_type" argument.')

  # deviance
  expect_vector(deviance(ssn_mod1))
  expect_equal(deviance(ssn_mod1), 1.599e-05, tolerance = 0.01)

  # fitted
  expect_vector(fitted(ssn_mod1))
  expect_equal(fitted(ssn_mod1)[1], c("1" = 14.918), tolerance = 0.01)
  expect_vector(fitted(ssn_mod1, type = "tailup"))
  expect_equal(fitted(ssn_mod1, type = "tailup")[1], c("1" = -0.01113), tolerance = 0.01)
  expect_vector(fitted(ssn_mod1, type = "taildown"))
  expect_equal(fitted(ssn_mod1, type = "taildown")[1], c("1" = 0.061), tolerance = 0.01)
  expect_vector(fitted(ssn_mod1, type = "euclid"))
  expect_equal(fitted(ssn_mod1, type = "euclid")[1], c("1" = 3.574e-09), tolerance = 0.01)
  expect_null(fitted(ssn_mod1, type = "randcov"))
  expect_error(fitted(ssn_mod1, type = "error"), 'Invalid type argument. The type argument must be "response", "tailup",  "taildown",  "euclid",  "nugget", or "randcov".')

  # fitted.values alias
  expect_vector(fitted.values(ssn_mod1))
  expect_equal(fitted.values(ssn_mod1)[1], c("1" = 14.918), tolerance = 0.01)
  expect_vector(fitted.values(ssn_mod1, type = "tailup"))
  expect_equal(fitted.values(ssn_mod1, type = "tailup")[1], c("1" = -0.01113), tolerance = 0.01)
  expect_vector(fitted.values(ssn_mod1, type = "taildown"))
  expect_equal(fitted.values(ssn_mod1, type = "taildown")[1], c("1" = 0.061), tolerance = 0.01)
  expect_vector(fitted.values(ssn_mod1, type = "euclid"))
  expect_equal(fitted.values(ssn_mod1, type = "euclid")[1], c("1" = 3.574e-09), tolerance = 0.01)
  expect_null(fitted.values(ssn_mod1, type = "randcov"))
  expect_error(fitted.values(ssn_mod1, type = "error"), 'Invalid type argument. The type argument must be "response", "tailup",  "taildown",  "euclid",  "nugget", or "randcov".')


  # formula
  expect_type(formula(ssn_mod1), "language")
  expect_equal(formula(ssn_mod1), form)

  # getCall
  expect_type(getCall(ssn_mod1), "language")
  expect_equal(getCall(ssn_mod1), ssn_mod1$call)

  # glance
  glance_ssn_mod1 <- glance(ssn_mod1)
  names_glance <- c("n", "p", "npar", "value", "AIC", "AICc", "BIC", "logLik", "deviance", "pseudo.r.squared")
  expect_s3_class(glance_ssn_mod1, "data.frame")
  expect_equal(dim(glance_ssn_mod1), c(1, 10))
  expect_identical(names_glance, names(glance_ssn_mod1))

  # glances
  expect_identical(glance_ssn_mod1, glances(ssn_mod1)[, -1])
  glance_ssn_mod12 <- glances(ssn_mod1, ssn_mod2)
  expect_s3_class(glances(ssn_mod1, ssn_mod2), "data.frame")
  expect_equal(dim(glance_ssn_mod12), c(2, 11))
  expect_identical(c("model", names_glance), names(glance_ssn_mod12))

  # hatvalues
  expect_vector(hatvalues(ssn_mod1))
  expect_equal(hatvalues(ssn_mod1)[1], c("1" = 0.06171), tolerance = 0.01)

  # influence
  infl_ssn_mod1 <- influence(ssn_mod1)
  names_infl <- c(".resid", ".hat", ".cooksd", ".std.resid")
  expect_s3_class(infl_ssn_mod1, "data.frame")
  expect_equal(dim(infl_ssn_mod1), c(45, 4))
  expect_identical(names_infl, names(infl_ssn_mod1))

  # labels
  expect_type(labels(ssn_mod1), "character")
  expect_identical(labels(ssn_mod1), "ELEV_DEM")

  # logLik
  expect_vector(as.numeric(logLik(ssn_mod1)))
  expect_equal(as.numeric(logLik(ssn_mod1)), -82.556, tolerance = 0.01)

  # loocv
  loocv_ssn_mod1 <- loocv(ssn_mod1)
  expect_s3_class(loocv_ssn_mod1, "data.frame")
  expect_identical(names(loocv_ssn_mod1), c("bias", "MSPE", "RMSPE", "RAV"))
  expect_equal(loocv_ssn_mod1$MSPE, 0.243, tolerance = 0.01)
  expect_identical(names(loocv(ssn_mod1, cv_predict = TRUE, se.fit = TRUE)), c("stats", "cv_predict", "se.fit"))
  expect_equal(loocv(ssn_mod1, cv_predict = TRUE, se.fit = TRUE)$cv_predict[1], 2.695, tolerance = 0.01)
  expect_equal(loocv(ssn_mod1, cv_predict = TRUE, se.fit = TRUE)$se.fit[1], 0.0285, tolerance = 0.01)
  expect_equal(loocv(ssn_mod1, cv_predict = TRUE, se.fit = TRUE, type = "response")$cv_predict[1], exp(2.695), tolerance = 0.01)
  expect_equal(loocv(ssn_mod1, cv_predict = TRUE, se.fit = TRUE, type = "response")$se.fit[1], 0.0285, tolerance = 0.01) # still on link scale

  # model.frame
  expect_equal(dim(model.frame(ssn_mod1)), c(45, 2))

  # model.matrix
  expect_equal(dim(model.matrix(ssn_mod1)), c(45, 2))

  # model.offset
  expect_null(model.offset(model.frame(ssn_mod1)))

  # model.response
  expect_length(model.response(model.frame(ssn_mod1)), 45)

  # plot
  expect_invisible(plot(ssn_mod1, which = 1))
  expect_invisible(plot(ssn_mod1, which = 2))
  expect_invisible(plot(ssn_mod1, which = 3))
  expect_invisible(plot(ssn_mod1, which = 4))
  expect_invisible(plot(ssn_mod1, which = 5))
  expect_invisible(plot(ssn_mod1, which = 6))

  # predict
  expect_equal(predict(ssn_mod1, newdata = "pred1km")[1],
    c("1" = 2.687),
    tolerance = 0.01
  )
  expect_equal(predict(ssn_mod1, newdata = "pred1km", type = "response")[1],
    c("1" = 14.688),
    tolerance = 0.01
  )
  expect_equal(length(predict(ssn_mod1, newdata = "pred1km")), 175)
  expect_equal(length(predict(ssn_mod1, newdata = "pred1km", type = "response")), 175)
  expect_identical(names(predict(ssn_mod1, newdata = "pred1km", se.fit = TRUE)), c("fit", "se.fit"))
  expect_identical(names(predict(ssn_mod1, newdata = "pred1km", type = "response", se.fit = TRUE)), c("fit", "se.fit"))
  expect_equal(predict(ssn_mod1, newdata = "pred1km", interval = "prediction", se.fit = TRUE)$fit[1, ],
    c("fit" = 2.687, "lwr" = 2.670, "upr" = 2.740),
    tolerance = 0.01
  )
  expect_equal(predict(ssn_mod1, newdata = "pred1km", interval = "prediction", type = "response", se.fit = TRUE)$fit[1, ],
    c("fit" = 14.688, "lwr" = 14.440, "upr" = 14.941),
    tolerance = 0.01
  )
  expect_equal(predict(ssn_mod1, newdata = "pred1km", interval = "confidence")[1, ],
    c("fit" = 2.635, "lwr" = 2.444, "upr" = 2.827),
    tolerance = 0.01
  )
  expect_equal(predict(ssn_mod1, newdata = "pred1km", type = "response", interval = "confidence")[1, ],
    c("fit" = 13.952, "lwr" = 11.519, "upr" = 16.898),
    tolerance = 0.01
  )
  expect_error(predict(ssn_mod1, newdata = "pred1km", interval = "error"))

  # print
  expect_output(print(ssn_mod1))
  expect_output(print(summary(ssn_mod1)))
  expect_output(print(anova(ssn_mod1)))

  # pseudoR2
  expect_vector(pseudoR2(ssn_mod1))
  expect_equal(pseudoR2(ssn_mod1), 0.162, tolerance = 0.01)

  # residuals
  expect_vector(residuals(ssn_mod1))
  expect_equal(length(residuals(ssn_mod1)), 45)
  expect_vector(residuals(ssn_mod1, type = "pearson"))
  expect_equal(residuals(ssn_mod1, type = "pearson")[1], c("1" = 0.0318), tolerance = 0.01)
  expect_equal(residuals(ssn_mod1, type = "standardized")[1], c("1" = 0.0001), tolerance = 0.01)
  expect_identical(residuals(ssn_mod1, type = "standardized"), rstandard(ssn_mod1))
  expect_error(residuals(ssn_mod1, type = "error"), "residuals must be deviance or response or pearson or standardized")
  # resid alias
  expect_vector(resid(ssn_mod1))
  expect_equal(length(resid(ssn_mod1)), 45)
  expect_vector(resid(ssn_mod1, type = "pearson"))
  expect_equal(resid(ssn_mod1, type = "pearson")[1], c("1" = 0.0318), tolerance = 0.01)
  expect_equal(resid(ssn_mod1, type = "standardized")[1], c("1" = 0.0001), tolerance = 0.01)
  expect_identical(resid(ssn_mod1, type = "standardized"), rstandard(ssn_mod1))
  expect_error(resid(ssn_mod1, type = "error"), "residuals must be deviance or response or pearson or standardized")

  # summary
  expect_type(summary(ssn_mod1), "list")
  expect_equal(length(ssn_mod1), 35)
  expect_equal(length(summary(ssn_mod1)), 8)
  expect_identical(ssn_mod1$family, "Gamma")

  # terms
  expect_type(terms(ssn_mod1), "language")
  expect_s3_class(terms(ssn_mod1), "terms")
  expect_s3_class(terms(ssn_mod1), "formula")

  # tidy
  expect_s3_class(tidy(ssn_mod1), "data.frame")
  expect_equal(dim(tidy(ssn_mod1)), c(2, 5))
  expect_s3_class(tidy(ssn_mod1, conf.int = TRUE, level = 0.9), "data.frame")
  expect_equal(dim(tidy(ssn_mod1, conf.int = TRUE, level = 0.9)), c(2, 7))
  expect_s3_class(tidy(ssn_mod1, effects = "ssn"), "data.frame")
  expect_equal(dim(tidy(ssn_mod1, effects = "ssn")), c(8, 4))
  expect_s3_class(tidy(ssn_mod1, effects = "tailup"), "data.frame")
  expect_equal(dim(tidy(ssn_mod1, effects = "tailup")), c(2, 3))
  expect_s3_class(tidy(ssn_mod1, effects = "taildown"), "data.frame")
  expect_equal(dim(tidy(ssn_mod1, effects = "taildown")), c(2, 3))
  expect_s3_class(tidy(ssn_mod1, effects = "euclid"), "data.frame")
  expect_equal(dim(tidy(ssn_mod1, effects = "euclid")), c(2, 3))
  expect_s3_class(tidy(ssn_mod1, effects = "nugget"), "data.frame")
  expect_equal(dim(tidy(ssn_mod1, effects = "nugget")), c(1, 3))
  expect_s3_class(tidy(ssn_mod1, effects = "dispersion"), "data.frame")
  expect_equal(dim(tidy(ssn_mod1, effects = "nugget")), c(1, 3))

  # update
  expect_s3_class(update(ssn_mod2), "ssn_glm")

  # varcomp
  expect_s3_class(varcomp(ssn_mod1), "data.frame")
  expect_equal(dim(varcomp(ssn_mod1)), c(5, 2))
  expect_equal(varcomp(ssn_mod1)$proportion, c(0.163, 0.317, 0.520, 3.3e-08, 4.33e-11), tolerance = 0.01)

  # vcov
  expect_equal(dim(vcov(ssn_mod1)), c(2, 2))
  expect_equal(diag(vcov(ssn_mod1)), c("(Intercept)" = 1.213, "ELEV_DEM" = 2.959e-07), tolerance = 0.01)
  expect_equal(dim(vcov(ssn_mod1, var_correct = FALSE)), c(2, 2))
  expect_equal(diag(vcov(ssn_mod1, var_correct = FALSE)), c("(Intercept)" = 1.194, "ELEV_DEM" = 2.911e-07), tolerance = 0.01)
})
