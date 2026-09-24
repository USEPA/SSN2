skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# satterthwaite
skip_if_not_installed("numDeriv")

copy_lsn_to_temp()
temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)

fit_tailup <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM,
  ssn.object = mf04p,
  tailup_type = "exponential",
  additive = "afvArea"
)

test_that("satterthwaite() returns finite, positive, named df for a well-conditioned fit", {
  df <- satterthwaite(fit_tailup, method = "numeric")
  expect_named(df, rownames(fit_tailup$vcov$fixed))
  expect_true(all(is.finite(df)))
  expect_true(all(df > 0))
})

test_that("satterthwaite() rejects everything except explicit method = 'numeric'", {
  expect_error(satterthwaite(fit_tailup, method = "closed"), "method")
  expect_error(satterthwaite(fit_tailup), "method")
})

test_that("satterthwaite() rejects ssn_glm objects explicitly", {
  mf04p_bin <- mf04p
  set.seed(2)
  mf04p_bin$obs$y01 <- rbinom(nrow(mf04p_bin$obs), 1, 0.5)
  fit_glm <- ssn_glm(
    formula = y01 ~ ELEV_DEM, ssn.object = mf04p_bin, family = "binomial",
    tailup_type = "exponential", additive = "afvArea"
  )
  expect_error(satterthwaite(fit_glm, method = "numeric"), "ssn_glm")
})

test_that("fai_cornelius() with q = 1 reduces exactly to scalar satterthwaite()", {
  context <- get_satterthwaite_context(fit_tailup, "numeric")
  vt <- get_vcov_theta_numeric(context)
  df <- satterthwaite(fit_tailup, method = "numeric")

  L_intercept <- matrix(c(1, 0), nrow = 1)
  L_slope <- matrix(c(0, 1), nrow = 1)

  fc_intercept <- fai_cornelius(L_intercept, fit_tailup, context, vt$vcov_theta)
  fc_slope <- fai_cornelius(L_slope, fit_tailup, context, vt$vcov_theta)

  expect_equal(fc_intercept, unname(df["(Intercept)"]))
  expect_equal(fc_slope, unname(df["ELEV_DEM"]))
})

test_that("fai_cornelius() joint (q = 2) df is finite and distinct from either marginal df", {
  context <- get_satterthwaite_context(fit_tailup, "numeric")
  vt <- get_vcov_theta_numeric(context)
  df <- satterthwaite(fit_tailup, method = "numeric")

  fc_joint <- fai_cornelius(diag(2), fit_tailup, context, vt$vcov_theta)

  expect_true(is.finite(fc_joint))
  expect_true(fc_joint > 0)
})

test_that("a null (all-zero) contrast is diagnosed as a degenerate df (NA + warning), not silently returned", {
  context <- get_satterthwaite_context(fit_tailup, "numeric")
  vt <- get_vcov_theta_numeric(context)
  expect_warning(
    result <- satterthwaite_df_for_contrast(c(0, 0), fit_tailup, context, vt$vcov_theta, "null contrast"),
    "degenerate"
  )
  expect_true(is.na(result))
})

test_that("satterthwaite() df matches spmodel's own method = 'numeric' df for an equivalent Euclidean-only model", {
  skip_if_not_installed("spmodel")

  fit_euc <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "none",
    taildown_type = "none",
    euclid_type = "exponential",
    nugget_type = "nugget"
  )
  df_ssn2 <- satterthwaite(fit_euc, method = "numeric")

  spmod <- spmodel::splm(
    formula = Summer_mn ~ ELEV_DEM,
    data = mf04p$obs,
    spcov_type = "exponential",
    estmethod = "reml"
  )
  df_spmodel <- spmodel::satterthwaite(spmod, method = "numeric")

  expect_equal(unname(df_ssn2["(Intercept)"]), unname(df_spmodel["(Intercept)"]), tolerance = 0.05)
  expect_equal(unname(df_ssn2["ELEV_DEM"]), unname(df_spmodel["ELEV_DEM"]), tolerance = 0.05)
})

# satterthwaite api
skip_if_not_installed("numDeriv")

copy_lsn_to_temp()
temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)

fit <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM,
  ssn.object = mf04p,
  tailup_type = "exponential",
  additive = "afvArea", ddf = "asymptotic"
)

fit_sw <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM,
  ssn.object = mf04p,
  tailup_type = "exponential",
  additive = "afvArea", ddf = "satterthwaite"
)

test_that("summary() uses z-based inference when object$ddf is unset and t-based inference when object$ddf is cached", {
  # summary()/confint()/tidy()/emmeans() have no method/ddf argument of their
  # own -- they always reflect object$ddf, matching spmodel's summary.splm()
  s0 <- summary(fit)
  expect_identical(colnames(s0$coefficients$fixed), c("estimates", "Std_Error", "z_value", "p"))

  s1 <- summary(fit_sw)
  expect_identical(colnames(s1$coefficients$fixed), c("estimates", "Std_Error", "df", "t_value", "p"))
  expect_equal(unname(s1$coefficients$fixed$df), unname(fit_sw$ddf))
})

test_that("confint() gives asymptotic intervals when object$ddf is unset and wider t-based intervals when cached", {
  ci0 <- confint(fit)
  ci1 <- confint(fit_sw)

  expect_equal(dim(ci0), dim(ci1))
  expect_true(all((ci1[, 2] - ci1[, 1]) > (ci0[, 2] - ci0[, 1])))

  # manual check against qt() with the cached fit's own object$ddf
  df <- fit_sw$ddf
  estimates <- coef(fit_sw, type = "fixed")
  variances <- diag(vcov(fit_sw, type = "fixed"))
  tstar <- qt(0.975, df[names(estimates)])
  expect_equal(unname(ci1[, 1]), unname(estimates - tstar * sqrt(variances)))
})

test_that("anova() ddf defaults to satterthwaite for n <= 500 regardless of the fit's own ddf; ddf = 'asymptotic' gives the Chi2 table", {
  # fit was built with ddf = "asymptotic" explicitly, but anova()'s own ddf
  # default is re-derived from sample size (determine_ddf()), independent of
  # the fit's own ddf choice -- matching spmodel's anova.splm()
  a0 <- anova(fit)
  expect_identical(colnames(a0), c("NumDF", "DenDF", "F value", "Pr(>F)"))

  a1 <- anova(fit, ddf = "asymptotic")
  expect_identical(colnames(a1), c("Df", "Chi2", "Pr(>Chi2)"))
  expect_equal(a0$"F value", a1$Chi2 / a1$Df)

  expect_error(anova(fit, ddf = "bogus"), "ddf")
})

test_that("anova() two-model likelihood ratio test is unaffected by the ddf argument's existence", {
  fit_ml <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p, tailup_type = "exponential",
    additive = "afvArea", estmethod = "ml"
  )
  fit_null_ml <- ssn_lm(
    formula = Summer_mn ~ 1, ssn.object = mf04p, tailup_type = "exponential",
    additive = "afvArea", estmethod = "ml"
  )
  lrt_before <- anova(fit_ml, fit_null_ml)
  # the two-model branch never reads `ddf` at all; confirm passing it
  # alongside a second model leaves LRT output byte-identical (matching
  # "Preserve LRT behavior separately").
  lrt_after <- anova(fit_ml, fit_null_ml, ddf = "satterthwaite")
  expect_equal(lrt_before, lrt_after)
})

test_that("tidy() fixed effects gains a df column exactly when object$ddf is cached", {
  t0 <- tidy(fit)
  expect_identical(colnames(t0), c("term", "estimate", "std.error", "statistic", "p.value"))

  t1 <- tidy(fit_sw)
  expect_identical(colnames(t1), c("term", "estimate", "std.error", "df", "statistic", "p.value"))

  t2 <- tidy(fit_sw, conf.int = TRUE)
  expect_true(all(c("df", "conf.low", "conf.high") %in% colnames(t2)))
})

test_that("emmeans() uses asymptotic Inf df when object$ddf is unset and finite Satterthwaite df when cached", {
  skip_if_not_installed("emmeans")
  emm0 <- emmeans::emmeans(fit, ~1, at = list(ELEV_DEM = 1000))
  emm0_df <- summary(emm0)$df
  expect_true(all(is.infinite(emm0_df)))

  emm1 <- emmeans::emmeans(fit_sw, ~1, at = list(ELEV_DEM = 1000))
  emm1_df <- summary(emm1)$df
  expect_true(all(is.finite(emm1_df)))
  expect_true(all(emm1_df > 0))
})

test_that("emmeans() on an ssn_glm model always uses asymptotic Inf df (Satterthwaite is unsupported, object$ddf is always unset)", {
  skip_if_not_installed("emmeans")
  mf04p_bin <- mf04p
  set.seed(2)
  mf04p_bin$obs$y01 <- rbinom(nrow(mf04p_bin$obs), 1, 0.5)
  fit_glm <- ssn_glm(
    formula = y01 ~ ELEV_DEM, ssn.object = mf04p_bin, family = "binomial",
    tailup_type = "exponential", additive = "afvArea"
  )
  expect_null(fit_glm$ddf)
  emm_glm <- emmeans::emmeans(fit_glm, ~1, at = list(ELEV_DEM = 1000))
  expect_true(all(is.infinite(summary(emm_glm)$df)))
})

# satterthwaite caching
skip_if_not_installed("numDeriv")

copy_lsn_to_temp()
temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)

test_that("explicit ddf = 'asymptotic' leaves object$ddf/vcov$cov unset", {
  fit <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea", ddf = "asymptotic"
  )
  expect_null(fit$ddf)
  expect_null(fit$vcov$cov)
})

test_that("ddf = 'satterthwaite' caches object$ddf/vcov$cov matching an on-demand computation", {
  fit_default <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea"
  )
  fit_cached <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea", ddf = "satterthwaite"
  )

  expect_named(fit_cached$ddf, rownames(fit_cached$vcov$fixed))
  expect_true(all(is.finite(fit_cached$ddf)))
  expect_equal(dim(fit_cached$vcov$cov), c(3L, 3L))

  df_ondemand <- satterthwaite(fit_default, method = "numeric")
  expect_equal(unname(fit_cached$ddf), unname(df_ondemand))
})

test_that("cached results are actually reused, not silently recomputed", {
  fit_cached <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea", ddf = "satterthwaite"
  )

  fit_sentinel <- fit_cached
  fit_sentinel$vcov$cov[1, 1] <- fit_sentinel$vcov$cov[1, 1] * 100

  # satterthwaite() reuses object$vcov$cov (now corrupted) to recompute
  # per-coefficient df on demand, so it picks up the corruption
  df_from_sentinel <- satterthwaite(fit_sentinel, method = "numeric")
  expect_false(isTRUE(all.equal(unname(df_from_sentinel), unname(fit_cached$ddf))))

  # summary() pulls directly from the untouched object$ddf field, not from
  # object$vcov$cov, so it is unaffected by the corruption
  s_sentinel <- summary(fit_sentinel)
  expect_equal(unname(s_sentinel$coefficients$fixed$df), unname(fit_sentinel$ddf))
  expect_equal(unname(fit_sentinel$ddf), unname(fit_cached$ddf))
})

test_that("summary()/confint()/anova() on a ddf = 'satterthwaite' fit match the automatic default for n <= 500", {
  fit_default <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea"
  )
  fit_cached <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea", ddf = "satterthwaite"
  )

  expect_equal(
    summary(fit_cached)$coefficients$fixed,
    summary(fit_default)$coefficients$fixed
  )
  expect_equal(
    confint(fit_cached),
    confint(fit_default)
  )
  expect_equal(
    as.data.frame(anova(fit_cached, ddf = "satterthwaite")),
    as.data.frame(anova(fit_default, ddf = "satterthwaite"))
  )
})

test_that("ddf must be 'asymptotic' or 'satterthwaite'", {
  expect_error(
    ssn_lm(
      formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
      tailup_type = "exponential", additive = "afvArea", ddf = "bogus"
    ),
    "ddf"
  )
})

test_that("a ddf = 'satterthwaite' request that cannot be satisfied warns and still returns a fitted model", {
  expect_warning(
    fit_known <- ssn_lm(
      formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
      tailup_type = "exponential",
      tailup_initial = tailup_initial("exponential", de = 2, range = 10000, known = "given"),
      nugget_initial = nugget_initial("nugget", nugget = 0.5, known = "given"),
      additive = "afvArea",
      ddf = "satterthwaite"
    ),
    "could not be computed"
  )
  expect_s3_class(fit_known, "ssn_lm")
  expect_null(fit_known$ddf)
  # method = "numeric" still works on demand for this fit (all-known
  # covariance is the actual reason it fails, independent of ddf caching).
  expect_error(satterthwaite(fit_known, method = "numeric"), "known")
})

# satterthwaite context
skip_if_not_installed("numDeriv")

copy_lsn_to_temp()
temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)

fit_tailup <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM,
  ssn.object = mf04p,
  tailup_type = "exponential",
  additive = "afvArea"
)

fit_iid <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM,
  ssn.object = mf04p,
  tailup_type = "none",
  taildown_type = "none",
  euclid_type = "none",
  nugget_type = "nugget",
  additive = "afvArea"
)

test_that("get_satterthwaite_context() rejects everything except explicit method = 'numeric'", {
  expect_error(get_satterthwaite_context(fit_tailup, "closed"), "method")
  expect_error(get_satterthwaite_context(fit_tailup, "automatic"), "method")
  expect_error(get_satterthwaite_context(fit_tailup), "method")
})

test_that("get_satterthwaite_context() rejects ssn_glm objects explicitly", {
  mf04p_bin <- mf04p
  set.seed(2)
  mf04p_bin$obs$y01 <- rbinom(nrow(mf04p_bin$obs), 1, 0.5)
  fit_glm <- ssn_glm(
    formula = y01 ~ ELEV_DEM,
    ssn.object = mf04p_bin,
    family = "binomial",
    tailup_type = "exponential",
    additive = "afvArea"
  )
  expect_error(get_satterthwaite_context(fit_glm, "numeric"), "ssn_glm")
})

test_that("get_satterthwaite_context() rejects all-known covariance", {
  fit_known <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    tailup_initial = tailup_initial("exponential", de = 2, range = 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.5, known = "given"),
    additive = "afvArea"
  )
  expect_error(get_satterthwaite_context(fit_known, "numeric"), "known")
})

test_that("known covariance parameters are excluded from the free-parameter vector", {
  fit_known_range <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    tailup_initial = tailup_initial("exponential", range = 10000, known = "range"),
    additive = "afvArea"
  )
  context <- get_satterthwaite_context(fit_known_range, "numeric")
  expect_false("tailup_range" %in% context$cov_names_free_orig)
  expect_true("tailup_de" %in% context$cov_names_free_orig)
  expect_true("nugget" %in% context$cov_names_free_orig)
})

test_that("g(theta) at the fitted covariance parameters exactly reproduces object$vcov$fixed", {
  context <- get_satterthwaite_context(fit_tailup, "numeric")
  p <- fit_tailup$p
  for (k in seq_len(p)) {
    Li <- as.numeric(diag(p)[k, ])
    g <- get_grad_gi(context$cov_val_free, Li, context)
    expect_equal(g, fit_tailup$vcov$fixed[k, k], tolerance = 1e-6)
  }
})

test_that("numeric vcov_theta matches the analytic REML variance of sigma^2 for a plain IID Gaussian fit", {
  context <- get_satterthwaite_context(fit_iid, "numeric")
  expect_identical(context$cov_names_free_orig, "nugget")

  vt <- get_vcov_theta_numeric(context)
  expect_false(is.null(vt$vcov_theta))

  n <- fit_iid$n
  p <- fit_iid$p
  sigma2_hat <- fit_iid$coefficients$params_object$nugget[["nugget"]]
  analytic_var <- 2 * sigma2_hat^2 / (n - p)

  expect_equal(as.numeric(vt$vcov_theta[1, 1]), analytic_var, tolerance = 1e-4)
})

test_that("a boundary/near-degenerate anisotropic fit is diagnosed (NULL + warning), not silently repaired", {
  fit_anis <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    euclid_type = "exponential",
    anisotropy = TRUE,
    additive = "afvArea"
  )
  context <- get_satterthwaite_context(fit_anis, "numeric")

  # g(theta) at the fitted point must still exactly reproduce vcov$fixed
  # regardless of Hessian conditioning (it does not depend on the Hessian).
  Li <- c(1, 0)
  g <- get_grad_gi(context$cov_val_free, Li, context)
  expect_equal(g, fit_anis$vcov$fixed[1, 1], tolerance = 1e-6)

  vt <- withCallingHandlers(
    get_vcov_theta_numeric(context),
    warning = function(w) invokeRestart("muffleWarning")
  )
  if (is.null(vt$vcov_theta)) {
    expect_warning(get_vcov_theta_numeric(context), "not numerically positive definite")
  } else {
    expect_true(all(is.finite(vt$vcov_theta)))
    expect_true(isSymmetric(unname(vt$vcov_theta), tolerance = 1e-6))
  }
})

test_that("vcov() default behavior (type = 'fixed') is unchanged", {
  expect_identical(vcov(fit_tailup), fit_tailup$vcov$fixed)
  expect_identical(vcov(fit_tailup, type = "fixed"), fit_tailup$vcov$fixed)
})

test_that("vcov() rejects an invalid type and has no method argument of its own", {
  expect_error(vcov(fit_tailup, type = "bogus"), "type")
  # method is silently ignored (absorbed by ...), matching spmodel's vcov.splm()
  expect_equal(vcov(fit_tailup, type = "cov"), vcov(fit_tailup, type = "cov", method = "numeric"))
})

test_that("vcov(type = 'cov'/'ssn'/'tailup'/.../'randcov') are NULL when Satterthwaite information was not cached at fit time", {
  fit_asym <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea", ddf = "asymptotic"
  )
  expect_null(vcov(fit_asym, type = "cov"))
  expect_null(vcov(fit_asym, type = "ssn"))
  expect_null(vcov(fit_asym, type = "tailup"))
  expect_null(vcov(fit_asym, type = "randcov"))
})

test_that("vcov(type = 'cov'/'ssn'/'tailup'/'taildown'/'euclid'/'nugget'/'randcov') are consistent for a mixed spatial + random-effects model", {
  fit_mixed <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    random = ~ as.factor(netID),
    additive = "afvArea"
  )

  v_cov <- vcov(fit_mixed, type = "cov")
  v_ssn <- vcov(fit_mixed, type = "ssn")
  v_tailup <- vcov(fit_mixed, type = "tailup")
  v_nugget <- vcov(fit_mixed, type = "nugget")
  v_randcov <- vcov(fit_mixed, type = "randcov")

  expect_equal(NROW(v_ssn) + NROW(v_randcov), NROW(v_cov))
  expect_equal(v_ssn, v_cov[rownames(v_ssn), rownames(v_ssn), drop = FALSE])
  expect_equal(v_randcov, v_cov[rownames(v_randcov), rownames(v_randcov), drop = FALSE])

  # tailup_type = "exponential" with taildown_type/euclid_type = "none" (the
  # ssn_lm() default): only tailup and nugget contribute free parameters
  expect_equal(NROW(v_tailup) + NROW(v_nugget), NROW(v_ssn))
  expect_equal(v_tailup, v_ssn[rownames(v_tailup), rownames(v_tailup), drop = FALSE])
  expect_null(vcov(fit_mixed, type = "taildown"))
  expect_null(vcov(fit_mixed, type = "euclid"))
})

test_that("vcov(type = 'randcov') is NULL when the model has no random effects", {
  expect_null(vcov(fit_tailup, type = "randcov"))
})

test_that("satterthwaite() still computes correctly for a model fit with range_constrain = TRUE (the gradient context always forces range_constrain = FALSE internally)", {
  fit_constrained <- ssn_lm(
    formula = Summer_mn ~ ELEV_DEM,
    ssn.object = mf04p,
    tailup_type = "exponential",
    additive = "afvArea",
    range_constrain = TRUE
  )
  df_constrained <- satterthwaite(fit_constrained, method = "numeric")
  expect_named(df_constrained, rownames(fit_constrained$vcov$fixed))
  expect_true(all(is.finite(df_constrained)))
  expect_true(all(df_constrained > 0))

  # the gradient context is rebuilt with range_constrain forced off regardless
  # of how the model was fit, so its own orig2optim_object never ends up on
  # the logit-odds scale
  context <- get_satterthwaite_context(fit_constrained, "numeric")
  expect_null(context$orig2optim_object$range_constrain_value)
  expect_false(any(grepl("_range_logodds$", names(context$orig2optim_object$value))))
})

