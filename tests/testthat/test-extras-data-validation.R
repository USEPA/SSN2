skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# validation dot formula
test_that("dot formula excludes SSN geometry/topology metadata but keeps real covariates", {
  obdata <- sf::st_drop_geometry(mf04p$obs)

  # no dot: data passed through unchanged
  r1 <- get_dot_formula_data(Summer_mn ~ ELEV_DEM, obdata, NULL)
  expect_identical(colnames(r1), colnames(obdata))

  # dot: topology columns dropped, ordinary covariates (including afvArea,
  # which is not part of the documented netgeom topology tuple) retained
  r2 <- get_dot_formula_data(Summer_mn ~ ., obdata, NULL)
  expect_false(any(c("netgeom", "netID", "rid", "upDist", "ratio", "pid", "locID") %in% colnames(r2)))
  expect_true(all(c("ELEV_DEM", "SLOPE", "afvArea") %in% colnames(r2)))

  # dot + a topology column explicitly named elsewhere in formula: retained
  r3 <- get_dot_formula_data(Summer_mn ~ . + rid, obdata, NULL)
  expect_true("rid" %in% colnames(r3))
  expect_false("netID" %in% colnames(r3)) # still unreferenced, still dropped

  # sf geometry list-column dropped when present and unreferenced
  r4 <- get_dot_formula_data(Summer_mn ~ ., mf04p$obs, "geometry")
  expect_false("geometry" %in% colnames(r4))
})

test_that("ssn_lm() fits a dot formula end to end, excluding topology columns", {
  mf04p_small <- mf04p
  keep <- c(
    "Summer_mn", "ELEV_DEM", "SLOPE", "AREAWTMAP",
    "rid", "pid", "ratio", "upDist", "afvArea", "locID", "netID", "netgeom",
    attributes(mf04p$obs)$sf_column
  )
  mf04p_small$obs <- mf04p$obs[, keep]

  fit <- ssn_lm(Summer_mn ~ ., mf04p_small, tailup_type = "exponential", additive = "afvArea")
  cn <- names(coef(fit))
  expect_true(all(c("ELEV_DEM", "SLOPE", "AREAWTMAP") %in% cn))
  expect_false(any(cn %in% c("netID", "rid", "upDist", "ratio", "pid", "locID", "netgeom")))

  # the resolved (dot-free) formula is usable downstream without a data
  # argument (e.g., summary()/formula() previously errored: "'.' in formula
  # and no 'data' argument", since data_object$formula retained the literal ".")
  expect_false("." %in% all.vars(formula(fit)))
  expect_s3_class(summary(fit), "summary.ssn_lm")
})

test_that("an explicitly named topology column still works outside a dot formula (unaffected baseline)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM + rid, mf04p, tailup_type = "exponential", additive = "afvArea")
  expect_true("rid" %in% names(coef(fit)))
})

test_that("local$index length is validated with an informative error", {
  n_obs <- nrow(mf04p$obs)

  expect_error(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
      tailup_type = "exponential", additive = "afvArea",
      local = list(index = seq_len(5), method = "kmeans")
    ),
    "local\\$index must have the same length"
  )

  expect_error(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
      tailup_type = "exponential", additive = "afvArea",
      local = list(index = seq_len(n_obs + 20), method = "kmeans")
    ),
    "local\\$index must have the same length"
  )

  # valid edge case: correct-length index still fits
  ssn_create_bigdist(mf04p, overwrite = TRUE)
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    local = list(index = rep(1:3, length.out = n_obs), method = "kmeans")
  )
  expect_equal(fit$n, n_obs)
})

test_that("random effect grouping variables with NA are rejected instead of silently dropped", {
  mf04p_na <- mf04p
  mf04p_na$obs$randgrp <- factor(rep(c("a", "b", "c"), length.out = nrow(mf04p_na$obs)))
  mf04p_na$obs$randgrp[3] <- NA

  expect_error(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p_na,
      tailup_type = "exponential", additive = "afvArea", random = ~randgrp
    ),
    "Missing values found in random effect"
  )

  # valid edge case: no NA, fits as before
  mf04p_na$obs$randgrp[3] <- "a"
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p_na,
    tailup_type = "exponential", additive = "afvArea", random = ~randgrp
  )
  expect_equal(fit$n, nrow(mf04p_na$obs))
})

test_that("partition_factor variables with NA are rejected instead of silently dropped", {
  mf04p_na <- mf04p
  mf04p_na$obs$partcol <- factor(rep(c("x", "y"), length.out = nrow(mf04p_na$obs)))
  mf04p_na$obs$partcol[5] <- NA

  expect_error(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p_na,
      tailup_type = "exponential", additive = "afvArea", partition_factor = ~partcol
    ),
    "Missing values found in partition_factor"
  )

  # valid edge case: no NA, fits as before
  mf04p_na$obs$partcol[5] <- "x"
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p_na,
    tailup_type = "exponential", additive = "afvArea", partition_factor = ~partcol
  )
  expect_equal(fit$n, nrow(mf04p_na$obs))
})

# nobs
test_that("nobs() dispatches without error and matches fit$n / logLik()'s nobs attribute, ssn_lm and ssn_glm", {
  fit_lm <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")
  expect_true(is.numeric(nobs(fit_lm)))
  expect_length(nobs(fit_lm), 1)
  expect_equal(nobs(fit_lm), fit_lm$n)
  expect_equal(nobs(fit_lm), attr(logLik(fit_lm), "nobs"))

  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fit_glm <- ssn_glm(count_response ~ ELEV_DEM, s, family = "poisson", tailup_type = "exponential", additive = "afvArea")
  expect_true(is.numeric(nobs(fit_glm)))
  expect_length(nobs(fit_glm), 1)
  expect_equal(nobs(fit_glm), fit_glm$n)
  expect_equal(nobs(fit_glm), attr(logLik(fit_glm), "nobs"))
})

test_that("nobs() agrees with an independently-recomputed AICc() (cross-check against a separate existing code path, not just object$n read back)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea", estmethod = "ml")
  n_est_param <- fit$npar + fit$p
  aicc_reference <- -2 * as.numeric(logLik(fit)) + 2 * nobs(fit) * n_est_param / (nobs(fit) - n_est_param - 1)
  expect_equal(unclass(AICc(fit)), aicc_reference, tolerance = 1e-10)
})

test_that("nobs() reflects response missingness (excludes NA-response rows), ssn_lm and ssn_glm", {
  s <- mf04p
  s$obs$Summer_mn[c(2, 5)] <- NA
  fit_lm <- ssn_lm(Summer_mn ~ ELEV_DEM, s, tailup_type = "exponential", additive = "afvArea")
  expect_equal(nobs(fit_lm), sum(!is.na(s$obs$Summer_mn)))
  expect_equal(nobs(fit_lm), nrow(s$obs) - 2)

  s2 <- mf04p
  s2$obs$count_response <- rpois(nrow(s2$obs), lambda = 5)
  s2$obs$count_response[c(3, 7, 11)] <- NA
  fit_glm <- ssn_glm(count_response ~ ELEV_DEM, s2, family = "poisson", tailup_type = "exponential", additive = "afvArea")
  expect_equal(nobs(fit_glm), sum(!is.na(s2$obs$count_response)))
  expect_equal(nobs(fit_glm), nrow(s2$obs) - 3)
})

test_that("nobs() for a two-column binomial response counts data rows/trial-groups, not total trials", {
  s <- mf04p
  s$obs$success <- rep(1:3, length.out = nrow(s$obs))
  s$obs$failure <- seq_len(nrow(s$obs)) + 4

  fit <- ssn_glm(cbind(success, failure) ~ ELEV_DEM, s,
    family = "binomial", tailup_type = "exponential", additive = "afvArea"
  )
  expect_equal(nobs(fit), nrow(s$obs))
  expect_false(isTRUE(all.equal(nobs(fit), sum(s$obs$success + s$obs$failure))))
})

# argument order
test_that("new conditional and kcv methods follow upstream positional conventions", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
                tailup_type = "exponential", additive = "afvArea")
  count <- ssn_glm(C16 ~ ELEV_DEM, mf04p, family = "poisson",
                   tailup_type = "exponential", additive = "afvArea")
  fit$ssn.object$preds$CapeHorn <- fit$ssn.object$preds$CapeHorn[1:3, ]
  count$ssn.object$preds$CapeHorn <- count$ssn.object$preds$CapeHorn[1:3, ]
  local <- list(approximation = "vecchia", method = "covariance", size = 4)
  set.seed(2)
  positional <- conditional(fit, "CapeHorn", "newdata", 2, local, FALSE)
  set.seed(2)
  named <- conditional(fit, "CapeHorn", samples = 2, local = local, simulate_covparams = FALSE)
  expect_equal(positional, named)
  set.seed(2)
  positional <- conditional(count, "CapeHorn", "newdata", "new", 2, local, NULL)
  set.seed(2)
  named <- conditional(count, "CapeHorn", type = "new", samples = 2, local = local, newdata_size = NULL)
  expect_equal(positional, named)

  folds <- rep(1:3, length.out = nrow(mf04p$obs))
  expect_equal(
    kcv(fit, 3, TRUE, TRUE, FALSE, "prediction", 0.9, folds),
    kcv(fit, k = 3, cv_predict = TRUE, se.fit = TRUE, local = FALSE,
        interval = "prediction", level = 0.9, folds_index = folds)
  )
  expect_equal(
    kcv(count, 3, TRUE, "response", TRUE, TRUE, FALSE, folds),
    kcv(count, k = 3, cv_predict = TRUE, type = "response", se.fit = TRUE,
        delta = TRUE, local = FALSE, folds_index = folds)
  )
})

test_that("decorrelation places learner controls before anisotropy and random effects", {
  testthat::local_mocked_bindings(
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) mean(y),
    predict_decorrelate_algorithm = function(fit, X, algorithm) rep(fit, NROW(X))
  )
  training <- list(training_index = 1:30, test_index = 31:45)
  positional <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, mf04p, "exponential", "none", "none", "nugget",
    NULL, NULL, NULL, NULL, "afvArea", "ranger", "RMSPE", training, TRUE, FALSE,
    NULL, NULL, NULL, "none", FALSE, NULL, FALSE
  )
  named <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea",
    algorithm = "ranger", statistic = "RMSPE", training = training, evaluate_test = TRUE,
    anisotropy = FALSE, ordering = "none", local = FALSE, dense_grid = FALSE
  )
  expect_equal(positional$grid, named$grid)
  expect_equal(positional$decorrelate_data$params, named$decorrelate_data$params)
})

test_that("forest summary printing accepts positional significance-star controls", {
  received <- NULL
  testthat::local_mocked_bindings(
    print.summary.ssn_lm = function(x, digits, signif.stars, ...) {
      received <<- list(digits = digits, signif.stars = signif.stars)
      invisible(x)
    }
  )
  x <- structure(list(ranger = NULL, ssn_lm = structure(list(), class = "summary.ssn_lm")),
                 class = "summary.ssn_lmRF")
  invisible(capture.output(print(x, 3, FALSE)))
  expect_identical(received, list(digits = 3, signif.stars = FALSE))
})

# tidy conf int order
test_that("tidy(conf.int = TRUE) preserves natural (nonalphabetic) coefficient order, ssn_lm", {
  s <- mf04p
  s$obs$Zpred <- s$obs$ELEV_DEM
  s$obs$Apred <- s$obs$SLOPE
  fit <- ssn_lm(Summer_mn ~ Zpred + Apred, s, tailup_type = "exponential", additive = "afvArea")

  natural_order <- names(coef(fit))
  expect_identical(tidy(fit)$term, natural_order)
  expect_identical(tidy(fit, conf.int = TRUE)$term, natural_order)
  # the bug would have sorted this alphabetically: (Intercept), Apred, Zpred
  expect_false(identical(tidy(fit, conf.int = TRUE)$term, sort(natural_order)))
})

test_that("tidy(conf.int = TRUE) preserves natural coefficient order, ssn_glm", {
  s <- mf04p
  s$obs$Zpred <- s$obs$ELEV_DEM
  s$obs$Apred <- s$obs$SLOPE
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)
  fit <- ssn_glm(count_response ~ Zpred + Apred, s,
    family = "poisson", tailup_type = "exponential", additive = "afvArea"
  )

  natural_order <- names(coef(fit))
  expect_identical(tidy(fit, conf.int = TRUE)$term, natural_order)
  expect_false(identical(tidy(fit, conf.int = TRUE)$term, sort(natural_order)))
})

test_that("tidy(conf.int = TRUE) estimate/CI values are unaffected by the ordering fix (only row order changed)", {
  s <- mf04p
  s$obs$Zpred <- s$obs$ELEV_DEM
  s$obs$Apred <- s$obs$SLOPE
  fit <- ssn_lm(Summer_mn ~ Zpred + Apred, s, tailup_type = "exponential", additive = "afvArea")

  td <- tidy(fit, conf.int = TRUE)
  ci_reference <- confint(fit, type = "fixed")
  for (term in td$term) {
    expect_equal(td$estimate[td$term == term], unname(coef(fit)[term]), tolerance = 1e-10)
    expect_equal(td$conf.low[td$term == term], unname(ci_reference[term, 1]), tolerance = 1e-10)
    expect_equal(td$conf.high[td$term == term], unname(ci_reference[term, 2]), tolerance = 1e-10)
  }
})

test_that("tidy(effects = 'ssn')/'randcov' effect and is_known columns and stream labels are unaffected", {
  s <- mf04p
  s$obs$Zgroup <- factor(rep(c("g1", "g2"), length.out = nrow(s$obs)))
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_type = "exponential", additive = "afvArea", random = ~Zgroup
  )
  td_ssn <- tidy(fit, effects = "ssn")
  expect_true(all(c("effect", "term", "estimate", "is_known") %in% names(td_ssn)))
  expect_true(all(c("tailup", "taildown", "euclid", "nugget") %in% td_ssn$effect))

  td_rc <- tidy(fit, effects = "randcov")
  expect_true(all(c("term", "estimate", "is_known") %in% names(td_rc)))
  expect_true(any(grepl("Zgroup", td_rc$term)))
})

