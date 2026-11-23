skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("fitting rejects formula variables absent from observation data", {
  outside_x <- mf04p$obs$ELEV_DEM
  outside_y <- mf04p$obs$Summer_mn
  outside_group <- factor(rep(c("a", "b"), length.out = length(outside_x)))
  cases <- list(
    list(formula = Summer_mn ~ outside_x),
    list(formula = outside_y ~ ELEV_DEM),
    list(formula = Summer_mn ~ ELEV_DEM + offset(outside_x)),
    list(formula = Summer_mn ~ ELEV_DEM, random = ~outside_group),
    list(formula = Summer_mn ~ ELEV_DEM, random = ~(outside_x | netID)),
    list(formula = Summer_mn ~ ELEV_DEM, partition_factor = ~outside_group),
    list(formula = Summer_mn ~ predictor_not_defined_anywhere)
  )
  for (local in list(FALSE, list(index = rep(1:3, length.out = length(outside_x))))) {
    for (case in cases) {
      args <- c(case, list(ssn.object = mf04p, local = local))
      expect_error(do.call(ssn_lm, c(args, list(ddf = "asymptotic"))),
                   "used in formula, random, or partition_factor not found in data", fixed = TRUE)
      expect_error(do.call(ssn_glm, c(args, list(family = "Gamma"))),
                   "used in formula, random, or partition_factor not found in data", fixed = TRUE)
    }
  }
  expect_error(ssn_glm(Summer_mn ~ outside_x, mf04p, family = "Gaussian"),
               "not found in data", fixed = TRUE)
})

test_that("formula validation preserves transformations and fixed-effect dot syntax", {
  data <- data.frame(y = 1:4, x = 2:5, g = factor(c("a", "b", "a", "b")))
  transform_x <- function(x) log(x)
  expect_null(check_formula_vars_in_data(y ~ transform_x(x) + offset(x), data,
                                        random = ~(x | g), partition_factor = ~g))
  expect_null(check_formula_vars_in_data(y ~ stats::poly(x, degree = 2), data))
  expect_null(check_formula_vars_in_data(y ~ ., data))
  for (fun in list(ssn_lm, function(...) ssn_glm(..., family = "Gamma"))) {
    expect_error(fun(Summer_mn ~ ELEV_DEM, mf04p, random = ~.),
                 "The `.` shorthand is not supported in random.", fixed = TRUE)
    expect_error(fun(Summer_mn ~ ELEV_DEM, mf04p, partition_factor = ~.),
                 "The `.` shorthand is not supported in partition_factor.", fixed = TRUE)
  }
  external_constant <- 2
  expect_error(check_formula_vars_in_data(y ~ I(x / external_constant), data),
               '"external_constant"', fixed = TRUE)
})

test_that("valid transformed formulas still fit exact and grouped models", {
  ssn <- mf04p
  ssn$obs$shift <- seq_len(nrow(ssn$obs)) / 100
  ssn$obs$Summer_mn[3] <- NA_real_
  transform_x <- function(x) x / 1000
  form <- Summer_mn ~ transform_x(ELEV_DEM) + offset(shift)
  for (local in list(FALSE, list(index = rep(1:3, length.out = nrow(ssn$obs) - 1)))) {
    lm_fit <- ssn_lm(form, ssn, local = local, ddf = "asymptotic",
                     nugget_initial = nugget_initial("nugget", nugget = 1, known = "given"))
    glm_fit <- ssn_glm(form, ssn, local = local,
                       dispersion_initial = dispersion_initial("Gamma", dispersion = 5, known = "given"),
                       nugget_initial = nugget_initial("nugget", nugget = 1, known = "given"))
    expect_equal(lm_fit$n, nrow(ssn$obs) - 1)
    expect_equal(glm_fit$n, nrow(ssn$obs) - 1)
    expect_identical(environment(formula(lm_fit)), environment(form))
    expect_identical(environment(formula(glm_fit)), environment(form))
    expect_true(all(is.finite(coef(lm_fit))))
    expect_true(all(is.finite(coef(glm_fit))))
  }
})

test_that("single-location polynomial blocks agree with point predictions", {
  ssn <- mf04p
  ssn$obs$px <- as.numeric(scale(ssn$obs$ELEV_DEM))
  ssn$obs$py <- seq_len(nrow(ssn$obs)) / nrow(ssn$obs)
  ssn$obs$shift <- seq_len(nrow(ssn$obs)) / 100
  ssn$preds$pred1km <- ssn$preds$pred1km[1:2, , drop = FALSE]
  ssn$preds$pred1km$px <- c(0.25, 0.35)
  ssn$preds$pred1km$py <- c(0.4, 0.6)
  ssn$preds$pred1km$shift <- c(2, 3)
  fit <- ssn_lm(Summer_mn ~ poly(px, py, degree = 2) + offset(shift), ssn,
                ddf = "asymptotic",
                euclid_initial = euclid_initial("exponential", de = 2, range = 10000, known = "given"),
                nugget_initial = nugget_initial("nugget", nugget = 1, known = "given"))
  grid <- fit$ssn.object$preds$pred1km
  frame <- model.frame(delete.response(terms(fit)), grid, xlev = fit$xlevels)
  design <- model.matrix(delete.response(terms(fit)), frame, contrasts = fit$contrasts)
  mean_design <- matrix(colMeans(design), nrow = 1)
  multi <- predict(fit, "pred1km", block = TRUE, interval = "confidence", se.fit = TRUE)
  expect_equal(unname(multi$fit[1, "fit"]),
               as.numeric(mean_design %*% coef(fit)) + mean(grid$shift))
  expect_equal(unname(multi$se.fit), sqrt(as.numeric(mean_design %*% vcov(fit) %*% t(mean_design))))

  fit$ssn.object$preds$pred1km <- grid[1, , drop = FALSE]
  for (interval in c("none", "prediction", "confidence")) {
    point <- predict(fit, "pred1km", interval = interval, se.fit = TRUE)
    block <- predict(fit, "pred1km", block = TRUE, interval = interval, se.fit = TRUE)
    expect_equal(unname(block$fit), unname(point$fit), tolerance = 1e-8)
    expect_equal(unname(block$se.fit), unname(point$se.fit), tolerance = 1e-8)
  }
  point_terms <- predict(fit, "pred1km", type = "terms", se.fit = TRUE)
  block_terms <- predict(fit, "pred1km", block = TRUE, type = "terms", se.fit = TRUE)
  expect_null(rownames(block_terms$fit))
  expect_null(rownames(block_terms$se.fit))
  rownames(point_terms$fit) <- rownames(point_terms$se.fit) <- NULL
  expect_equal(block_terms, point_terms, tolerance = 1e-8)
  expect_equal(as.numeric(SSN2:::get_newdata_model_matrix(fit, grid[1, , drop = FALSE])$newdata_model),
               as.numeric(design[1, ]))
})
