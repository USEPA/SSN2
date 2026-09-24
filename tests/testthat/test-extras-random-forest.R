skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# ssn_lmRF
skip_if_not_installed("ranger")

rf_initials <- function() {
  list(
    tailup = tailup_initial("exponential", de = 1, range = 10000, known = c("de", "range")),
    euclid = euclid_initial("exponential", de = 1, range = 10000, known = c("de", "range")),
    nugget = nugget_initial("nugget", nugget = 0.1, known = "given")
  )
}

fit_ssn_lmRF <- function(ssn, response = "Summer_mn", local, ...) {
  initials <- rf_initials()
  args <- list(
    formula = stats::as.formula(paste(response, "~ ELEV_DEM")),
    ssn.object = ssn,
    tailup_type = "exponential",
    tailup_initial = initials$tailup,
    nugget_initial = initials$nugget,
    additive = "afvArea",
    num.trees = 80,
    seed = 2,
    num.threads = 1
  )
  if (!missing(local)) args$local <- local
  do.call(ssn_lmRF, c(args, list(...)))
}

test_that("ssn_lmRF uses OOB residuals and direct residual kriging", {
  ssn <- mf04p
  obs_before <- ssn$obs
  fit <- fit_ssn_lmRF(ssn)

  forest_reference <- ranger::ranger(
    Summer_mn ~ ELEV_DEM, data = sf::st_drop_geometry(ssn$obs),
    num.trees = 80, seed = 2, num.threads = 1,
    write.forest = TRUE, oob.error = TRUE
  )
  expect_s3_class(fit, "ssn_lmRF")
  expect_equal(fit$ranger$predictions, forest_reference$predictions, tolerance = 0)
  expect_equal(fit$ssn_lm$ssn.object$obs[[fit$residual_name]],
    ssn$obs$Summer_mn - fit$ranger$predictions,
    tolerance = 1e-12
  )
  expect_identical(ssn$obs, obs_before)
  expect_false(fit$residual_name %in% names(ssn$obs))

  residual_model <- fit$ssn_lm
  residual <- as.numeric(model.response(model.frame(residual_model)))
  covariance <- covmatrix(residual_model)
  cross_covariance <- covmatrix(residual_model, "CapeHorn")
  one <- matrix(1, nrow = length(residual), ncol = 1)
  beta <- as.numeric(solve(crossprod(one, solve(covariance, one)), crossprod(one, solve(covariance, residual))))
  residual_prediction <- as.numeric(beta + cross_covariance %*% solve(covariance, residual - beta))
  forest_prediction <- stats::predict(
    fit$ranger,
    data = as.data.frame(sf::st_drop_geometry(residual_model$ssn.object$preds$CapeHorn)),
    type = "response"
  )$predictions
  expect_equal(
    unname(predict(fit, "CapeHorn")), forest_prediction + residual_prediction,
    tolerance = 1e-10
  )
  expect_equal(fit$ranger$num.trees, 80)
  expect_s3_class(summary(fit), "summary.ssn_lmRF")
  expect_true(is.function(getS3method("print", "ssn_lmRF", optional = TRUE)))
  expect_true(is.function(getS3method("print", "summary.ssn_lmRF", optional = TRUE)))
  expect_output(print(fit), "Random forest mean model")
  expect_output(print(summary(fit)), "SSN residual model")
})

test_that("ssn_lmRF preserves missing rows and supports all prediction sets", {
  ssn <- mf04p
  ssn$obs$rf_response <- ssn$obs$Summer_mn
  missing_rows <- c(3, 17)
  missing_pid <- ssn_get_netgeom(ssn$obs, "pid")$pid[missing_rows]
  ssn$obs$rf_response[missing_rows] <- NA_real_
  ssn$obs$.ssn_lmRF_residual <- -999
  original_residual <- ssn$obs$.ssn_lmRF_residual

  fit <- fit_ssn_lmRF(ssn, response = "rf_response")
  expect_identical(ssn$obs$.ssn_lmRF_residual, original_residual)
  expect_identical(fit$missing_prediction_name, ".missing")
  expect_identical(fit$residual_name, ".ssn_lmRF_residual.1")
  expect_equal(
    ssn_get_netgeom(fit$ssn_lm$ssn.object$preds$.missing, "pid")$pid,
    missing_pid
  )
  expect_equal(predict(fit), predict(fit, ".missing"), tolerance = 1e-12)
  expect_equal(names(predict(fit)), as.character(missing_rows))
  all_predictions <- predict(fit, "all")
  expect_named(all_predictions, c("pred1km", "CapeHorn", ".missing"))
  expect_equal(all_predictions$.missing, predict(fit, ".missing"), tolerance = 1e-12)
})

test_that("ssn_lmRF respects the observed PID order", {
  reference <- fit_ssn_lmRF(mf04p)
  ssn <- mf04p
  permutation <- c(seq.int(2, nrow(ssn$obs)), 1L)
  ssn$obs <- ssn$obs[permutation, ]
  fit <- fit_ssn_lmRF(ssn)

  expect_identical(
    ssn_get_netgeom(fit$ssn_lm$ssn.object$obs, "pid")$pid,
    ssn_get_netgeom(ssn$obs, "pid")$pid
  )
  expect_equal(
    covmatrix(fit$ssn_lm),
    covmatrix(reference$ssn_lm)[permutation, permutation],
    tolerance = 1e-12
  )
  expect_true(all(is.finite(predict(fit, "CapeHorn"))))
})

test_that("ssn_lmRF routes supported arguments and rejects unsupported contracts", {
  ssn <- mf04p
  ssn$obs$rf_offset <- seq_len(nrow(ssn$obs)) / 100
  fit <- fit_ssn_lmRF(ssn)
  routed <- ssn_lmRF(
    Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential",
    tailup_initial = tailup_initial("exponential", de = 1, range = 10000, known = c("de", "range")),
    additive = "afvArea", num.trees = 40, seed = 2, num.threads = 1,
    control = list(reltol = 1e-4, maxit = 500)
  )
  expect_equal(routed$ssn_lm$optim$control$reltol, 1e-4)
  unquoted <- ssn_lmRF(
    Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential",
    tailup_initial = tailup_initial("exponential", de = 1, range = 10000, known = c("de", "range")),
    additive = afvArea, num.trees = 40, seed = 2, num.threads = 1
  )
  expect_equal(predict(unquoted, "CapeHorn"), predict(routed, "CapeHorn"), tolerance = 1e-12)
  estimated <- ssn_lmRF(
    Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential",
    tailup_initial = tailup_initial("exponential", de = 1, range = 10000, known = c("de", "range")),
    nugget_initial = nugget_initial("nugget", nugget = 0.1),
    additive = "afvArea", num.trees = 40, seed = 2, num.threads = 1
  )
  expect_identical(estimated$ssn_lm$optim$method, "Brent")

  expect_error(
    fit_ssn_lmRF(ssn, write.forest = FALSE),
    "requires write.forest = TRUE"
  )
  expect_error(
    fit_ssn_lmRF(ssn, oob.error = FALSE),
    "requires oob.error = TRUE"
  )
  expect_error(
    fit_ssn_lmRF(ssn, unsupported_argument = 1),
    "Unsupported ssn_lmRF argument"
  )
  expect_error(
    ssn_lmRF(
      Summer_mn ~ ELEV_DEM + offset(rf_offset), ssn,
      num.trees = 20, seed = 2
    ),
    "Offsets are not supported"
  )
  expect_error(predict(fit, "CapeHorn", se.fit = TRUE), "se.fit is not supported")
  expect_error(predict(fit, "CapeHorn", interval = "confidence"), "interval must be")
  expect_error(predict(fit, "CapeHorn", type = "terms"), "type must be")
  expect_error(predict(fit, "CapeHorn", block = TRUE), "block prediction")
  expect_error(predict(fit, "CapeHorn", scale = 1), "Unsupported ssn_lmRF prediction")
  expect_error(predict(fit, "CapeHorn", unsupported_argument = 1), "Unsupported ssn_lmRF prediction")
})

test_that("ssn_lmRF retains random, partition, and one-group local contracts", {
  ssn <- mf04p
  ssn$obs$rf_group <- factor(rep(1:3, length.out = nrow(ssn$obs)))
  ssn$obs$rf_partition <- factor(rep(1:2, length.out = nrow(ssn$obs)))
  ssn$preds$CapeHorn$rf_group <- factor(
    rep(1:3, length.out = nrow(ssn$preds$CapeHorn)), levels = levels(ssn$obs$rf_group)
  )
  ssn$preds$CapeHorn$rf_partition <- factor(
    rep(1:2, length.out = nrow(ssn$preds$CapeHorn)), levels = levels(ssn$obs$rf_partition)
  )
  initials <- rf_initials()
  common <- list(
    tailup_type = "exponential", tailup_initial = initials$tailup,
    nugget_initial = initials$nugget, additive = "afvArea",
    random = ~rf_group,
    randcov_initial = randcov_initial(rf_group = 0.1, known = "given"),
    partition_factor = ~rf_partition,
    num.trees = 80, seed = 2, num.threads = 1
  )
  exact <- do.call(ssn_lmRF, c(list(Summer_mn ~ ELEV_DEM, ssn), common))
  ssn_create_bigdist(ssn, predpts = "CapeHorn", overwrite = TRUE, no_cores = 1, verbose = FALSE)
  local <- do.call(
    ssn_lmRF,
    c(list(Summer_mn ~ ELEV_DEM, ssn), common,
      list(local = list(index = rep(1, nrow(ssn$obs)), parallel = FALSE)))
  )
  expect_equal(
    predict(local, "CapeHorn", local = FALSE), predict(exact, "CapeHorn"),
    tolerance = 1e-10
  )
})

