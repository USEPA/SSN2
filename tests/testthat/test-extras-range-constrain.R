skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# get_range_constrain_setup() ------------------------------------------------

test_that("get_range_constrain_setup() returns range_constrain_value = NULL and every component FALSE when range_constrain = FALSE", {
  initial_object_val <- get_initial_object(
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential", nugget_type = "nugget",
    tailup_initial = NULL, taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL
  )
  setup <- get_range_constrain_setup(mf04p$obs, initial_object_val, FALSE)
  expect_null(setup$range_constrain_value)
  expect_false(setup$tailup_range_constrain)
  expect_false(setup$taildown_range_constrain)
  expect_false(setup$euclid_range_constrain)
})

test_that("get_range_constrain_setup() constrains every active range parameter to 4x the bounding box diagonal when range_constrain = TRUE and no initial ranges are supplied", {
  initial_object_val <- get_initial_object(
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential", nugget_type = "nugget",
    tailup_initial = NULL, taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL
  )
  setup <- get_range_constrain_setup(mf04p$obs, initial_object_val, TRUE)

  bbox <- sf::st_bbox(mf04p$obs)
  bbox_dist <- sqrt((bbox[["xmax"]] - bbox[["xmin"]])^2 + (bbox[["ymax"]] - bbox[["ymin"]])^2)
  expect_equal(setup$range_constrain_value, 4 * bbox_dist)
  expect_true(setup$tailup_range_constrain)
  expect_true(setup$taildown_range_constrain)
  expect_true(setup$euclid_range_constrain)
})

test_that("get_range_constrain_setup() skips constraining a range parameter whose covariance type is 'none'", {
  initial_object_val <- get_initial_object(
    tailup_type = "none", taildown_type = "none", euclid_type = "exponential", nugget_type = "nugget",
    tailup_initial = NULL, taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL
  )
  setup <- get_range_constrain_setup(mf04p$obs, initial_object_val, TRUE)
  expect_false(setup$tailup_range_constrain)
  expect_false(setup$taildown_range_constrain)
  expect_true(setup$euclid_range_constrain)
})

test_that("get_range_constrain_setup() skips constraining a range parameter that is already known (fixed)", {
  initial_object_val <- get_initial_object(
    tailup_type = "exponential", taildown_type = "none", euclid_type = "none", nugget_type = "nugget",
    tailup_initial = tailup_initial("exponential", range = 5000, known = "range"),
    taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL
  )
  setup <- get_range_constrain_setup(mf04p$obs, initial_object_val, TRUE)
  expect_false(setup$tailup_range_constrain)
})

test_that("get_range_constrain_setup() skips constraining a range parameter whose own initial value already exceeds the bound", {
  initial_object_val <- get_initial_object(
    tailup_type = "exponential", taildown_type = "none", euclid_type = "none", nugget_type = "nugget",
    tailup_initial = tailup_initial("exponential", range = 1e9),
    taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL
  )
  setup <- get_range_constrain_setup(mf04p$obs, initial_object_val, TRUE)
  expect_false(setup$tailup_range_constrain)
})

test_that("get_range_constrain_setup() constrains freely when an initial range is explicitly NA (not yet known, no comparison basis)", {
  initial_object_val <- get_initial_object(
    tailup_type = "exponential", taildown_type = "none", euclid_type = "none", nugget_type = "nugget",
    tailup_initial = tailup_initial("exponential", range = NA),
    taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL
  )
  setup <- get_range_constrain_setup(mf04p$obs, initial_object_val, TRUE)
  expect_true(setup$tailup_range_constrain)
})

test_that("get_range_constrain_setup() errors cleanly on non-TRUE/FALSE input", {
  initial_object_val <- get_initial_object(
    tailup_type = "exponential", taildown_type = "none", euclid_type = "none", nugget_type = "nugget",
    tailup_initial = NULL, taildown_initial = NULL, euclid_initial = NULL, nugget_initial = NULL
  )
  expect_error(get_range_constrain_setup(mf04p$obs, initial_object_val, "yes"), "range_constrain must be TRUE or FALSE")
  expect_error(get_range_constrain_setup(mf04p$obs, initial_object_val, NA), "range_constrain must be TRUE or FALSE")
  expect_error(get_range_constrain_setup(mf04p$obs, initial_object_val, c(TRUE, FALSE)), "range_constrain must be TRUE or FALSE")
})

# orig2optim_range_component() / optim2orig_range_component() round-trip ----

test_that("orig2optim_range_component() / optim2orig_range_component() round-trip on the log scale when constrain = FALSE", {
  component_initial <- list(initial = c(de = 1, range = 8000), is_known = c(de = FALSE, range = FALSE))
  transformed <- orig2optim_range_component(component_initial, FALSE, NULL, "tailup")
  expect_equal(names(transformed$value), "tailup_range_log")
  expect_equal(unname(transformed$value), log(8000))

  recovered <- optim2orig_range_component(transformed$value, NULL, "tailup")
  expect_equal(recovered, 8000)
})

test_that("orig2optim_range_component() / optim2orig_range_component() round-trip on the logit-odds scale when constrain = TRUE", {
  component_initial <- list(initial = c(de = 1, range = 8000), is_known = c(de = FALSE, range = FALSE))
  transformed <- orig2optim_range_component(component_initial, TRUE, 100000, "euclid")
  expect_equal(names(transformed$value), "euclid_range_logodds")
  expect_equal(unname(transformed$value), logit(8000 / 100000))

  recovered <- optim2orig_range_component(transformed$value, 100000, "euclid")
  expect_equal(recovered, 8000, tolerance = 1e-8)
})

# orig2optim_ssn_components() / optim2orig_ssn_components() with a mix of constrained/unconstrained ranges ----

test_that("orig2optim_ssn_components() / optim2orig_ssn_components() round-trip correctly when only some ranges are constrained", {
  initial_object <- list(
    tailup_initial = list(initial = c(de = 2.5, range = 8000), is_known = c(de = FALSE, range = FALSE)),
    taildown_initial = list(initial = c(de = 1.2, range = 5000), is_known = c(de = FALSE, range = FALSE)),
    euclid_initial = list(
      initial = c(de = 0.8, range = 3000, rotate = 1.1, scale = 0.4),
      is_known = c(de = FALSE, range = FALSE, rotate = FALSE, scale = FALSE)
    ),
    nugget_initial = list(initial = c(nugget = 0.3), is_known = c(nugget = FALSE))
  )
  data_object <- list(
    range_constrain_value = 100000,
    tailup_range_constrain = TRUE,
    taildown_range_constrain = FALSE,
    euclid_range_constrain = TRUE
  )

  optim_val <- orig2optim_ssn_components(initial_object, data_object)
  expect_true("tailup_range_logodds" %in% names(optim_val$value))
  expect_true("taildown_range_log" %in% names(optim_val$value))
  expect_true("euclid_range_logodds" %in% names(optim_val$value))
  expect_false("tailup_range_log" %in% names(optim_val$value))
  expect_false("euclid_range_log" %in% names(optim_val$value))

  recovered <- optim2orig_ssn_components(optim_val$value, range_constrain_value = data_object$range_constrain_value)
  expect_equal(unname(recovered[["tailup_range"]]), 8000, tolerance = 1e-8)
  expect_equal(unname(recovered[["taildown_range"]]), 5000)
  expect_equal(unname(recovered[["euclid_range"]]), 3000, tolerance = 1e-8)
})

test_that("orig2optim_ssn_components() defaults to the unconstrained log scale when data_object is omitted", {
  initial_object <- list(
    tailup_initial = list(initial = c(de = 2.5, range = 8000), is_known = c(de = FALSE, range = FALSE)),
    taildown_initial = list(initial = c(de = 1.2, range = 5000), is_known = c(de = FALSE, range = FALSE)),
    euclid_initial = list(
      initial = c(de = 0.8, range = 3000, rotate = 1.1, scale = 0.4),
      is_known = c(de = FALSE, range = FALSE, rotate = FALSE, scale = FALSE)
    ),
    nugget_initial = list(initial = c(nugget = 0.3), is_known = c(nugget = FALSE))
  )
  optim_val <- orig2optim_ssn_components(initial_object)
  expect_true(all(c("tailup_range_log", "taildown_range_log", "euclid_range_log") %in% names(optim_val$value)))
})

# ssn_lm()/ssn_glm() integration -------------------------------------------

test_that("ssn_lm(range_constrain = TRUE) matches range_constrain = FALSE (up to optimizer tolerance) when the unconstrained optimum is well within the bound", {
  fit_unconstrained <- ssn_lm(
    Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    euclid_type = "exponential", additive = "afvArea",
    range_constrain = FALSE
  )
  fit_constrained <- ssn_lm(
    Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    euclid_type = "exponential", additive = "afvArea",
    range_constrain = TRUE
  )

  bbox <- sf::st_bbox(mf04p$obs)
  bbox_dist <- sqrt((bbox[["xmax"]] - bbox[["xmin"]])^2 + (bbox[["ymax"]] - bbox[["ymin"]])^2)
  unconstrained_range <- unname(coef(fit_unconstrained, type = "euclid")[["range"]])
  # sanity check that this fit is not vacuously within the bound
  expect_true(unconstrained_range < 4 * bbox_dist)

  expect_equal(
    unname(coef(fit_unconstrained, type = "euclid")[c("de", "range")]),
    unname(coef(fit_constrained, type = "euclid")[c("de", "range")]),
    tolerance = 0.05
  )
  expect_equal(as.numeric(logLik(fit_unconstrained)), as.numeric(logLik(fit_constrained)), tolerance = 1e-2)
})

test_that("ssn_lm(range_constrain = TRUE) does not error when a range parameter is known (fixed), and matches range_constrain = FALSE exactly", {
  fixed_fit <- function(range_constrain) {
    ssn_lm(
      Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
      tailup_type = "exponential", additive = "afvArea",
      tailup_initial = tailup_initial("exponential", range = 5000, known = "range"),
      range_constrain = range_constrain
    )
  }
  fit_f <- fixed_fit(FALSE)
  fit_t <- fixed_fit(TRUE)
  expect_equal(coef(fit_f, type = "tailup"), coef(fit_t, type = "tailup"))
  expect_equal(as.numeric(logLik(fit_f)), as.numeric(logLik(fit_t)), tolerance = 1e-6)
})

test_that("ssn_lm(range_constrain = TRUE) auto-disables constraining (without erroring) when the supplied initial range already exceeds the bound", {
  fit <- ssn_lm(
    Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea",
    tailup_initial = tailup_initial("exponential", range = 1e9),
    range_constrain = TRUE
  )
  expect_s3_class(fit, "ssn_lm")
})

test_that("ssn_lm() rejects a non-logical range_constrain", {
  expect_error(
    ssn_lm(
      Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
      tailup_type = "exponential", additive = "afvArea",
      range_constrain = "yes"
    ),
    "range_constrain must be TRUE or FALSE"
  )
})

test_that("ssn_glm(range_constrain = TRUE) threads through without erroring", {
  s <- mf04p
  s$obs$y_pois <- round(s$obs$Summer_mn)
  fit <- ssn_glm(
    y_pois ~ ELEV_DEM, ssn.object = s, family = "poisson",
    euclid_type = "exponential", additive = "afvArea",
    range_constrain = TRUE
  )
  expect_s3_class(fit, "ssn_glm")
})
