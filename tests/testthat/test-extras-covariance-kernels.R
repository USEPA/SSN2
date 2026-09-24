skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# euclid circular
test_that("circular covariance has the expected two-dimensional support", {
  coordinates <- rbind(c(0, 0), c(0, 5), c(10, 0))
  h <- Matrix::Matrix(as.matrix(dist(coordinates)), sparse = TRUE)
  params <- euclid_params("circular", de = 2, range = 10)
  expected <- diag(2, 3)
  expected[1, 2] <- expected[2, 1] <- 2 * (2 / 3 - sqrt(3) / (2 * pi))

  expect_s3_class(params, "euclid_circular")
  expect_equal(unname(as.matrix(cov_matrix(params, list(euclid_mat = h), FALSE))), expected)
  expect_equal(unname(as.matrix(cov_vector(params, list(euclid_pred_mat = h), FALSE))), expected)
  expect_equal(get_effective_range(params), 10)
  reference <- spmodel::spcov_params("circular", de = 2, ie = 0.5, range = 10)
  reference_cov <- getFromNamespace("spcov_matrix", "spmodel")(reference, as.matrix(h))
  expect_equal(unname(reference_cov), expected + diag(0.5, 3))
})

test_that("constructors reject one-dimensional covariance names", {
  for (constructor in list(euclid_params, euclid_initial)) {
    expect_error(constructor("cosine", de = 1, range = 10), "Use 'circular'")
    expect_error(constructor("triangular", de = 1, range = 10), "not a valid Euclidean")
  }
  expect_error(ssn_lm(Summer_mn ~ 1, mf04p, euclid_type = "cosine"), "Use 'circular'")
  expect_error(ssn_decorrelate_grid(Summer_mn ~ 1, mf04p, euclid_type = "cosine"), "Use 'circular'")
})

test_that("circular models, simulation, and decorrelation use the same covariance", {
  network <- mf04p
  network$preds$CapeHorn <- network$preds$CapeHorn[1:4, ]
  initial <- euclid_initial("circular", de = 0.4, range = 12000, rotate = 0, scale = 1, known = "given")
  nugget <- nugget_initial("nugget", nugget = 0.2, known = "given")
  fit <- ssn_lm(Summer_mn ~ 1, network, euclid_type = "circular",
    euclid_initial = initial, nugget_initial = nugget)
  sigma <- as.matrix(covmatrix(fit))
  weights <- solve(sigma, rep(1, nobs(fit)))
  expect_equal(unname(coef(fit)), sum(weights * network$obs$Summer_mn) / sum(weights))
  expect_equal(unname(vcov(fit)[1, 1]), 1 / sum(weights))
  expect_s3_class(fit$coefficients$params_object$euclid, "euclid_circular")
  expect_equal(covmatrix(fit, "CapeHorn", cov_type = "obs.pred"),
    t(covmatrix(fit, "CapeHorn", cov_type = "pred.obs")))
  expect_true(all(is.finite(predict(fit, "CapeHorn", se.fit = TRUE)$se.fit)))

  fit_glm <- ssn_glm(C16 ~ 1, network, family = "poisson", euclid_type = "circular",
    euclid_initial = initial, nugget_initial = nugget)
  expect_equal(as.matrix(covmatrix(fit_glm)), sigma)
  expect_true(all(is.finite(predict(fit_glm, "CapeHorn", se.fit = TRUE)$se.fit)))

  set.seed(2)
  draws <- ssn_simulate(family = "Gaussian", ssn.object = network, tailup_params = tailup_params("none"),
    taildown_params = taildown_params("none"), euclid_params = euclid_params("circular", 0.4, 12000),
    nugget_params = nugget_params("nugget", 0.2), samples = 2)
  set.seed(2)
  expected_draws <- t(chol(sigma)) %*% matrix(rnorm(nobs(fit) * 2), nobs(fit), 2)
  expect_equal(unname(as.matrix(draws)), unname(expected_draws))

  for (local in list(FALSE, list(method = "covariance", size = nobs(fit)))) {
    transformed <- ssn_decorrelate_data(Summer_mn ~ 1, network,
      euclid_params = euclid_params("circular", de = 0.4, range = 12000, rotate = 0, scale = 1),
      nugget_params = nugget_params("nugget", nugget = 0.2), ordering = "none", local = local)
    expect_equal(as.numeric(crossprod(transformed$tX)), 0.6 * sum(weights), tolerance = 1e-8)
    expect_equal(as.numeric(crossprod(transformed$tX, transformed$ty)),
      0.6 * sum(weights * network$obs$Summer_mn), tolerance = 1e-8)
    newdata <- ssn_decorrelate_newdata(transformed, "CapeHorn")
    expect_true(all(is.finite(newdata$tX_newdata)))
  }
  grid <- ssn_decorrelate_grid(Summer_mn ~ 1, network, euclid_type = "circular", dense_grid = FALSE)
  expect_setequal(grid$euclid_type, c("circular", "none"))
})

# euclid extra
test_that("Euclidean extra-family covariance formulas and identities hold", {
  h <- Matrix::Matrix(c(0, 2, 2, 0), nrow = 2, sparse = TRUE)
  observed <- list(euclid_mat = h)
  predicted <- list(euclid_pred_mat = h)

  matern <- euclid_params("matern", de = 3, range = 5, extra = 0.5)
  cauchy <- euclid_params("cauchy", de = 3, range = 5, extra = 1)
  pexponential <- euclid_params("pexponential", de = 3, range = 5, extra = 1)

  expect_equal(
    as.matrix(SSN2:::cov_matrix(matern, observed, FALSE)),
    as.matrix(SSN2:::cov_matrix(euclid_params("exponential", 3, 5), observed, FALSE))
  )
  expect_equal(
    as.matrix(SSN2:::cov_matrix(cauchy, observed, FALSE)),
    as.matrix(SSN2:::cov_matrix(euclid_params("rquad", 3, 5), observed, FALSE))
  )
  expect_equal(
    as.matrix(SSN2:::cov_matrix(pexponential, observed, FALSE)),
    as.matrix(SSN2:::cov_matrix(euclid_params("exponential", 3, 5), observed, FALSE))
  )

  pexponential2 <- euclid_params("pexponential", de = 3, range = 25, extra = 2)
  expect_equal(
    as.matrix(SSN2:::cov_vector(pexponential2, predicted, FALSE)),
    as.matrix(SSN2:::cov_vector(euclid_params("gaussian", 3, 5), predicted, FALSE))
  )
  expect_equal(diag(as.matrix(SSN2:::cov_matrix(matern, observed, FALSE))), c(3, 3))
})

test_that("Euclidean extra constructors preserve existing types and validate domains", {
  matern <- euclid_params("matern", de = 1, range = 2, extra = 1)
  expect_s3_class(matern, "euclid_matern")
  expect_named(matern, c("de", "range", "extra", "rotate", "scale"))
  expect_error(euclid_params("exponential", de = 1, range = 2, extra = 1), "extra")
  expect_error(euclid_params("matern", de = 1, range = 2, extra = 0.1), "extra")
  expect_error(euclid_params("cauchy", de = 1, range = 2, extra = 0), "extra")
  expect_error(euclid_params("pexponential", de = 1, range = 2, extra = 2.1), "extra")

  initial <- euclid_initial("matern", extra = NA)
  initial_na <- SSN2:::euclid_initial_NA(initial, list(anisotropy = FALSE))
  expect_named(initial_na$initial, c("de", "range", "extra", "rotate", "scale"))
  expect_false(initial_na$is_known[["extra"]])
})

test_that("Euclidean extra transforms and spmodel range starts are family aware", {
  for (type in c("matern", "cauchy", "pexponential")) {
    initial <- list(
      tailup_initial = list(initial = c(de = 1, range = 1), is_known = c(de = TRUE, range = TRUE)),
      taildown_initial = list(initial = c(de = 1, range = 1), is_known = c(de = TRUE, range = TRUE)),
      euclid_initial = SSN2:::euclid_initial_NA(
        euclid_initial(type, de = 1, range = 10, extra = 1, rotate = 0, scale = 1, known = "given"),
        list(anisotropy = FALSE)
      ),
      nugget_initial = list(initial = c(nugget = 1), is_known = c(nugget = TRUE))
    )
    transformed <- SSN2:::orig2optim_ssn_components(initial)
    restored <- SSN2:::optim2orig_ssn_components(transformed$value, type)
    expect_equal(restored[["euclid_extra"]], 1)
  }

  ranges <- rep(c(25, 75), 5)
  starts <- SSN2:::get_euclid_start_values(
    euclid_initial("pexponential"), ranges, euclid_max = 100
  )
  expect_equal(unique(starts$extra), c(0.4, 1.6))
  expect_equal(starts$range[[1]], 25 / 3)

  fixed_starts <- SSN2:::get_euclid_start_values(
    euclid_initial("pexponential", extra = 1), ranges, euclid_max = 100
  )
  expect_equal(unique(fixed_starts$extra), 1)
  expect_equal(fixed_starts$range[[1]], 25 / 3)
})

test_that("added Euclidean families work through fitted covariance requests", {
  copy_lsn_to_temp()
  ssn.object <- ssn_import(
    file.path(tempdir(), "MiddleFork04.ssn"), predpts = "CapeHorn", overwrite = TRUE
  )

  for (type in c("matern", "cauchy", "pexponential")) {
    range <- if (type == "pexponential") 1e8 else 10000
    euclid_initial_val <- euclid_initial(
      type, de = 1, range = range, extra = 1, known = "given"
    )
    fit <- ssn_lm(
      Summer_mn ~ ELEV_DEM, ssn.object,
      tailup_type = "none", taildown_type = "none", euclid_type = type,
      nugget_type = "nugget", euclid_initial = euclid_initial_val,
      nugget_initial = nugget_initial("nugget", nugget = 0.5, known = "given")
    )
    pred_obs <- covmatrix(fit, "CapeHorn", cov_type = "pred.obs")
    obs_pred <- covmatrix(fit, "CapeHorn", cov_type = "obs.pred")
    pred_pred <- covmatrix(fit, "CapeHorn", cov_type = "pred.pred")

    expect_true(all(is.finite(pred_obs)))
    expect_equal(obs_pred, t(pred_obs))
    expect_true(all(is.finite(pred_pred)))
  }
})

test_that("Satterthwaite keeps a free Euclidean extra on the numeric path", {
  skip_if_not_installed("numDeriv")
  copy_lsn_to_temp()
  ssn.object <- ssn_import(file.path(tempdir(), "MiddleFork04.ssn"), overwrite = TRUE)
  fit <- suppressWarnings(ssn_lm(
    Summer_mn ~ ELEV_DEM, ssn.object,
    tailup_type = "none", taildown_type = "none", euclid_type = "matern",
    nugget_type = "nugget",
    euclid_initial = euclid_initial(
      "matern", de = 1, range = 10000, extra = 1, known = c("de", "range")
    ),
    nugget_initial = nugget_initial("nugget", nugget = 0.5, known = "given"),
    control = list(maxit = 15)
  ))
  context <- SSN2:::get_satterthwaite_context(fit, "numeric")
  expect_identical(context$cov_names_free_orig, "euclid_extra")

  covariance <- vcov(fit, type = "cov")
  if (!is.null(covariance)) expect_identical(rownames(covariance), "euclid_extra")
})

# effective range
test_that("closed-form effective ranges agree with spmodel conventions", {
  target <- 0.05
  for (range in c(1e-8, 1, 1e8)) {
    for (constructor in list(tailup_params, taildown_params, euclid_params)) {
      expect_equal(get_effective_range(constructor("exponential", de = 2, range = range)), -log(target) * range)
    }
    expect_equal(get_effective_range(euclid_params("gaussian", de = 2, range = range)), sqrt(-log(target)) * range)
    for (type in c("gravity", "rquad", "magnetic", "cauchy")) {
      p <- switch(type, gravity = 0.5, rquad = 1, magnetic = 1.5, cauchy = 0.7)
      params <- make_euclid_params(type, 2, range, 0, 1, p)
      expect_equal(get_effective_range(params), range * sqrt(target^(-1 / p) - 1))
    }
    params <- euclid_params("pexponential", de = 2, range = range, extra = 1.6)
    expect_equal(get_effective_range(params), (-log(target) * range)^(1 / 1.6))
    expect_equal(get_effective_range(euclid_params("matern", de = 0, range = range, extra = 0.5)),
                 -log(target) * range, tolerance = 1e-8)
  }
})

test_that("stream kernels use flow-connected distances and SSN Mariah scaling", {
  for (constructor in list(tailup_params, taildown_params)) {
    for (type in c("exponential", "mariah", "gaussian")) {
      params <- constructor(type, de = 3, range = 250)
      for (target in c(0.05, 0.2)) {
        distance <- get_effective_range(params, target)
        r <- distance / 250
        correlation <- switch(type,
          exponential = exp(-r), mariah = log1p(90 * r) / (90 * r),
          gaussian = 2 * exp(-r^2) * pnorm(-sqrt(2) * r)
        )
        expect_equal(correlation, target, tolerance = 1e-8)
      }
    }
    for (type in c("linear", "spherical", "epa")) {
      expect_equal(get_effective_range(constructor(type, de = 1, range = 250)), 250)
    }
  }
  for (type in c("spherical", "circular", "cubic", "pentaspherical")) {
    expect_no_warning(expect_equal(get_effective_range(euclid_params(type, de = 1, range = 250)), 250))
  }
})

test_that("range inversion respects powered-exponential and Bessel parameter units", {
  for (type in c("exponential", "gaussian", "matern", "cauchy", "pexponential")) {
    params <- make_euclid_params(type, 0, 1, 0.5, 0.4, 0.7)
    distances <- c(10, 7000)
    ranges <- get_range_from_effective(params, distances)
    result <- vapply(ranges, function(value) {
      params[["range"]] <- value
      get_effective_range(params)
    }, numeric(1))
    expect_equal(result, distances, tolerance = 1e-8)
  }
  wave <- euclid_params("wave", de = 1, range = 7)
  expect_warning(expect_equal(get_effective_range(wave), 140), "not well defined")
  bessel <- euclid_params("jbessel", de = 1, range = 7)
  expect_warning(expect_equal(get_effective_range(bessel), 2 / (pi * 0.05^2 * 7)), "not well defined")
  expect_warning(ranges <- get_range_from_effective(bessel, c(10, 100)), "not well defined")
  expect_equal(ranges[[1]] / ranges[[2]], 10)
})

test_that("effective-range validation rejects undefined input", {
  params <- tailup_params("exponential", de = 1, range = 10)
  for (target in list(NA_real_, NaN, Inf, 0, 1, c(0.05, 0.1), "0.05")) {
    expect_error(get_effective_range(params, target), "target must")
  }
  expect_error(get_effective_range(c(range = 1)), "parameter object")
  expect_error(get_effective_range(tailup_initial("exponential", range = 1)), "parameter object")
  # Bypass constructor validation to exercise the effective-range check itself.
  for (range in c(NA, Inf, -1, 0)) {
    bad_params <- structure(c(de = 1, range = range), class = "tailup_exponential")
    expect_error(get_effective_range(bad_params), "finite and positive")
  }
  expect_equal(get_effective_range(tailup_params("none")), 0)
  expect_equal(get_effective_range(euclid_params("none")), 0)
  expect_equal(get_effective_range(nugget_params("nugget", nugget = 3)), 0)
})

test_that("decorrelation shares stream effective distances and separates Euclidean distances", {
  grid <- SSN2:::ssn_decorrelate_grid_internal(
    Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential",
    taildown_type = "spherical", euclid_type = "gaussian", additive = "afvArea",
    dense_grid = FALSE, add_iid = FALSE, candidates = NULL
  )
  values <- tidy(grid)
  stream <- sort(unique(values$tailup_range * -log(0.05)))
  expect_equal(stream, sort(unique(values$taildown_range)))
  bbox <- sf::st_bbox(mf04p$obs)
  euclid_max <- sqrt((bbox[["xmax"]] - bbox[["xmin"]])^2 + (bbox[["ymax"]] - bbox[["ymin"]])^2)
  expect_equal(sort(unique(values$euclid_range * sqrt(-log(0.05)))), c(0.25, 0.75) * euclid_max)
  expect_false(isTRUE(all.equal(stream, c(0.25, 0.75) * euclid_max)))

  given <- SSN2:::ssn_decorrelate_grid_internal(
    Summer_mn ~ ELEV_DEM, mf04p,
    euclid_type = "pexponential", euclid_params = c(extra = 1.3),
    dense_grid = FALSE, add_iid = FALSE, candidates = NULL
  )
  values <- tidy(given)
  expect_true(all(values$euclid_extra == 1.3))
  expect_equal(sort(unique((-log(0.05) * values$euclid_range)^(1 / 1.3))), c(0.25, 0.75) * euclid_max)
})

