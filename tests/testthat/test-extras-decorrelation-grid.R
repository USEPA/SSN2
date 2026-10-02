skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# decorrelate api
test_that("decorrelation exposes local and resolves conditional defaults internally", {
  for (fun in list(ssn_decorrelate, ssn_decorrelate_data, ssn_decorrelate_newdata,
                   predict.ssn_decorrelate)) {
    expect_false("conditioning" %in% names(formals(fun)))
    expect_identical(formals(fun)["local"], as.list(formals(function(local) NULL)))
  }
  expect_identical(formals(ssn_decorrelate)[c("training", "grid")],
                   as.list(formals(function(training, grid) NULL)))
  for (fun in list(ssn_decorrelate, ssn_decorrelate_grid)) {
    expect_identical(formals(fun)$dense_grid, FALSE)
  }
  expect_identical(formals(ssn_decorrelate)$algorithm, "ranger")
  for (fun in list(ssn_decorrelate, ssn_decorrelate_data)) {
    expect_identical(formals(fun)["ordering"], as.list(formals(function(ordering) NULL)))
  }
  expect_identical(SSN2:::get_decorrelate_ordering(NULL), "pid")
  expect_error(SSN2:::get_decorrelate_ordering("bogus"), "ordering must be")
  expect_identical(get_decorrelate_local(NULL, 5000), list(method = "all"))
  expect_identical(get_decorrelate_local(NULL, 5001), list(method = "covariance", size = 30L))
  expect_identical(get_decorrelate_local(FALSE, 5001), list(method = "all"))
  expect_identical(get_decorrelate_local(list(size = 4), 45), list(method = "covariance", size = 4L))
  expect_identical(get_decorrelate_local(list(method = "covariance"), 45),
                   get_decorrelate_local(TRUE, 45))
  expect_identical(get_decorrelate_local(list(method = "all", size = 4), 45),
                   get_decorrelate_local(FALSE, 45))
  expect_error(get_decorrelate_local(list(size = 1.5), 45), "positive integer")
  expect_error(get_decorrelate_local(list(method = "exact"), 45), "local\\$method")
  expect_error(get_decorrelate_local(list(szie = 4), 45), "named list")
  expect_error(ssn_decorrelate(Summer_mn ~ ELEV_DEM, mf04p, conditioning = FALSE),
               "conditioning is not an argument")
})

test_that("pid, none, and seeded random orderings preserve their meanings", {
  rows <- data.frame(pid = c(3, 1, 2), NetworkID = c("a", "b", "a"), row = 1:3)
  expect_identical(get_decorrelate_order(rows, "pid"), c(2L, 3L, 1L))
  expect_identical(get_decorrelate_order(rows, "none"), 1:3)
  set.seed(2)
  expected <- sample.int(3)
  set.seed(2)
  expect_identical(get_decorrelate_order(rows, "random"), expected)
})

test_that("editable grids and tidied grids round trip without changing parameters", {
  ssn <- mf04p
  ssn$obs$group <- rep(1:3, length.out = nrow(ssn$obs))
  grid <- ssn_decorrelate_grid(
    Summer_mn ~ ELEV_DEM, ssn, tailup_type = "exponential",
    euclid_type = "cauchy", additive = "afvArea", anisotropy = TRUE,
    random = ~ group, dense_grid = FALSE
  )
  expect_s3_class(grid, "data.frame")
  expect_false("candidate" %in% names(grid))
  expect_equal(tidy(SSN2:::get_decorrelate_grid_table(grid)), tidy(grid))
  expect_equal(tidy(SSN2:::get_decorrelate_grid_table(tidy(grid))), tidy(grid))
  invalid <- grid[1, , drop = FALSE]
  invalid$tailup_range <- NULL
  expect_error(SSN2:::get_decorrelate_grid_table(invalid), "tailup_range")
  expect_error(SSN2:::get_decorrelate_grid_table(grid[FALSE, ]), "grid must have rows")
})

test_that("an edited single-row grid is evaluated and overrides covariance types", {
  testthat::local_mocked_bindings(
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) mean(y),
    predict_decorrelate_algorithm = function(fit, X, algorithm) rep(fit, NROW(X))
  )
  grid <- SSN2:::ssn_decorrelate_grid_internal(
    Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea",
    dense_grid = FALSE, add_iid = FALSE, candidates = NULL
  )[1, , drop = FALSE]
  grid$tailup_range <- 1234
  fit <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "spherical",
    grid = grid, additive = "afvArea", local = list(size = 4),
    training = list(training_index = 1:30, test_index = 31:45)
  )
  expect_equal(nrow(fit$grid), 1L)
  expect_equal(fit$grid$tailup_range, 1234)
  expect_equal(fit$grid$tailup_type, "exponential")
  expect_equal(unname(fit$decorrelate_data$params$tailup["range"]), 1234)
  expect_identical(fit$decorrelate_data$local, list(method = "covariance", size = 4L))
  expect_true(is.finite(fit$test$RMSPE))
  expect_length(predict(fit, "CapeHorn", local = FALSE), nrow(mf04p$preds$CapeHorn))
})

test_that("prediction local controls can override the observed transformation", {
  args <- list(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p, additive = "afvArea",
    tailup_params = tailup_params("exponential", de = 1, range = 10000),
    nugget_params = nugget_params("nugget", nugget = 0.1),
    ordering = "none"
  )
  exact <- do.call(ssn_decorrelate_data, c(args, list(local = FALSE)))
  approximate <- do.call(ssn_decorrelate_data, c(args, list(local = list(size = 4))))
  expected <- ssn_decorrelate_newdata(exact, "CapeHorn")
  overridden <- ssn_decorrelate_newdata(approximate, "CapeHorn", local = FALSE)
  expect_equal(overridden$tX_newdata, expected$tX_newdata)
  expect_equal(overridden$yoffset, expected$yoffset)
  expect_equal(overridden$yscale, expected$yscale)
  expected_local <- ssn_decorrelate_newdata(approximate, "CapeHorn")
  overridden_local <- ssn_decorrelate_newdata(exact, "CapeHorn", local = list(size = 4))
  expect_equal(overridden_local$tX_newdata, expected_local$tX_newdata)
  expect_equal(overridden_local$yoffset, expected_local$yoffset)
  expect_equal(overridden_local$yscale, expected_local$yscale)
  expect_equal(ssn_decorrelate_newdata(exact, "CapeHorn", local = list(method = "all"))$tX_newdata,
               expected$tX_newdata)
})

# decorrelate grid density
test_that("large decorrelation grids use observed geometry for stream ranges", {
  large <- mf04p
  large$obs <- large$obs[rep(1:2, length.out = 5002L), ]
  points <- rep(list(sf::st_point(c(0, 0)), sf::st_point(c(300, 400))), 2501L)
  points[[5002L]] <- sf::st_point(c(1e6, 1e6))
  sf::st_geometry(large$obs) <- sf::st_sfc(points, crs = sf::st_crs(mf04p$obs))
  large$obs$Summer_mn[5002L] <- NA_real_
  testthat::local_mocked_bindings(
    get_dist_object = function(...) stop("Large grids must not load full distances")
  )
  defaults <- ssn_decorrelate_grid(Summer_mn ~ 1, large,
    tailup_type = "exponential", taildown_type = "exponential",
    euclid_type = "exponential", additive = "afvArea")
  sparse <- ssn_decorrelate_grid(Summer_mn ~ 1, large,
    tailup_type = "exponential", taildown_type = "exponential",
    euclid_type = "exponential", additive = "afvArea", dense_grid = FALSE)
  expect_equal(defaults, sparse)
  expect_equal(nrow(defaults), 19L)
  for (dense in c(FALSE, TRUE)) {
    grid <- SSN2:::ssn_decorrelate_grid_internal(
      Summer_mn ~ 1, large, tailup_type = "exponential",
      taildown_type = "exponential", euclid_type = "exponential",
      additive = "afvArea", dense_grid = dense, add_iid = FALSE, candidates = NULL
    )
    for (component in c("tailup", "taildown", "euclid")) {
      ranges <- sort(unique(grid[[paste0(component, "_range")]]))
      expect_equal(-log(0.05) * ranges, c(125, 375))
    }
  }
})

test_that("both grid densities expand over all active SSN components", {
  for (dense in c(FALSE, TRUE)) {
    expected <- if (dense) c(11L, 65L, 137L) else c(5L, 15L, 19L)
    for (nspatial in 1:3) {
      grid <- ssn_decorrelate_grid(
        Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential",
        taildown_type = if (nspatial >= 2) "exponential" else "none",
        euclid_type = if (nspatial == 3) "exponential" else "none",
        additive = "afvArea", dense_grid = dense
      )
      values <- tidy(grid)
      expect_equal(NROW(values), expected[[nspatial]])
      expect_equal(sum(values$tailup_type == "no transformation"), 1L)
      spatial <- values[values$tailup_type != "no transformation", ]
      active <- c("tailup", "taildown", "euclid")[seq_len(nspatial)]
      range_columns <- paste0(active, "_range")
      expect_equal(NROW(unique(spatial[range_columns])), if (dense) 2^nspatial else 2L)
      for (column in range_columns) expect_length(unique(spatial[[column]]), 2L)
      variance_columns <- c(paste0(active, "_de"), "nugget_nugget")
      totals <- rowSums(spatial[variance_columns])
      expect_equal(totals, rep(totals[[1]], length(totals)))
      if (nspatial == 3L) {
        proportions <- as.matrix(spatial[variance_columns]) / totals
        for (column in seq_along(variance_columns)) {
          expect_true(any(abs(proportions[, column] - 0.95) < 1e-12))
        }
      }
    }
  }
})

test_that("dense joint grids retain the non-dense candidates", {
  args <- list(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", taildown_type = "exponential",
    euclid_type = "exponential", additive = "afvArea"
  )
  dense <- as.data.frame(tidy(do.call(ssn_decorrelate_grid, c(args, list(dense_grid = TRUE)))))
  sparse <- as.data.frame(tidy(do.call(ssn_decorrelate_grid, c(args, list(dense_grid = FALSE)))))
  parameters <- names(dense)
  expect_equal(NROW(merge(dense[parameters], sparse[parameters])), NROW(sparse))
  expect_gt(NROW(dense), NROW(sparse))
  defaults <- tidy(do.call(ssn_decorrelate_grid, args))
  expect_equal(defaults, tidy(do.call(ssn_decorrelate_grid, c(args, list(dense_grid = FALSE)))))
})

test_that("omitted and explicit sparse defaults evaluate the same candidates", {
  testthat::local_mocked_bindings(
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) mean(y),
    predict_decorrelate_algorithm = function(fit, X, algorithm) rep(fit, NROW(X))
  )
  args <- list(formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", additive = "afvArea",
    training = list(training_index = 1:30, test_index = 31:45))
  defaults <- do.call(ssn_decorrelate, args)
  sparse <- do.call(ssn_decorrelate, c(args, list(dense_grid = FALSE)))
  expect_equal(defaults$grid, sparse$grid)
  expect_equal(defaults$test, sparse$test)
  expect_equal(nrow(defaults$grid), 5L)
})

test_that("random-effect grid allocations retain total variance and short random-dominant ranges", {
  ssn <- mf04p
  ssn$obs$group1 <- rep(1:3, length.out = nrow(ssn$obs))
  ssn$obs$group2 <- rep(1:5, length.out = nrow(ssn$obs))
  args <- list(formula = Summer_mn ~ ELEV_DEM, ssn.object = ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    additive = "afvArea", dense_grid = FALSE)
  dense <- do.call(ssn_decorrelate_grid,
    utils::modifyList(args, list(random = ~ group1, dense_grid = TRUE)))
  expect_equal(nrow(dense), 146L)
  ns2 <- 1.2 * sum(lm.fit(model.matrix(Summer_mn ~ ELEV_DEM, ssn$obs), ssn$obs$Summer_mn)$residuals^2) /
    (nrow(ssn$obs) - 2L)
  for (k in 1:2) {
    random <- if (k == 1L) ~ group1 else ~ group1 + group2
    g <- as.data.frame(do.call(ssn_decorrelate_grid, c(args, list(random = random))))
    expect_equal(nrow(g), if (k == 1L) 22L else 24L)
    core <- g$tailup_type != "none"
    spatial_names <- c("tailup_de", "taildown_de", "euclid_de", "nugget_nugget")
    random_names <- names(g)[startsWith(names(g), "randcov_")]
    expect_length(random_names, k)
    expect_equal(unname(rowSums(g[core, c(spatial_names, random_names)])), rep(ns2, sum(core)))
    random_share <- rowSums(g[core, random_names, drop = FALSE]) / ns2
    expect_equal(sort(unique(round(random_share, 10))), round(sort(c(.1, k / (4 + k), .9)), 10))
    random_rows <- which(core)[abs(random_share - .9) < 1e-10]
    expect_length(random_rows, if (k == 1L) 1L else 3L)
    for (name in c("tailup_range", "taildown_range", "euclid_range")) {
      expect_true(all(g[random_rows, name] == min(g[core, name])))
    }
    expect_true(all(g[!core, random_names, drop = FALSE] == 0))
  }
  for (setting in list(list(anisotropy = TRUE), list(euclid_type = "matern"),
      list(anisotropy = TRUE, random = ~ group1))) {
    a <- utils::modifyList(args, setting)
    g <- as.data.frame(do.call(ssn_decorrelate_grid, a))
    expected <- if (!is.null(setting$random)) 64L else if (!is.null(setting$euclid_type)) 37L else 55L
    expect_equal(nrow(g), expected)
    expect_equal(nrow(unique(g)), nrow(g))
    expect_true(all(g$euclid_rotate[g$euclid_scale == 1] == 0))
    dense <- as.data.frame(do.call(ssn_decorrelate_grid, utils::modifyList(a, list(dense_grid = TRUE))))
    expect_equal(nrow(merge(g, dense)), nrow(g))
  }
})

test_that("expanded compact grids preserve supplied random and spatial parameters", {
  ssn <- mf04p
  ssn$obs$group <- rep(1:3, length.out = nrow(ssn$obs))
  args <- list(formula = Summer_mn ~ ELEV_DEM, ssn.object = ssn,
    tailup_type = "exponential", tailup_params = c(range = 7000),
    euclid_type = "matern", euclid_params = c(extra = 1, rotate = .2, scale = .7),
    additive = "afvArea", anisotropy = TRUE, random = ~ group,
    randcov_params = randcov_params(group = 2))
  set.seed(2)
  before <- .Random.seed
  g <- as.data.frame(do.call(ssn_decorrelate_grid, args))
  expect_identical(.Random.seed, before)
  g <- g[g$tailup_type != "none", ]
  expect_true(all(g$tailup_range == 7000))
  expect_true(all(g$euclid_extra == 1 & g$euclid_rotate == .2 & g$euclid_scale == .7))
  expect_true(all(g[["randcov_1 | group"]] == 2))
  expect_equal(nrow(unique(g)), nrow(g))
})

test_that("compact extension counts have fewer active components", {
  ssn <- mf04p
  ssn$obs$group <- factor(rep(1:3, length.out = nrow(ssn$obs)))
  for (nspatial in 1:3) {
    g <- ssn_decorrelate_grid(Summer_mn ~ ELEV_DEM, ssn,
      tailup_type = "exponential", taildown_type = if (nspatial > 1) "exponential" else "none",
      euclid_type = if (nspatial > 2) "exponential" else "none",
      additive = "afvArea", random = ~ group)
    expect_equal(nrow(g), c(8L, 18L, 22L)[nspatial])
  }
  g <- ssn_decorrelate_grid(Summer_mn ~ ELEV_DEM, ssn, random = ~ group)
  expect_equal(nrow(g), 4L)
  g <- ssn_decorrelate_grid(Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential", additive = "afvArea", nugget_type = "none", random = ~ group)
  expect_equal(nrow(g), 6L)
})

# decorrelate output
test_that("grid output hides identifiers and formats the untransformed baseline", {
  grid <- ssn_decorrelate_grid(
    Summer_mn ~ ELEV_DEM, mf04p, euclid_type = "exponential", dense_grid = FALSE
  )
  values <- tidy(grid)
  expect_false("candidate" %in% names(values))
  expect_false(any(c("euclid_rotate", "euclid_scale") %in% names(values)))
  baseline <- values[values$euclid_type == "no transformation", ]
  expect_equal(NROW(baseline), 1L)
  expect_true(all(is.na(baseline[c("euclid_de", "euclid_range", "nugget_nugget")])))
  printed <- capture.output(returned <- withVisible(print(grid)))
  expect_false(any(grepl("candidate|\\$grid|\\$iid", printed)))
  expect_true(any(grepl("euclid_type", printed)))
  expect_identical(returned$value, grid)
  expect_false(returned$visible)
  expect_equal(tidy(grid, sort_by = "euclid_range")$euclid_range,
               sort(values$euclid_range, na.last = TRUE))
})

test_that("evaluated grid output retains the winning parameters and supports sorting", {
  testthat::local_mocked_bindings(
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) mean(y),
    predict_decorrelate_algorithm = function(fit, X, algorithm) rep(fit, NROW(X))
  )
  fit <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea",
    training = list(training_index = 1:30, test_index = 31:45), dense_grid = FALSE
  )
  expect_false("candidate" %in% names(fit$grid))
  expect_false("candidate" %in% names(tidy(fit)))
  expect_equal(tidy(fit), tidy(fit$grid))
  expect_equal(fit$test$RMSPE, min(fit$grid$RMSPE))
  expect_equal(remove_covtype(class(fit$decorrelate_data$params$tailup)),
               fit$grid$tailup_type[[1]])
  expect_equal(unname(fit$decorrelate_data$params$nugget), fit$grid$nugget_nugget[[1]],
               ignore_attr = TRUE)
  expect_equal(tidy(fit, sort_by = "RMSPE", decreasing = TRUE)$RMSPE,
               sort(fit$grid$RMSPE, decreasing = TRUE))
  expect_equal(tidy(fit, sort_by = "bias")$bias, sort(fit$grid$bias))
  expect_error(tidy(fit, sort_by = "candidate"), "sort_by must be a variable in x.", fixed = TRUE)
  expect_error(tidy(fit, decreasing = NA), "decreasing must be TRUE or FALSE.", fixed = TRUE)
  expect_error(tidy(fit, decreasing = 1), "decreasing must be TRUE or FALSE.", fixed = TRUE)
  expect_false(any(grepl("candidate", capture.output(print(fit$grid)))))
})

test_that("grid sorting defaults to its stored evaluation statistic", {
  grid <- structure(data.frame(
    tailup_type = rep("exponential", 3), tailup_de = 1,
    taildown_type = "none", euclid_type = "none", nugget_nugget = 0.1,
    bias = c(-2, 1, 0), MSPE = c(4, 1, 9),
    RMSPE = c(2, 1, 3), cor2 = c(0.5, 0.9, 0.2)
  ), class = c("ssn_decorrelate_grid", "data.frame"))
  for (statistic in c("bias", "MSPE", "RMSPE", "cor2")) {
    attr(grid, "statistic") <- statistic
    expected <- sort(grid[[statistic]], decreasing = statistic == "cor2")
    expect_equal(tidy(grid)[[statistic]], expected)
    expect_equal(tidy(grid, sort_by = NULL)$RMSPE, grid$RMSPE)
  }
  attr(grid, "statistic") <- NULL
  expect_equal(tidy(grid)$RMSPE, grid$RMSPE)
  expect_equal(tidy(grid, sort_by = "RMSPE")$RMSPE, sort(grid$RMSPE))
})

test_that("model printing uses spmodel headings and selected SSN covariance parameters", {
  params <- list(
    tailup = tailup_params("exponential", de = 1.234567, range = 5000),
    taildown = taildown_params("spherical", de = 2, range = 6000),
    euclid = euclid_params("exponential", de = 3, range = 7000),
    nugget = nugget_params("nugget", nugget = 0.5),
    randcov = c(group = 0.25)
  )
  fit <- structure(list(
    call = quote(ssn_decorrelate(y ~ x, stream)),
    test = list(bias = 0.1, MSPE = 4, RMSPE = 2, cor2 = 0.8),
    decorrelate_data = list(params = params, covariance_fit = list(anisotropy = FALSE))
  ), class = "ssn_decorrelate")
  printed <- capture.output(returned <- withVisible(print(fit, digits = 3)))
  headings <- printed[grepl("^Call:|^stats:|^Coefficients", printed)]
  expect_identical(headings, c(
    "Call:", "stats:", "Coefficients (exponential tailup covariance):",
    "Coefficients (spherical taildown covariance):",
    "Coefficients (exponential Euclidean covariance):",
    "Coefficients (nugget):", "Coefficients (random effects):"
  ))
  expect_true(any(grepl("1.23", printed, fixed = TRUE)))
  expect_false(any(grepl("rotate|scale|Machine-learning algorithm|Holdout statistics", printed)))
  expect_identical(returned$value, fit)
  expect_false(returned$visible)
  fit$decorrelate_data$covariance_fit$anisotropy <- TRUE
  expect_output(print(fit), "rotate")
  fit$test <- NULL
  fit$decorrelate_data$params$tailup <- tailup_params("none")
  fit$decorrelate_data$params$randcov <- NULL
  printed <- capture.output(print(fit))
  expect_false(any(grepl("stats:|tailup covariance|random effects", printed)))
  expect_true(any(grepl("Euclidean covariance", printed)))
})

