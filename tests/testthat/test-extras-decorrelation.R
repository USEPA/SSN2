skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# decorrelate
decorrelate_initials <- function(range = 10000) {
  list(
    tailup_initial = tailup_initial("exponential", de = 1, range = range, known = "given"),
    taildown_initial = taildown_initial("none"),
    euclid_initial = euclid_initial("none"),
    nugget_initial = nugget_initial("nugget", nugget = 0.1, known = "given")
  )
}

# ssn_decorrelate_data()'s known-parameter analog of decorrelate_initials(),
# used wherever a *_params() object (not a *_initial()) is required
decorrelate_params <- function(range = 10000) {
  list(
    tailup_params = tailup_params("exponential", de = 1, range = range),
    taildown_params = taildown_params("none"),
    euclid_params = euclid_params("none"),
    nugget_params = nugget_params("nugget", nugget = 0.1)
  )
}

# Converts one ssn_decorrelate_grid()/get_decorrelate_grid_candidates()
# candidate's *_initial()-shaped elements into ssn_decorrelate_data()'s
# tailup_params/taildown_params/euclid_params/nugget_params + randcov_initial,
# mirroring the package's own internal get_decorrelate_data_args() conversion
candidate_to_params_args <- function(candidate) {
  list(
    tailup_params = SSN2:::get_params_from_initial(candidate$tailup_initial, tailup_params),
    taildown_params = SSN2:::get_params_from_initial(candidate$taildown_initial, taildown_params),
    euclid_params = SSN2:::get_params_from_initial(candidate$euclid_initial, euclid_params),
    nugget_params = SSN2:::get_params_from_initial(candidate$nugget_initial, nugget_params),
    randcov_params = SSN2:::get_randcov_params_from_initial(candidate$randcov_initial)
  )
}

test_that("default decorrelation evaluates a parameter grid and an untransformed baseline", {
  testthat::local_mocked_bindings(
    ssn_lm = function(...) stop("No covariance model may be fitted"),
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) mean(y),
    predict_decorrelate_algorithm = function(fit, X, algorithm) rep(fit, NROW(X))
  )
  fit <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea",
    training = list(training_index = 1:30, test_index = 31:45)
  )
  grid <- tidy(fit)
  expect_equal(NROW(grid), 5L)
  expect_equal(sum(grid$tailup_type == "no transformation"), 1L)
  expect_length(unique(grid$tailup_range[grid$tailup_type == "exponential"]), 2L)
  expect_length(unique(grid$tailup_de[grid$tailup_type == "exponential"]), 2L)
  expect_true(all(is.finite(grid$RMSPE)))
  expected <- sqrt(mean((mean(mf04p$obs$Summer_mn[1:30]) - mf04p$obs$Summer_mn[31:45])^2))
  expect_equal(grid$RMSPE[grid$tailup_type == "no transformation"], expected)

  sparse <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea",
    training = list(training_index = 1:30, test_index = 31:45),
    dense_grid = FALSE, local = list(method = "covariance", size = 4),
    evaluate_test = FALSE
  )
  expect_equal(NROW(tidy(sparse)), 5L)
  expect_equal(sparse$decorrelate_data$local$size, 4L)
})

test_that("generated grids preserve supplied parameters and include a true IID transformation", {
  grid <- ssn_decorrelate_grid(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p,
    tailup_type = "exponential", tailup_params = c(range = 7000),
    additive = "afvArea"
  )
  expect_equal(NROW(tidy(grid)), 3L)
  expect_true(all(tidy(grid)$tailup_range[tidy(grid)$tailup_type == "exponential"] == 7000))
  iid <- do.call(ssn_decorrelate_data, c(list(Summer_mn ~ ELEV_DEM, mf04p),
    candidate_to_params_args(get_decorrelate_grid_candidates(grid)[[nrow(grid)]])))
  expect_equal(as.numeric(iid$tX), as.numeric(iid$X))
  expect_equal(iid$ty, iid$y)
  prediction <- ssn_decorrelate_newdata(iid, "CapeHorn")
  # tX_newdata is a computed transform, not a raw model matrix, so (unlike
  # X_newdata) it carries no "assign" attribute
  expect_equal(unname(prediction$tX_newdata), unname(prediction$X_newdata), ignore_attr = "assign")
  expect_equal(prediction$yoffset, rep(0, NROW(prediction$X_newdata)))
  expect_equal(prediction$yscale, rep(1, NROW(prediction$X_newdata)))
  expect_error(ssn_decorrelate_grid(formula = Summer_mn ~ ELEV_DEM, ssn.object = mf04p, dense_grid = NA), "dense_grid")
})

test_that("automatic grids support joint covariance, anisotropy, shape, and random effects", {
  ssn <- mf04p
  ssn$obs$decorr_group <- factor(rep(1:3, length.out = NROW(ssn$obs)))
  grid <- ssn_decorrelate_grid(
    formula = Summer_mn ~ ELEV_DEM, ssn.object = ssn,
    tailup_type = "exponential", taildown_type = "spherical", euclid_type = "cauchy",
    additive = "afvArea", anisotropy = TRUE, random = ~ decorr_group, dense_grid = FALSE
  )
  parameters <- tidy(grid)
  spatial <- parameters[parameters$tailup_type != "no transformation", ]
  expect_true(all(is.finite(spatial$euclid_range)))
  expect_equal(sort(unique(spatial$euclid_extra)), c(0.5, 2))
  expect_true(all(spatial$euclid_rotate[spatial$euclid_scale == 1] == 0))
  expect_true(all(spatial[["randcov_1 | decorr_group"]] > 0))
  # random and anisotropy are no longer ssn_decorrelate_data() arguments --
  # both are rederived from randcov_params/euclid_params, which
  # candidate_to_params_args() already supplies from the candidate itself
  iid <- do.call(ssn_decorrelate_data, c(
    list(Summer_mn ~ ELEV_DEM, ssn),
    candidate_to_params_args(get_decorrelate_grid_candidates(grid)[[nrow(grid)]])
  ))
  expect_equal(as.numeric(iid$tX), as.numeric(iid$X))
  expect_equal(iid$ty, iid$y)
  candidate <- do.call(ssn_decorrelate_data, c(
    list(Summer_mn ~ ELEV_DEM, ssn, additive = "afvArea"),
    candidate_to_params_args(get_decorrelate_grid_candidates(grid)[[1]])
  ))
  expect_true(all(is.finite(candidate$ty)))
})

test_that("large automatic grids avoid full distance matrices", {
  testthat::local_mocked_bindings(
    get_dist_object = function(...) stop("A full distance matrix is not allowed")
  )
  ssn <- mf04p
  ssn$obs <- ssn$obs[rep(seq_len(NROW(ssn$obs)), length.out = 5001L), ]
  grid <- ssn_decorrelate_grid(
    Summer_mn ~ ELEV_DEM, ssn, tailup_type = "exponential", additive = "afvArea"
  )
  expect_equal(nrow(grid), 5L)
})

test_that("partially supplied random variances still require a grid search", {
  ssn <- mf04p
  ssn$obs$decorr_group <- factor(rep(1:3, length.out = NROW(ssn$obs)))
  ssn$obs$decorr_group2 <- factor(rep(1:4, length.out = NROW(ssn$obs)))
  initial <- get_initial_object(
    "none", "none", "none", "nugget", NULL, NULL, NULL,
    nugget_initial("nugget", nugget = 1, known = "given")
  )
  random <- ~ decorr_group + decorr_group2
  randcov <- randcov_initial(decorr_group = 0.2, known = "given")
  expect_false(get_decorrelate_covariance_known(initial, random, randcov))
  grid <- SSN2:::ssn_decorrelate_grid_internal(
    Summer_mn ~ ELEV_DEM, ssn, nugget_params = nugget_params("nugget", nugget = 1),
    random = random, randcov_params = randcov_params(decorr_group = 0.2),
    add_iid = FALSE, candidates = NULL
  )
  values <- tidy(grid)
  expect_true(all(values[["randcov_1 | decorr_group"]] == 0.2))
  expect_gt(length(unique(values[["randcov_1 | decorr_group2"]])), 1L)
})

fit_decorrelate_data <- function(ssn, ...) {
  params <- decorrelate_params()
  ssn_decorrelate_data(
    Summer_mn ~ ELEV_DEM, ssn, additive = "afvArea",
    tailup_params = params$tailup_params,
    taildown_params = params$taildown_params,
    euclid_params = params$euclid_params,
    nugget_params = params$nugget_params,
    ...
  )
}

test_that("exact stream decorrelation reproduces GLS quadratic forms", {
  # method = "all" conditions each observation on every earlier-ordered
  # observation sequentially (matching spmodel's decorrelate_data_internal()/
  # get_decorrelated_value()); a valid whitening transform tX = L^-1 P X (any
  # permutation P, any conditioning order) must satisfy tX'tX = X'Sigma^-1 X
  # and likewise for ty, regardless of the internal ordering, so this checks
  # the transform algebraically without needing to inspect any per-row
  # conditioning detail
  transformed <- fit_decorrelate_data(mf04p)
  expect_s3_class(transformed, "ssn_decorrelate_data")

  Sigma <- covmatrix(transformed$covariance_fit) / transformed$total_var
  X <- transformed$X
  y <- transformed$y
  tX <- transformed$tX
  ty <- transformed$ty

  expect_equal(unname(crossprod(tX)), unname(crossprod(X, solve(Sigma, X))), tolerance = 1e-8)
  expect_equal(
    unname(as.numeric(crossprod(tX, ty))), unname(as.numeric(crossprod(X, solve(Sigma, y)))),
    tolerance = 1e-8
  )
  expect_equal(as.numeric(crossprod(ty)), as.numeric(crossprod(y, solve(Sigma, y))), tolerance = 1e-8)

  prediction <- ssn_decorrelate_newdata(transformed, "CapeHorn")
  recolored <- ssn_recorrelate_newdata(prediction, rep(0, NROW(prediction$tX_newdata)))
  expect_equal(unname(recolored), prediction$yoffset + prediction$offset, tolerance = 1e-12)
})

test_that("every ordering value produces a valid permutation and reproduces GLS quadratic forms", {
  # GLS quadratic forms are invariant to a consistent row/column permutation.
  Sigma <- NULL
  orderings <- c("pid", "none", "random", "maxmin", "middleout", "outsidein", "coordinate", "grts")
  for (ord in orderings) {
    transformed <- fit_decorrelate_data(mf04p, ordering = ord)
    if (is.null(Sigma)) Sigma <- covmatrix(transformed$covariance_fit) / transformed$total_var
    X <- transformed$X
    tX <- transformed$tX
    expect_equal(
      unname(crossprod(tX)), unname(crossprod(X, solve(Sigma, X))),
      tolerance = 1e-6, info = ord
    )
  }
  expect_error(fit_decorrelate_data(mf04p, ordering = "bogus"), "ordering must be")
})

test_that("decorrelation applies a formula offset exactly once", {
  ssn <- mf04p
  ssn$obs$decorr_offset <- seq_len(NROW(ssn$obs)) / 100
  ssn$preds$CapeHorn$decorr_offset <- seq_len(NROW(ssn$preds$CapeHorn)) / 100
  params <- decorrelate_params()
  transformed <- ssn_decorrelate_data(
    Summer_mn ~ ELEV_DEM + offset(decorr_offset), ssn, additive = "afvArea",
    tailup_params = params$tailup_params,
    taildown_params = params$taildown_params,
    euclid_params = params$euclid_params,
    nugget_params = params$nugget_params
  )
  prediction <- ssn_decorrelate_newdata(transformed, "CapeHorn")
  expect_equal(
    unname(ssn_recorrelate_newdata(prediction, rep(0, NROW(prediction$tX_newdata)))),
    prediction$yoffset + prediction$offset,
    tolerance = 1e-12
  )
})

test_that("decorrelation predicts response-missing observations without among-prediction distances", {
  ssn <- mf04p
  ssn$obs$decorrelation_response <- ssn$obs$Summer_mn
  ssn$obs$decorrelation_response[c(3, 7)] <- NA_real_
  params <- decorrelate_params()
  transformed <- ssn_decorrelate_data(
    decorrelation_response ~ ELEV_DEM, ssn, additive = "afvArea",
    tailup_params = params$tailup_params,
    taildown_params = params$taildown_params,
    euclid_params = params$euclid_params,
    nugget_params = params$nugget_params
  )
  prediction <- ssn_decorrelate_newdata(transformed, ".missing")
  expect_equal(NROW(prediction$tX_newdata), 2)
  expect_true(all(is.finite(ssn_recorrelate_newdata(prediction, rep(0, 2)))))
})

test_that("covariance graph converges to the exact transform", {
  ssn <- mf04p
  ssn_create_bigdist(ssn, predpts = "CapeHorn", overwrite = TRUE, no_cores = 1, verbose = FALSE)
  exact <- fit_decorrelate_data(ssn)
  local <- fit_decorrelate_data(
    ssn,
    local = list(method = "covariance", size = NROW(ssn$obs))
  )

  expect_equal(unname(local$tX), unname(exact$tX), tolerance = 1e-8)
  expect_equal(local$ty, exact$ty, tolerance = 1e-8)

  local$covariance_fit$ssn.object$preds$CapeHorn <-
    local$covariance_fit$ssn.object$preds$CapeHorn[seq_len(3), , drop = FALSE]
  prediction <- ssn_decorrelate_newdata(local, "CapeHorn")
  expect_true(all(is.finite(ssn_recorrelate_newdata(prediction, rep(0, 3)))))
})

test_that("learner wrapper evaluates joint covariance candidates and predicts", {
  skip_if_not_installed("ranger")
  candidates <- SSN2:::get_decorrelate_grid_table(list(
    short = decorrelate_initials(5000),
    long = decorrelate_initials(10000)
  ))
  fit <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    grid = candidates, training = list(test_index = seq_len(10)),
    num.trees = 20, seed = 2, num.threads = 1
  )

  expect_s3_class(fit, "ssn_decorrelate")
  expect_s3_class(fit$grid, "ssn_decorrelate_grid")
  expect_equal(NROW(fit$grid), 2)
  expect_equal(length(predict(fit, "CapeHorn")), NROW(mf04p$preds$CapeHorn))
  expect_s3_class(tidy(fit), "tbl")
  expect_output(print(fit), "stats:")
})

test_that("recoloring treats ty_newdata purely positionally, ignoring any names", {
  transformed <- fit_decorrelate_data(mf04p)
  prediction <- ssn_decorrelate_newdata(transformed, "CapeHorn")
  values <- seq_len(NROW(prediction$tX_newdata))

  unnamed_result <- ssn_recorrelate_newdata(prediction, values)
  expect_equal(
    unnamed_result,
    prediction$yscale * values + prediction$yoffset + prediction$offset
  )

  # names -- whether valid row identities, garbage, or duplicated -- are
  # silently ignored; there is no longer a rows$key validation/reorder step
  named_result <- ssn_recorrelate_newdata(prediction, setNames(values, paste0("x", values)))
  expect_equal(named_result, unnamed_result)

  # a genuinely reordered vector is applied positionally rather than
  # realigned -- the accepted tradeoff of matching spmodel's own
  # recorrelate_newdata(), which is likewise purely positional
  reversed_result <- ssn_recorrelate_newdata(prediction, rev(values))
  expect_equal(
    reversed_result,
    prediction$yscale * rev(values) + prediction$yoffset + prediction$offset
  )
  expect_false(isTRUE(all.equal(reversed_result, unnamed_result)))
})

test_that("decorrelation training uses exact fold fields and original row indices", {
  response_index <- setdiff(seq_len(12), c(2, 9))
  folds <- rep(c("upstream", "downstream"), length.out = 12)
  set.seed(2)
  rng_before <- .Random.seed
  training <- SSN2:::get_decorrelate_training(
    list(method = "cv", folds = 5, folds_index = folds),
    response_index, 12
  )

  expect_identical(.Random.seed, rng_before)
  expect_identical(training$folds_index, folds[response_index])
  expect_identical(
    training$splits$upstream$test_index,
    which(folds[response_index] == "upstream")
  )

  alias <- SSN2:::get_decorrelate_training(
    list(method = "cv", fold = folds), response_index, 12
  )
  expect_identical(alias$folds_index, training$folds_index)
  expect_error(
    SSN2:::get_decorrelate_training(
      list(method = "cv", folds_index = folds, fold = folds), response_index, 12
    ),
    "Supply only one"
  )

  split <- SSN2:::get_decorrelate_training(
    list(method = "split", test_index = c(11, 12)), response_index, 12
  )
  expect_identical(split$splits[[1]]$test_index, match(c(11L, 12L), response_index))
  expect_error(
    SSN2:::get_decorrelate_training(
      list(method = "split", test_index = 2), response_index, 12
    ),
    "response is missing"
  )
})

test_that("learner evaluation honors an explicit non-exhaustive training subset", {
  ssn <- mf04p
  ssn$obs$decorrelation_response <- ssn$obs$Summer_mn
  ssn$obs$decorrelation_response[[1]] <- NA_real_
  training_rows <- 2:11
  test_rows <- 44:45
  unused_rows <- setdiff(
    which(!is.na(ssn$obs$decorrelation_response)),
    c(training_rows, test_rows)
  )
  ssn$obs$evaluation_group <- rep(c("a", "b"), length.out = NROW(ssn$obs))
  ssn$obs$evaluation_group[unused_rows] <- "unused"
  ssn$obs$evaluation_group <- factor(ssn$obs$evaluation_group)
  fitted_responses <- list()
  testthat::local_mocked_bindings(
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) {
      fitted_responses[[length(fitted_responses) + 1L]] <<- y
      list(mean = mean(y))
    },
    predict_decorrelate_algorithm = function(fit, X, algorithm) {
      rep(fit$mean, NROW(X))
    },
    .package = "SSN2"
  )
  fit_once <- function(object) {
    ssn_decorrelate(
      decorrelation_response ~ ELEV_DEM + evaluation_group, object,
      nugget_params = nugget_params("nugget", nugget = 1),
      training = list(
        training_index = training_rows,
        test_index = test_rows
      ),
      evaluate_test = TRUE
    )
  }

  fit <- fit_once(ssn)
  first_evaluation_response <- fitted_responses[[1]]
  expect_length(first_evaluation_response, length(training_rows))
  expect_equal(first_evaluation_response, ssn$obs$decorrelation_response[training_rows])
  expect_length(fitted_responses[[2]], sum(!is.na(ssn$obs$decorrelation_response)))

  ssn_changed <- ssn
  ssn_changed$obs$decorrelation_response[unused_rows] <-
    ssn_changed$obs$decorrelation_response[unused_rows] + 1000
  fit_changed <- fit_once(ssn_changed)
  expect_equal(
    unclass(fit$grid[c("bias", "MSPE", "RMSPE", "cor2")]),
    unclass(fit_changed$grid[c("bias", "MSPE", "RMSPE", "cor2")])
  )

  fit_local <- ssn_decorrelate(
    decorrelation_response ~ ELEV_DEM + evaluation_group, ssn,
    nugget_params = nugget_params("nugget", nugget = 1),
    training = list(
      training_index = training_rows,
      test_index = test_rows
    ),
    evaluate_test = TRUE,
    local = list(method = "covariance", size = 4)
  )
  expect_s3_class(fit_local, "ssn_decorrelate")
  expect_length(fitted_responses[[5]], length(training_rows))
  expect_equal(fitted_responses[[5]], ssn$obs$decorrelation_response[training_rows])
})

test_that("learner evaluation supports dot formulas on SSN observations", {
  ssn <- mf04p
  keep <- c(
    "Summer_mn", "ELEV_DEM", "netgeom", "netID", "rid", "upDist",
    "ratio", "pid", "locID"
  )
  ssn$obs <- ssn$obs[, intersect(keep, names(ssn$obs))]
  testthat::local_mocked_bindings(
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) {
      list(mean = mean(y))
    },
    predict_decorrelate_algorithm = function(fit, X, algorithm) {
      rep(fit$mean, NROW(X))
    },
    .package = "SSN2"
  )
  fit <- ssn_decorrelate(
    Summer_mn ~ ., ssn,
    nugget_params = nugget_params("nugget", nugget = 1),
    training = list(test_index = 41:45), evaluate_test = TRUE
  )

  expect_s3_class(fit, "ssn_decorrelate")
  expect_identical(colnames(fit$decorrelate_data$X), c("(Intercept)", "ELEV_DEM"))
  expect_true(all(is.finite(unlist(fit$test[c("bias", "MSPE", "RMSPE")]))))
})

test_that("learner evaluation maps holdouts after missing responses are removed", {
  skip_if_not_installed("ranger")
  ssn <- mf04p
  ssn$obs$decorrelation_response <- ssn$obs$Summer_mn
  ssn$obs$decorrelation_response[[1]] <- NA_real_
  candidates <- SSN2:::get_decorrelate_grid_table(list(
    short = decorrelate_initials(5000),
    long = decorrelate_initials(10000)
  ))
  fit <- ssn_decorrelate(
    decorrelation_response ~ ELEV_DEM, ssn,
    tailup_type = "exponential", additive = "afvArea",
    grid = candidates, training = list(test_index = c(44, 45)),
    num.trees = 10, seed = 2, num.threads = 1
  )

  expect_identical(
    fit$training$splits[[1]]$test_index,
    match(c(44L, 45L), fit$training$response_index)
  )
  expect_identical(fit$newdata, ".missing")
  expect_length(predict(fit), 1)
  expect_true(is.finite(predict(fit)[[1]]))
})

test_that("local candidate evaluation avoids exact and full covariance requests", {
  skip_if_not_installed("ranger")
  ssn_create_bigdist(mf04p, overwrite = TRUE, no_cores = 1, verbose = FALSE)
  candidates <- SSN2:::get_decorrelate_grid_table(list(
    short = decorrelate_initials(5000),
    long = decorrelate_initials(10000)
  ))
  requested <- integer()
  observed_covariance <- SSN2:::get_decorrelate_observed_covariance
  testthat::local_mocked_bindings(
    get_decorrelate_observed_covariance = function(covariance_fit, data) {
      requested <<- c(requested, NROW(data))
      observed_covariance(covariance_fit, data)
    },
    .package = "SSN2"
  )

  fit <- ssn_decorrelate(
    Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    grid = candidates, training = list(test_index = seq_len(5)),
    local = list(method = "covariance", size = 4),
    num.trees = 10, seed = 2, num.threads = 1
  )
  expect_s3_class(fit, "ssn_decorrelate")
  expect_true(length(requested) > 0)
  expect_lte(max(requested), 4)
})

test_that("candidate ranking only requires its selected statistic", {
  evaluation <- data.frame(
    candidate = c("one", "two"), split = "1",
    bias = c(1, 2), MSPE = c(1, 4), RMSPE = c(1, 2),
    cor2 = c(NA_real_, 0.5), failure = NA_character_
  )
  ranked <- SSN2:::get_decorrelate_grid_summary(evaluation, "RMSPE")
  expect_identical(ranked$candidate, c("one", "two"))
  expect_true(is.na(ranked$cor2[[1]]))

  failed <- rbind(
    evaluation,
    data.frame(
      candidate = "one", split = "2", bias = 1, MSPE = 1,
      RMSPE = 1, cor2 = 0.4, failure = "learner failed"
    )
  )
  ranked_after_failure <- SSN2:::get_decorrelate_grid_summary(failed, "RMSPE")
  expect_identical(ranked_after_failure$candidate, "two")

  evaluation$cor2 <- NA_real_
  expect_error(
    SSN2:::get_decorrelate_grid_summary(evaluation, "cor2"),
    "finite cor2"
  )
})

test_that("exact and local decorrelation retain scale equivariance", {
  for (variance in c(1e-10, 1e10)) {
    ssn <- mf04p
    ssn$obs$scaled_response <- sqrt(variance) * ssn$obs$Summer_mn
    params <- nugget_params("nugget", nugget = variance)
    exact <- ssn_decorrelate_data(
      scaled_response ~ ELEV_DEM, ssn,
      nugget_params = params
    )
    local <- ssn_decorrelate_data(
      scaled_response ~ ELEV_DEM, ssn,
      nugget_params = params,
      local = list(method = "covariance", size = NROW(ssn$obs)),
      ordering = "none"
    )
    expect_true(all(is.finite(local$ty)))
    expect_equal(local$ty, exact$ty, tolerance = 1e-9)
    expect_equal(unname(local$tX), unname(exact$tX), tolerance = 1e-9)
  }
})

# decorrelate ties
test_that("decorrelation covariance ties prefer later candidates", {
  expect_identical(get_decorrelate_covariance_neighbors(c(0.8, -0.8, 0, 0), 3), c(2L, 1L, 4L))
  expect_identical(get_decorrelate_covariance_neighbors(rep(0, 5), 2), c(5L, 4L))
  expect_identical(get_decorrelate_covariance_neighbors(c(0.2, -0.9, 0.4), 2), c(2L, 3L))
  expect_identical(get_decorrelate_covariance_neighbors(c(0.5, 0.5), 9), c(2L, 1L))
  expect_identical(get_decorrelate_covariance_neighbors(numeric(), 3), integer())
})

test_that("compact-kernel local decorrelation and recoloring match spmodel", {
  skip_if_not("decorrelate_data" %in% getNamespaceExports("spmodel"))
  network <- mf04p
  network$preds$CapeHorn <- network$preds$CapeHorn[c(1, 7, 15, 23), ]
  observed <- sf::st_drop_geometry(network$obs)
  observed$cx <- sf::st_coordinates(network$obs)[, 1]
  observed$cy <- sf::st_coordinates(network$obs)[, 2]
  newdata <- sf::st_drop_geometry(network$preds$CapeHorn)
  newdata$cx <- sf::st_coordinates(network$preds$CapeHorn)[, 1]
  newdata$cy <- sf::st_coordinates(network$preds$CapeHorn)[, 2]
  for (kernel in c("spherical", "circular", "cubic", "pentaspherical")) {
    a <- ssn_decorrelate_data(Summer_mn ~ ELEV_DEM, network,
      euclid_params = euclid_params(kernel, 0.4, 12000),
      nugget_params = nugget_params("nugget", 0.6),
      ordering = "none", local = list(method = "covariance", size = 12))
    b <- spmodel::decorrelate_data(Summer_mn ~ ELEV_DEM, observed,
      spcov_params = spmodel::spcov_params(kernel, de = 0.4, ie = 0.6, range = 12000),
      xcoord = "cx", ycoord = "cy", ordering = "none", local = list(method = "covariance", size = 12))
    expect_equal(unname(a$tX), unname(b$tX), tolerance = 1e-8)
    expect_equal(unname(a$ty), unname(b$ty), tolerance = 1e-8)
    an <- ssn_decorrelate_newdata(a, "CapeHorn")
    bn <- spmodel::decorrelate_newdata(b, newdata)
    expect_equal(unname(an$tX_newdata), unname(bn$tX_newdata), tolerance = 1e-8)
    beta <- qr.solve(a$tX, a$ty)
    expect_equal(unname(ssn_recorrelate_newdata(an, as.numeric(an$tX_newdata %*% beta))),
      unname(spmodel::recorrelate_newdata(bn, as.numeric(bn$tX_newdata %*% beta))), tolerance = 1e-8)
  }
})

