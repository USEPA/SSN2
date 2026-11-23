skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

order_cache_fixture <- function(method = "split") {
  ssn <- mf04p
  ssn$obs$cache_response <- ssn$obs$Summer_mn
  ssn$obs$cache_response[c(1, 9, 20)] <- NA_real_
  formula <- cache_response ~ ELEV_DEM
  specification <- list(training_index = rev(setdiff(2:16, 9L)),
                        test_index = c(45L, 18L, 30L))
  if (method == "cv") {
    ssn$obs$cache_response[21:45] <- NA_real_
    specification <- list(method = "cv", folds_index = rep(1:3, length.out = 45))
  }
  response_index <- SSN2:::get_decorrelate_response_index(formula, ssn)
  candidate <- function(stream, euclid) list(
    tailup_initial = tailup_initial("exponential", de = stream, range = 10000, known = "given"),
    taildown_initial = taildown_initial("none"),
    euclid_initial = euclid_initial("exponential", de = euclid, range = 5000,
                                    rotate = 0, scale = 1, known = "given"),
    nugget_initial = nugget_initial("nugget", nugget = 0.1, known = "given")
  )
  list(ssn = ssn, formula = formula, specification = specification,
       training = SSN2:::get_decorrelate_training(specification, response_index, 45L),
       candidates = list(stream = candidate(1, 0.001), euclid = candidate(0.001, 1)))
}

order_cache_args <- function(fixture, candidate = fixture$candidates[[1]]) {
  SSN2:::get_decorrelate_data_args(
    fixture$formula, fixture$ssn, candidate$tailup_initial,
    candidate$taildown_initial, candidate$euclid_initial, candidate$nugget_initial,
    "afvArea", candidate$randcov_initial, NULL
  )
}

order_cache_evaluate <- function(fixture, ordering = "random", local = list(size = 3)) {
  SSN2:::get_decorrelate_grid_evaluation(
    fixture$candidates, fixture$training, fixture$formula, fixture$ssn,
    "exponential", "none", "exponential", "nugget", "afvArea", FALSE, NULL,
    NULL, ordering, "ranger", "RMSPE", local, list()
  )
}

order_cache_instrument <- function(.local_envir = parent.frame()) {
  seen <- new.env(parent = emptyenv())
  seen$orders <- seen$preparations <- seen$graphs <- seen$predictions <- list()
  seen$contexts <- seen$fits <- 0L
  seen$events <- character()
  seen$active <- NULL
  prepare <- SSN2:::get_decorrelate_preparation
  order_rows <- SSN2:::get_decorrelate_order
  context <- SSN2:::get_decorrelate_context
  graph <- SSN2:::get_decorrelate_covariance_graph
  neighbors <- SSN2:::get_decorrelate_covariance_neighbors
  data_args <- SSN2:::get_decorrelate_data_args
  newdata_original <- ssn_decorrelate_newdata
  testthat::local_mocked_bindings(
    get_decorrelate_preparation = function(formula, ssn.object, ordering) {
      value <- prepare(formula, ssn.object, ordering)
      seen$preparations[[length(seen$preparations) + 1L]] <- value
      value
    },
    get_decorrelate_order = function(rows, ordering = NULL, data = NULL) {
      seen$events <- c(seen$events, "order")
      value <- order_rows(rows, ordering, data)
      seen$orders[[length(seen$orders) + 1L]] <- list(rows = rows, permutation = value)
      value
    },
    get_decorrelate_data_args = function(...) {
      seen$events <- c(seen$events, "candidate")
      data_args(...)
    },
    get_decorrelate_context = function(...) {
      seen$contexts <- seen$contexts + 1L
      context(...)
    },
    get_decorrelate_covariance_graph = function(covariance_fit, X, y, permutation, local) {
      i <- length(seen$graphs) + 1L
      seen$graphs[[i]] <- list(observed_index = covariance_fit$observed_index,
        missing_index = covariance_fit$missing_index,
        rows = SSN2:::get_decorrelate_rows(covariance_fit$ssn.object$obs),
        missing_rows = SSN2:::get_decorrelate_rows(covariance_fit$ssn.object$preds$.missing),
        permutation = permutation, neighbors = list(), X = X, y = y)
      seen$active <- i
      on.exit(seen$active <- NULL)
      value <- graph(covariance_fit, X, y, permutation, local)
      seen$graphs[[i]]$transformed <- value
      value
    },
    get_decorrelate_covariance_neighbors = function(...) {
      value <- neighbors(...)
      if (!is.null(seen$active)) {
        i <- seen$active
        j <- length(seen$graphs[[i]]$neighbors) + 1L
        seen$graphs[[i]]$neighbors[[j]] <- value
      }
      value
    },
    ssn_decorrelate_newdata = function(object, newdata, ...) {
      seen$predictions[[length(seen$predictions) + 1L]] <-
        SSN2:::get_decorrelate_rows(object$covariance_fit$ssn.object$preds[[newdata]])
      do.call(newdata_original, c(list(object = object, newdata = newdata), list(...)))
    },
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) {
      seen$fits <- seen$fits + 1L
      mean(y)
    },
    predict_decorrelate_algorithm = function(fit, X, algorithm) rep(fit, NROW(X)),
    .package = "SSN2", .env = .local_envir
  )
  seen
}

order_cache_check_reuse <- function(ordering, method) {
  set.seed(2)
  fixture <- order_cache_fixture(method)
  set.seed(2)
  expected <- lapply(fixture$training$splits, function(split) {
    SSN2:::get_decorrelate_split_preparation(fixture$ssn, fixture$formula,
      fixture$training$response_index, split, ordering)$preparation
  })
  expected_rng <- .Random.seed
  seen <- order_cache_instrument()
  set.seed(2)
  result <- order_cache_evaluate(fixture, ordering)
  nsplits <- length(expected)
  expect_true(all(is.na(result$failure)))
  expect_length(seen$orders, nsplits)
  expect_length(seen$preparations, nsplits)
  expect_length(seen$graphs, 2L * nsplits)
  expect_identical(seen$contexts, 2L * nsplits)
  expect_identical(seen$fits, 2L * nsplits)
  expect_identical(seen$events, c(rep("order", nsplits), rep("candidate", 2L)))
  expect_identical(.Random.seed, expected_rng)
  for (i in seq_len(nsplits)) {
    split <- fixture$training$splits[[i]]
    train <- sort(fixture$training$response_index[split$training_index])
    test <- fixture$training$response_index[split$test_index]
    missing <- setdiff(seq_len(NROW(fixture$ssn$obs)), train)
    preparation <- seen$preparations[[i]]
    expect_named(preparation, c("observed_index", "rows", "ordering", "permutation"),
                 ignore.order = TRUE)
    expect_identical(preparation$observed_index, train)
    expect_identical(sort(preparation$permutation), seq_along(train))
    expect_identical(preparation$rows, SSN2:::get_decorrelate_rows(fixture$ssn$obs[train, ]))
    expect_identical(preparation$permutation, expected[[i]]$permutation)
    for (j in c(i, i + nsplits)) {
      graph <- seen$graphs[[j]]
      expect_identical(graph$permutation, preparation$permutation)
      expect_identical(graph$observed_index, train)
      expect_identical(graph$missing_index, missing)
      expect_identical(graph$rows, preparation$rows)
      expect_identical(graph$missing_rows, SSN2:::get_decorrelate_rows(fixture$ssn$obs[missing, ]))
      expect_identical(seen$predictions[[j]], SSN2:::get_decorrelate_rows(fixture$ssn$obs[test, ]))
      expect_equal(graph$y, fixture$ssn$obs$cache_response[train])
      expect_length(graph$neighbors, length(train) - 1L)
    }
  }
  expect_true(any(vapply(seq_len(nsplits), function(i) {
    !identical(seen$graphs[[i]]$neighbors, seen$graphs[[i + nsplits]]$neighbors)
  }, logical(1))))
  expect_false(isTRUE(all.equal(seen$graphs[[1]]$transformed,
                               seen$graphs[[1 + nsplits]]$transformed)))
}

test_that("random ordering is cached per holdout with fresh covariance neighborhoods", {
  order_cache_check_reuse("random", "split")
})

test_that("random ordering is cached separately for each CV fold", {
  order_cache_check_reuse("random", "cv")
})

test_that("actual GRTS ordering is cached per holdout", {
  skip_if_not_installed("spsurvey")
  order_cache_check_reuse("grts", "split")
})

test_that("actual GRTS ordering is cached separately for each CV fold", {
  skip_if_not_installed("spsurvey")
  order_cache_check_reuse("grts", "cv")
})

test_that("fixed ordering matches uncached candidate statistics for split and CV", {
  set.seed(2)
  testthat::local_mocked_bindings(
    fit_decorrelate_algorithm = function(X, y, algorithm, dots) mean(y),
    predict_decorrelate_algorithm = function(fit, X, algorithm) rep(fit, NROW(X))
  )
  for (method in c("split", "cv")) {
    fixture <- order_cache_fixture(method)
    for (ordering in c("pid", "none")) {
      for (local in list(FALSE, list(size = 3))) {
        baseline <- do.call(rbind, lapply(names(fixture$candidates), function(name) {
          args <- order_cache_args(fixture, fixture$candidates[[name]])
          do.call(rbind, lapply(names(fixture$training$splits), function(split_name) {
            statistics <- SSN2:::get_decorrelate_candidate_statistics(
              args, fixture$ssn, fixture$formula, fixture$training$response_index,
              fixture$training$splits[[split_name]], ordering, local, "ranger", list()
            )
            data.frame(candidate = name, split = split_name, statistics, row.names = NULL)
          }))
        }))
        cached <- order_cache_evaluate(fixture, ordering, local)
        expect_true(all(is.na(cached$failure)))
        expect_equal(cached, baseline, tolerance = 1e-10)
        expect_equal(SSN2:::get_decorrelate_grid_summary(cached, "RMSPE"),
                     SSN2:::get_decorrelate_grid_summary(baseline, "RMSPE"), tolerance = 1e-10)
      }
    }
  }
})

test_that("invalid candidates cannot redraw split orders or poison later candidates", {
  set.seed(2)
  fixture <- order_cache_fixture("cv")
  invalid <- fixture$candidates[[1]]
  invalid$tailup_initial$initial[["de"]] <- -1
  fixture$candidates <- c(list(invalid = invalid), fixture$candidates)
  seen <- order_cache_instrument()
  set.seed(2)
  result <- order_cache_evaluate(fixture)
  nsplits <- length(fixture$training$splits)
  expect_length(seen$orders, nsplits)
  expect_length(seen$preparations, nsplits)
  expect_true(all(!is.na(result$failure[result$candidate == "invalid"])))
  expect_true(all(is.na(result$failure[result$candidate != "invalid"])))
  expect_true(all(is.finite(result$RMSPE[result$candidate != "invalid"])))
  expect_length(seen$graphs, 2L * nsplits)
  for (i in seq_len(nsplits)) {
    expect_identical(seen$graphs[[i]]$permutation, seen$preparations[[i]]$permutation)
    expect_identical(seen$graphs[[i + nsplits]]$permutation, seen$preparations[[i]]$permutation)
  }
})

test_that("a split preparation failure is reported for every candidate without redraw", {
  set.seed(2)
  fixture <- order_cache_fixture("cv")
  seen <- order_cache_instrument()
  order_rows <- SSN2:::get_decorrelate_order
  attempts <- 0L
  testthat::local_mocked_bindings(
    get_decorrelate_order = function(rows, ordering = NULL, data = NULL) {
      attempts <<- attempts + 1L
      value <- order_rows(rows, ordering, data)
      if (attempts == 1L) stop("Split preparation failed after ordering", call. = FALSE)
      value
    }
  )
  set.seed(2)
  result <- order_cache_evaluate(fixture)
  nsplits <- length(fixture$training$splits)
  failed <- result$split == names(fixture$training$splits)[[1]]
  expect_identical(attempts, nsplits)
  expect_length(seen$orders, nsplits)
  expect_equal(result$failure[failed], rep("Split preparation failed after ordering", 2))
  expect_true(all(is.na(result$RMSPE[failed])))
  expect_true(all(is.na(result$failure[!failed])))
  expect_true(all(is.finite(result$RMSPE[!failed])))
  expect_length(seen$graphs, 2L * (nsplits - 1L))
  for (i in seq_len(nsplits - 1L)) {
    expect_identical(seen$graphs[[i]]$permutation, seen$orders[[i + 1L]]$permutation)
    expect_identical(seen$graphs[[i + nsplits - 1L]]$permutation,
                     seen$orders[[i + 1L]]$permutation)
  }
})

test_that("split preparation masks every response variable outside training", {
  set.seed(2)
  fixture <- order_cache_fixture()
  fixture$ssn$obs$cache_second <- fixture$ssn$obs$Summer_mn / 10
  fixture$formula <- I(cache_response + cache_second) ~ ELEV_DEM
  before <- fixture$ssn
  split <- fixture$training$splits[[1]]
  train <- sort(fixture$training$response_index[split$training_index])
  preparation <- SSN2:::get_decorrelate_split_preparation(fixture$ssn, fixture$formula,
    fixture$training$response_index, split, "random")
  excluded <- setdiff(fixture$training$response_index, train)
  expect_identical(fixture$ssn, before)
  expect_named(preparation, c("ssn.object", "preparation"))
  for (variable in all.vars(fixture$formula[[2]])) {
    expect_true(all(is.na(preparation$ssn.object$obs[[variable]][excluded])))
    expect_identical(preparation$ssn.object$obs[[variable]][train], before$obs[[variable]][train])
  }
  expect_identical(preparation$preparation$observed_index, train)
  expect_identical(sort(preparation$preparation$permutation), seq_along(train))
})

test_that("preparation rejects stale row maps and ordering", {
  fixture <- order_cache_fixture()
  split <- fixture$training$splits[[1]]
  prepared <- SSN2:::get_decorrelate_split_preparation(fixture$ssn, fixture$formula,
    fixture$training$response_index, split, "pid")
  args <- order_cache_args(fixture)
  args$ssn.object <- prepared$ssn.object
  args$ordering <- "pid"
  args$local <- list(size = 3)
  args$preparation <- prepared$preparation
  expect_s3_class(do.call(SSN2:::ssn_decorrelate_data_internal, args), "ssn_decorrelate_data")
  variants <- list(args$preparation, args$preparation, args$preparation)
  variants[[1]]$rows$pid <- rev(variants[[1]]$rows$pid)
  variants[[2]]$observed_index <- rev(variants[[2]]$observed_index)
  variants[[3]]$ordering <- "none"
  for (preparation in variants) {
    args$preparation <- preparation
    expect_error(do.call(SSN2:::ssn_decorrelate_data_internal, args), "aligned|match|preparation")
  }
  testthat::local_mocked_bindings(
    get_decorrelate_order = function(...) stop("Cached order must be reused")
  )
  args$preparation <- prepared$preparation
  expect_s3_class(do.call(SSN2:::ssn_decorrelate_data_internal, args), "ssn_decorrelate_data")
  for (fun in list(ssn_decorrelate, ssn_decorrelate_data)) {
    expect_false(any(c("preparation", "permutation") %in% names(formals(fun))))
  }
})

test_that("direct data prepares once and the user learner adds one final fit order", {
  set.seed(2)
  fixture <- order_cache_fixture("cv")
  seen <- order_cache_instrument()
  args <- order_cache_args(fixture)
  args$ordering <- "random"
  args$local <- list(size = 3)
  direct <- do.call(ssn_decorrelate_data, args)
  expect_length(seen$orders, 1L)
  expect_length(seen$preparations, 1L)
  expect_identical(seen$preparations[[1]]$observed_index, fixture$training$response_index)
  graph <- seen$graphs[[1]]
  expect_equal(direct$ty, graph$transformed$ty[order(graph$permutation)])
  expect_equal(direct$tX, graph$transformed$tX[order(graph$permutation), , drop = FALSE])
  grid <- SSN2:::get_decorrelate_grid_table(fixture$candidates)
  fit <- ssn_decorrelate(fixture$formula, fixture$ssn, grid = grid, additive = "afvArea",
    training = fixture$specification, ordering = "random", local = list(size = 3))
  nsplits <- length(fixture$training$splits)
  expect_s3_class(fit, "ssn_decorrelate")
  expect_length(seen$orders, nsplits + 2L)
  expect_length(seen$preparations, nsplits + 2L)
  expect_length(seen$graphs, 2L * nsplits + 2L)
  expect_identical(seen$fits, 2L * nsplits + 1L)
  final <- tail(seen$preparations, 1)[[1]]
  expect_identical(final$observed_index, fixture$training$response_index)
  expect_identical(final$permutation, tail(seen$graphs, 1)[[1]]$permutation)
  expect_true(is.finite(fit$test$RMSPE))
})
