skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_context_fixture <- function(stream = TRUE, structured = FALSE, missing = FALSE) {
  network <- mf04p
  take <- round(seq(1, NROW(network$obs), length.out = 12))
  network$obs <- network$obs[take, ]
  network$preds <- list(pred1km = network$preds$pred1km[c(5, 1, 8, 3), ])
  network$obs$x <- seq(-1, 1, length.out = 12)
  network$obs$y <- seq_len(12)
  network$preds$pred1km$x <- c(-.5, .2, .8, 0)
  if (stream) {
    root <- tempfile("conditional-context-")
    dir.create(root)
    stopifnot(file.copy(network$path, root, recursive = TRUE, copy.mode = FALSE))
    network$path <- file.path(root, basename(network$path))
    ssn_create_distmat(network, predpts = "pred1km", among_predpts = TRUE,
                       overwrite = TRUE)
  }
  if (missing) network$obs$y[c(2, 7)] <- NA
  fit <- ssn_lm(y ~ x, network, additive = "afvArea", ddf = "asymptotic",
    tailup_initial = if (stream) tailup_initial("exponential", .4, 10000, known = "given") else tailup_initial("none"),
    taildown_initial = if (stream) taildown_initial("spherical", .3, 10000, known = "given") else taildown_initial("none"),
    euclid_initial = euclid_initial("exponential", .2, 12000, rotate = .3, scale = .6, known = "given"),
    nugget_initial = nugget_initial("nugget", .1, known = "given"), anisotropy = TRUE)
  if (structured) {
    obs <- fit$ssn.object$obs
    pred <- fit$ssn.object$preds$pred1km
    obs$g <- factor(rep(c("a", "b"), length.out = NROW(obs)))
    obs$h <- factor(rep(c("c", "d", "c"), length.out = NROW(obs)))
    obs$part <- factor(rep(c("one", "one", "two"), length.out = NROW(obs)))
    pred$g <- factor(c("a", "b", "new", "a"))
    pred$h <- factor(c("c", "d", "c", "d"))
    pred$part <- factor(c("one", "two", "new", "one"))
    fit$ssn.object$obs <- obs
    fit$ssn.object$preds$pred1km <- pred
    fit$random <- ~ g + (0 + x | g:h)
    fit$coefficients$params_object$randcov <- c("1 | g" = .15, "x | g:h" = .1)
    fit$random_xlev <- get_randcov_xlev(names(fit$coefficients$params_object$randcov), obs)
    fit$partition_factor <- ~part
    fit$partition_xlev <- get_partition_xlev(fit$partition_factor, obs)
  }
  fit
}

conditional_context_dense <- function(fit, name) {
  obs <- covmatrix(fit)
  cross <- covmatrix(fit, name, cov_type = "obs.pred")
  pred <- covmatrix(fit, name, cov_type = "pred.pred")
  rbind(cbind(obs, cross), cbind(t(cross), pred))
}

test_that("selected covariance blocks match dense stream and Euclidean covariance", {
  for (stream in c(FALSE, TRUE)) for (structured in c(FALSE, TRUE)) {
    fit <- conditional_context_fixture(stream, structured)
    before <- fit
    reference <- conditional_context_dense(fit, "pred1km")
    context <- get_conditional_vecchia_covariance(fit, "pred1km")
    rows <- c(15, 9, 2, 13)
    cols <- c(6, 14, 1)
    for (r in list(rows, 1L, 14L)) for (c in list(cols, 3L, 15L)) {
      expect_equal(unname(conditional_vecchia_covariance(context, r, c)),
                   unname(reference[r, c, drop = FALSE]), tolerance = 1e-10)
    }
    expect_equal(unname(conditional_vecchia_covariance(context, rows)),
                 unname(reference[rows, rows]), tolerance = 1e-10)
    expect_equal(dim(conditional_vecchia_covariance(context, integer(), cols)), c(0, 3))
    spatial <- get_conditional_vecchia_covariance(fit, "pred1km", spatial = TRUE)
    expected <- conditional_context_dense(conditional_spatial_object(fit), "pred1km")
    expect_equal(unname(conditional_vecchia_covariance(spatial, rows)),
                 unname(expected[rows, rows]), tolerance = 1e-10)
    expect_identical(fit, before)
    expect_error(conditional_vecchia_covariance(context, c(1, NA)), "Covariance indices")
  }
})

test_that("observed-only and missing-response supports preserve identity and nuggets", {
  fit <- conditional_context_fixture(missing = TRUE)
  context <- get_conditional_vecchia_covariance(fit, ".missing")
  reference <- conditional_context_dense(fit, ".missing")
  rows <- c(fit$n + 1L, 4L, fit$n + 2L)
  expect_equal(unname(conditional_vecchia_covariance(context, rows)),
               unname(reference[rows, rows]), tolerance = 1e-10)
  fit <- conditional_context_fixture(stream = FALSE)
  fit$ssn.object$path <- tempfile("no-distance-files-")
  sf::st_geometry(fit$ssn.object$obs)[2] <- sf::st_geometry(fit$ssn.object$obs)[1]
  fit$ssn.object$preds <- NULL
  context <- get_conditional_vecchia_covariance(fit)
  block <- conditional_vecchia_covariance(context, c(2, 1), c(1, 2))
  expect_equal(unname(block), unname(covmatrix(fit)[c(2, 1), c(1, 2)]))
  expect_equal(block[1, 2] - block[1, 1], .1)
  expect_identical(context$backends, list(observed = "none"))
})

test_that("random designs preserve tuple identities and base conditioning membership", {
  observed <- data.frame(g = factor(c("a:b", "a", "a:b", "a")),
    h = factor(c("c", "b:c", "c", "b:c")), p = factor(c("x", "x", "y", "y")),
    x = c(0, 2, -1, 3))
  newdata <- observed[c(2, 3), ]
  newdata$g <- c("new", "a:b")
  params <- c("1 | g:h" = .4, "x | g:h" = .2, "1 | g" = 0)
  full <- get_randcov_design(observed, params, ~p, newdata)
  data <- rbind(observed, newdata)
  group <- outer(as.character(data$g), as.character(data$g), "==") &
    outer(as.character(data$h), as.character(data$h), "==")
  partition <- outer(as.character(data$p), as.character(data$p), "==")
  reference <- (.4 + .2 * outer(data$x, data$x)) * group * partition
  Z <- rbind(full$observed, full$newdata)
  expect_equal(as.matrix(Z %*% Matrix::Diagonal(x = full$variance) %*% t(Z)),
               reference, ignore_attr = TRUE)
  expect_identical(anyDuplicated(rownames(full$info)), 0L)
  base <- get_randcov_design(observed, params, ~p, newdata, index = 1:2)
  expect_true(any(base$info$observed & !base$info$conditioned))
  expect_true(any(!base$info$observed & !base$info$conditioned))
  expect_equal(NROW(base$observed), 2)
  expect_true(all(c("level_tuple", "partition_tuple") %in% names(base$info)))
  expect_equal(NCOL(get_randcov_design(observed, NULL)$observed), 0)
  observed$g[1] <- NA
  expect_error(get_randcov_design(observed, params), "must not be missing")
})

test_that("score vectors are reused and respect true partitions", {
  context <- get_conditional_vecchia_covariance(conditional_context_fixture(structured = TRUE))
  scores <- conditional_vecchia_scores(context, 1, "covariance")
  expect_true(is.infinite(scores[1]))
  expect_true(all(is.infinite(scores[context$partition != context$partition[1]])))
  local_mocked_bindings(conditional_vecchia_scores = function(...) stop("Scores recomputed"))
  pool <- seq_len(context$n)
  expect_equal(conditional_vecchia_neighbors(context, 1, pool, 2, "covariance", scores),
               pool[order(scores, pool)][1:2])
  expect_identical(conditional_vecchia_neighbors(context, 1, pool, 0, "covariance", scores), integer())
  expect_identical(conditional_vecchia_neighbors(context, 1, pool, 0, "all"), pool[-1])
})

test_that("selected bmat reads and reusable contexts match dense blocks", {
  fit <- conditional_context_fixture()
  ssn_create_bigdist(fit$ssn.object, predpts = "pred1km", among_predpts = TRUE,
    overwrite = TRUE, no_cores = 1, verbose = FALSE)
  context <- get_conditional_vecchia_covariance(fit, "pred1km")
  expect_true(all(unlist(context$backends) == "bigdata"))
  reference <- conditional_context_dense(fit, "pred1km")
  original_pred <- get_dist_pred_object
  original_open <- fm.open
  selected_sizes <- list()
  opened <- 0L
  local_mocked_bindings(
    get_conditional_random_design = function(...) stop("Design rebuilt"),
    get_block_pred_backend = function(...) stop("Backend reselected"),
    get_decorrelate_observed_backend = function(...) stop("Backend reselected"),
    get_dist_pred_object = function(object, newdata_name, initial_object, backend = NULL) {
      selected_sizes[[length(selected_sizes) + 1L]] <<- c(NROW(object$ssn.object$obs),
        NROW(object$ssn.object$preds[[newdata_name]]))
      original_pred(object, newdata_name, initial_object, backend)
    },
    fm.open = function(...) { opened <<- opened + 1L; original_open(...) }
  )
  rows <- c(14, 2, 13)
  cols <- c(9, 15, 1)
  expect_equal(unname(conditional_vecchia_covariance(context, rows, cols)),
               unname(reference[rows, cols]), tolerance = 1e-10)
  expect_true(opened > 0)
  expect_true(all(vapply(selected_sizes, function(x) all(x <= 2), logical(1))))
})

test_that("observed scoring reads no prediction distances and limits candidate columns", {
  fit <- conditional_context_fixture()
  ssn_create_bigdist(fit$ssn.object, predpts = "pred1km", among_predpts = TRUE,
    overwrite = TRUE, no_cores = 1, verbose = FALSE)
  context <- get_conditional_vecchia_covariance(fit, "pred1km")
  original <- get_dist_object_bigdata_cross
  square_reads <- character()
  local_mocked_bindings(
    get_dist_pred_object = function(...) stop("Unexpected prediction-distance read"),
    get_dist_object_bigdata_cross = function(..., ext, backend) {
      square_reads <<- c(square_reads, ext)
      original(..., ext = ext, backend = backend)
    })
  scores <- conditional_vecchia_scores(context, 1, "covariance", seq_len(context$n))
  expect_equal(length(scores), context$n)
  expect_true(length(square_reads) > 0)
  expect_true(all(square_reads == "obs"))
  expect_identical(names(scores), as.character(seq_len(context$n)))
  local_mocked_bindings(conditional_vecchia_scores = function(...) stop("Scores recomputed"))
  expect_equal(length(conditional_vecchia_neighbors(context, 1, 1:12, 3, "covariance", scores)), 3)
  expect_equal(length(conditional_vecchia_neighbors(context, 1, 1:5, 2, "covariance", scores)), 2)
})

test_that("prediction scoring reads only observed and earlier prediction candidates", {
  context <- get_conditional_vecchia_covariance(conditional_context_fixture(), "pred1km")
  columns_read <- NULL
  original <- conditional_vecchia_covariance
  local_mocked_bindings(conditional_vecchia_covariance = function(context, rows, columns) {
    columns_read <<- columns
    original(context, rows, columns)
  })
  pool <- seq_len(context$n + 1L)
  conditional_vecchia_neighbors(context, context$n + 2L, pool, 3, "covariance")
  expect_identical(columns_read, pool)
})

test_that("eligible-only scores preserve neighbors and exclude partitions before reads", {
  context <- get_conditional_vecchia_covariance(conditional_context_fixture(structured = TRUE), "pred1km")
  original <- conditional_vecchia_covariance
  columns_read <- integer()
  local_mocked_bindings(conditional_vecchia_covariance = function(context, rows, columns) {
    columns_read <<- columns
    original(context, rows, columns)
  })
  all_rows <- seq_len(context$n + NROW(context$newdata))
  for (i in c(1L, context$n + 4L)) {
    columns <- if (i <= context$n) seq_len(context$n) else seq_len(i - 1L)
    full <- setNames(-abs(as.numeric(original(context, i, all_rows))), as.character(all_rows))
    full[context$partition != context$partition[i] | all_rows == i] <- Inf
    scores <- conditional_vecchia_scores(context, i, "covariance", columns)
    expect_identical(columns_read,
      columns[columns != i & context$partition[columns] == context$partition[i]])
    expect_equal(scores, full[as.character(columns)], tolerance = 1e-12)
    for (pool in list(columns, rev(columns[seq_len(min(5, length(columns)))]))) {
      expect_identical(conditional_vecchia_neighbors(context, i, pool, 2, "covariance", scores),
                       conditional_vecchia_neighbors(context, i, pool, 2, "covariance", full))
    }
  }
})

test_that("new conditional resolver uses exact names and caller-specific defaults", {
  set.seed(2)
  fit <- conditional_context_fixture()
  pred <- fit$ssn.object$preds$pred1km
  settings <- get_local_list_conditional(TRUE, fit, pred)
  expect_identical(settings$approximation, "vecchia")
  expect_identical(settings$size, 30L)
  expect_identical(settings$ordering, "pid")
  expect_identical(settings$chunk_size, 1000L)
  expect_equal(settings$order[settings$inv_order], seq_len(NROW(pred)))
  glm <- fit
  class(glm) <- c("ssn_glm", "ssn_lm")
  expect_identical(get_local_list_conditional(TRUE, glm, pred)$size, 60L)
  expect_identical(get_local_list_conditional(TRUE, glm)$ordering_observed, "pid")
  expect_true(get_local_list_conditional(FALSE, fit, pred)$exact)
  expect_true(get_local_list_conditional(list(method = "all"), fit, pred)$exact)
  expect_true(get_local_list_conditional(list(approximation = "low-rank"), fit, pred)$exact)
  expect_true(get_local_list_conditional(list(approximation = "low-rank",
    method = "covariance", size = 2, ordering = "none", size_response = 2), fit, pred)$exact)
  expect_error(get_local_list_conditional(list(method_new = "subset"), fit, pred), "Base and block settings")
  expect_error(get_local_list_conditional(list(size_response = 2), fit, pred), "size_response is not supported")
  expect_error(get_local_list_conditional(list(size = 1.5), fit, pred), "positive integer")
  expect_error(get_local_list_conditional(list(chunk_size = Inf), fit, pred), "positive integer")
  expect_warning(ignored <- get_local_list_conditional(list(parallel = TRUE), fit, pred), "Ignoring", fixed = TRUE)
  expect_false(ignored$parallel)
  expect_error(get_local_list_conditional(list(size_random = 2), fit, pred), "only used for Vecchia")
  fit$n <- 5001L
  expect_message(auto <- get_local_list_conditional(NULL, fit, pred), "exceeds 5,000")
  expect_identical(auto$approximation, "vecchia")
})

test_that("effect blocks precede prediction randomness and validate support", {
  set.seed(2)
  fit <- conditional_context_fixture(structured = TRUE)
  pred <- fit$ssn.object$preds$pred1km
  local <- list(groups_random = 3, method_random = "random", ordering = "random")
  set.seed(2)
  a <- get_local_list_conditional(local, fit, pred)
  set.seed(2)
  b <- get_local_list_conditional(local, fit, pred[1, ])
  expect_identical(a$random_blocks, b$random_blocks)
  expect_equal(sort(unlist(a$random_blocks)), seq_len(fit$n))
  expect_equal(a$size_random, 4)
  expect_equal(get_local_list_conditional(list(size_random = 6, groups_random = 3), fit, pred)$groups_random, 2)
  expect_error(get_local_list_conditional(list(method = "all", size_random = 3), fit, pred), "only used for Vecchia")
  sf::st_geometry(fit$ssn.object$obs) <- rep(sf::st_geometry(fit$ssn.object$obs)[1], fit$n)
  expect_error(get_local_list_conditional(list(groups_random = 2), fit, pred), "Too few distinct")
})
