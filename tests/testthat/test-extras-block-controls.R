skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("block controls resolve exact names before applying defaults", {
  withr::local_options(warnPartialMatchDollar = TRUE)
  resolve <- SSN2:::get_local_list_prediction_block
  for (method_new in c("subset", "basis")) {
    for (controls in list(list(method_new = method_new),
                         list(method_new = method_new, size_new = 3))) {
      expect_warning(actual <- resolve(controls), NA)
      expect_identical(actual, resolve(c(controls, list(method = "covariance"))))
      expect_identical(actual$method, "covariance")
      expect_identical(actual$size, 4000L)
    }
  }
  expect_identical(resolve(TRUE), resolve(list()))
  expect_identical(resolve(FALSE), resolve(list(method = "all", method_new = "basis", size_new = Inf)))
  expect_null(resolve(list(method = "all", size_new = 3))[["size"]])
  expect_identical(resolve(list(size = 5, size_new = 3))$size, 5L)
  expect_error(resolve(list(method = "subset")), "Invalid local method")
  expect_error(resolve(list(method_new = "covariance")), "local\\$method_new")
  expect_error(resolve(list(method_new = "subset", size_new = 0)), "local\\$size_new")
})

test_that("exported block controls preserve predictions and selected rows", {
  set.seed(2)
  network <- mf04p
  network$obs <- network$obs[seq_len(24), ]
  network$preds$CapeHorn <- network$preds$CapeHorn[seq_len(12), ]
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, network,
    tailup_type = "none", taildown_type = "none",
    euclid_initial = euclid_initial("exponential", 1, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.2, known = "given"))
  withr::local_options(warnPartialMatchDollar = TRUE)
  cross_cov <- as.matrix(covmatrix(fit, "CapeHorn"))
  selected <- NULL
  local_mocked_bindings(coefficients = function(object, ...) {
    frame <- parent.frame()
    if (exists("covariance_rows", frame, inherits = FALSE)) {
      selected <<- list(
        nodes = frame$nodes, grid = frame$covariance_rows,
        obs = if (exists("cov_index", frame, inherits = FALSE)) frame$cov_index else seq_len(object$n)
      )
    }
    stats::coefficients(object, ...)
  }, .package = "SSN2")
  run <- function(controls) {
    expect_warning(result <- predict(fit, "CapeHorn", block = TRUE,
      se.fit = TRUE, local = controls), NA)
    expect_true(all(is.finite(c(result$fit, result$se.fit))))
    list(result = result, selected = selected)
  }

  for (method_new in c("subset", "basis")) {
    controls <- list(method_new = method_new, size = 6, size_new = 3,
                     ordering = "none", chunk_size = 4)
    omitted <- run(controls)
    explicit <- run(c(list(method = "covariance"), controls))
    expect_equal(omitted$result, explicit$result, tolerance = 1e-12)
    expect_identical(omitted$selected, explicit$selected)
    rows <- if (method_new == "subset") 1:3 else 1:12
    c0 <- colMeans(cross_cov[rows, , drop = FALSE])
    expected_obs <- order(c0, decreasing = TRUE)[seq_len(6)]
    expect_identical(omitted$selected$nodes, 1:3)
    expect_identical(omitted$selected$grid, rows)
    expect_identical(omitted$selected$obs, expected_obs)

    all <- run(c(list(method = "all"), controls))
    expect_identical(all$selected$nodes, 1:3)
    expect_identical(all$selected$grid, rows)
    expect_identical(all$selected$obs, seq_len(fit$n))
    full <- run(modifyList(controls, list(size = fit$n)))
    expect_equal(all$result, full$result, tolerance = 1e-10)
  }
  expect_equal(run(FALSE), run(list(method = "all", size_new = Inf)))
  expect_equal(run(TRUE), run(list()))
  expect_equal(run(TRUE)$result, run(FALSE)$result, tolerance = 1e-10)
})

test_that("block controls leave Gaussian and GLM point prediction unchanged", {
  set.seed(2)
  network <- mf04p
  network$obs <- network$obs[seq_len(24), ]
  network$obs$count <- rpois(24, 5)
  network$preds$CapeHorn <- network$preds$CapeHorn[seq_len(3), ]
  common <- list(ssn.object = network, tailup_type = "none", taildown_type = "none",
    euclid_type = "none", nugget_initial = nugget_initial("nugget", 0.2, known = "given"))
  fits <- list(
    do.call(ssn_lm, c(list(formula = Summer_mn ~ ELEV_DEM), common)),
    do.call(ssn_glm, c(list(formula = count ~ ELEV_DEM, family = "poisson"), common))
  )
  withr::local_options(warnPartialMatchDollar = TRUE)
  for (fit in fits) {
    for (method_new in c("subset", "basis")) {
      controls <- list(method_new = method_new, size = 6, size_new = 2)
      expect_warning(omitted <- predict(fit, "CapeHorn", se.fit = TRUE,
        local = controls), NA)
      expect_warning(explicit <- predict(fit, "CapeHorn", se.fit = TRUE,
        local = c(list(method = "covariance"), controls)), NA)
      expect_equal(omitted, explicit)
    }
  }
  expect_error(predict(fits[[2]], "CapeHorn", block = TRUE,
    local = list(method_new = "subset")), "Block prediction is not supported")
})

block_control_fixture <- function() {
  network <- mf04p
  for (nm in c("obs", "CapeHorn")) {
    data <- if (nm == "obs") network$obs else network$preds[[nm]]
    data$group <- factor(rep(letters[1:3], length.out = nrow(data)))
    data$part <- factor(rep(1:2, length.out = nrow(data)))
    if (nm == "obs") network$obs <- data else network$preds[[nm]] <- data
  }
  network
}

block_index_reference <- function(fit, name, keep, nodes, method_new = "subset") {
  grid <- fit$ssn.object$preds[[name]]
  nm <- get_newdata_model_matrix(fit, grid)
  x0 <- colMeans(nm$newdata_model)
  grid_rows <- if (method_new == "subset") nodes else seq_len(nrow(grid))
  c0 <- colMeans(as.matrix(covmatrix(fit, name, cov_type = "pred.obs"))[grid_rows, keep, drop = FALSE])
  s0 <- mean(as.matrix(covmatrix(fit, name, cov_type = "pred.pred"))[grid_rows, nodes, drop = FALSE])
  params <- fit$coefficients$params_object
  nugget <- max(params$nugget[["nugget"]], fit$diagtol,
    1e-4 * sum(params$tailup[["de"]], params$taildown[["de"]], params$euclid[["de"]]))
  if (method_new == "subset") s0 <- s0 + nugget * (1 / nrow(grid) - 1 / length(nodes))
  weights <- solve(as.matrix(covmatrix(fit))[keep, keep, drop = FALSE], c0)
  X <- model.matrix(fit)[keep, , drop = FALSE]
  frame <- model.frame(fit)
  y <- model.response(frame)[keep]
  offset <- model.offset(frame)
  if (!is.null(offset)) y <- y - offset[keep]
  h <- x0 - drop(crossprod(weights, X))
  mean <- sum(x0 * coef(fit)) + sum(weights * (y - X %*% coef(fit)))
  if (!is.null(nm$offset)) mean <- mean + base::mean(nm$offset)
  variance <- s0 - sum(c0 * weights) + drop(h %*% vcov(fit) %*% h)
  list(fit = c("1" = mean), se.fit = c("1" = sqrt(variance)))
}

test_that("block subset reweights only the nugget and supports explicit rows", {
  set.seed(2)
  network <- block_control_fixture()
  tailup <- ssn_lm(Summer_mn ~ ELEV_DEM, network, additive = "afvArea",
    tailup_initial = tailup_initial("exponential", .8, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", .1, known = "given"))
  mixed <- update(tailup,
    taildown_initial = taildown_initial("spherical", .4, 10000, known = "given"),
    euclid_initial = euclid_initial("exponential", .2, 5000, known = "given"))
  random <- update(mixed, random = ~group, randcov_initial = randcov_initial(group = .3, known = "given"))
  partition <- update(mixed, partition_factor = ~part)
  nodes <- c(1L, 6L, 10L)
  for (fit in list(tailup, mixed, random, partition)) {
    controls <- list(method = "all", method_new = "subset", ordering = "none", size_new = 3)
    actual <- predict(fit, "CapeHorn", block = TRUE, se.fit = TRUE, local = controls)
    ref <- block_index_reference(fit, "CapeHorn", seq_len(fit$n), 1:3)
    expect_equal(actual$se.fit^2, ref$se.fit^2, tolerance = 1e-10)
    for (method_new in c("subset", "basis")) {
      keep <- c(1L, 4L, 7L, 10L)
      ref <- block_index_reference(fit, "CapeHorn", keep, nodes, method_new)
      for (method in c("all", "covariance")) {
        actual <- predict(fit, "CapeHorn", block = TRUE, se.fit = TRUE, local = list(
          method = method, size = 1, index = rev(keep), index_new = rev(nodes),
          size_new = Inf, method_new = method_new, ordering = "none"))
        expect_equal(actual, ref, tolerance = 1e-12)
      }
    }
  }
  for (fit in list(mixed, random, partition)) {
    context <- get_conditional_vecchia_covariance(fit)
    for (keep in list(2L, c(10L, 1L, 7L, 4L))) {
      expect_equal(unname(conditional_vecchia_covariance(context, keep)),
        unname(as.matrix(covmatrix(fit))[keep, keep, drop = FALSE]), tolerance = 1e-12)
    }
  }
})

test_that("observation submatrices preserve the nugget floor and fitted row order", {
  set.seed(2)
  network <- block_control_fixture()
  missing <- c(2L, 12L, 22L)
  network$obs$Summer_mn[missing] <- NA
  taildown <- ssn_lm(Summer_mn ~ offset(ELEV_DEM / 100) + SLOPE, network,
    taildown_initial = taildown_initial("exponential", .8, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", .001, known = "given"))
  taildown$diagtol <- .05
  iid <- update(taildown, taildown_initial = taildown_initial("none"))
  expect_gt(taildown$diagtol, taildown$coefficients$params_object$nugget[["nugget"]])
  expect_equal(unname(diag(as.matrix(covmatrix(taildown)))), rep(.85, taildown$n))
  for (fit in list(taildown, iid)) {
    expect_identical(fit$ssn.object$obs$pid, network$obs$pid[-missing])
    expect_identical(fit$ssn.object$preds$.missing$pid, network$obs$pid[missing])
    expect_identical(rownames(model.frame(fit)), rownames(fit$ssn.object$obs))
    for (keep in list(2L, c(10L, 1L, 7L, 4L))) {
      actual <- conditional_vecchia_covariance(get_conditional_vecchia_covariance(fit), keep)
      expect_equal(unname(actual), unname(as.matrix(covmatrix(fit))[keep, keep, drop = FALSE]), tolerance = 1e-12)
    }
    keep <- c(2L, 7L, 11L)
    ref <- block_index_reference(fit, ".missing", keep, 2L)
    controls <- list(method = "all", index = rev(keep), method_new = "subset", index_new = 2L)
    actual <- predict(fit, ".missing", block = TRUE, se.fit = TRUE, local = controls)
    expect_equal(actual, ref, tolerance = 1e-12)
    controls$index <- seq_len(fit$n) %in% keep
    controls$index_new <- c(FALSE, TRUE, FALSE)
    expect_identical(predict(fit, ".missing", block = TRUE, se.fit = TRUE, local = controls), actual)
  }
})

test_that("explicit block rows replace selection but not mean-only predictions", {
  set.seed(2)
  network <- block_control_fixture()
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, network,
    euclid_initial = euclid_initial("exponential", 1, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", .01, known = "given"))
  run <- function(local, ...) predict(fit, "CapeHorn", block = TRUE, se.fit = TRUE, local = local, ...)
  grid <- fit$ssn.object$preds$CapeHorn
  c0 <- colMeans(as.matrix(covmatrix(fit, "CapeHorn")))
  keep <- tail(order(c0), 6)
  expect_equal(run(list(size = 6)), run(list(index = keep)), tolerance = 1e-12)
  for (method_new in c("subset", "basis")) {
    nodes <- get_block_nodes(grid, 4, "pid")
    expect_equal(run(list(method_new = method_new, size_new = 4)),
      run(list(method_new = method_new, index_new = nodes)), tolerance = 1e-12)
  }
  expect_identical(run(list(method = "all", index_new = rev(seq_len(nrow(grid))), size_new = 1)), run(FALSE))
  for (name in c("index", "index_new")) {
    for (invalid in list(FALSE, c(1, NA), 0, 10000, 1.5, Inf, c(1, 1), integer(), "1", factor(1), matrix(1))) {
      expect_error(run(setNames(list(invalid), name)), paste0("local\\$", name))
    }
  }
  expect_error(run(list(index_new = 1, size_new = -1)), "local\\$size_new")
  expect_error(run(list(index_new = 1, ordering = "bad")), "ordering")
  expect_identical(run(list(index = 0, index_new = 0), type = "terms"), run(FALSE, type = "terms"))
  expect_identical(run(list(index = 0, index_new = 0), interval = "confidence"), run(FALSE, interval = "confidence"))
  expect_error(predict(fit, c("CapeHorn", "pred1km"), block = TRUE, local = list(index_new = 1)),
    "local\\$index_new can only be used with a single prediction set")
  expect_length(predict(fit, c("CapeHorn", "pred1km"), block = TRUE, local = list(index = keep)), 2)
})

test_that("dense subset prediction retains the nugget variance floor", {
  set.seed(2)
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    euclid_initial = euclid_initial("exponential", 1, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", .01, known = "given"))
  G <- nrow(fit$ssn.object$preds$CapeHorn)
  controls <- list(method = "all", method_new = "subset", size_new = 3, ordering = "none")
  run <- function(local) predict(fit, "CapeHorn", block = TRUE, se.fit = TRUE, local = local)
  variance <- unname(run(controls)$se.fit^2)
  expect_equal(variance, .2278926334921867^2, tolerance = 1e-10)
  expect_gte(variance, .01 / G * (1 - 1e-8))
  controls$size_new <- G
  expect_identical(run(controls), run(FALSE))
})
