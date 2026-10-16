skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_iid_fixture <- function(family = NULL, missing = FALSE, local = FALSE,
                                    variance = .35, estimated = FALSE, random = FALSE) {
  network <- mf04p
  network$obs <- network$obs[round(seq(1, NROW(network$obs), length.out = 16)), ]
  network$preds <- list(iid = network$preds$pred1km[c(5, 1, 8, 3), ])
  network$obs$x <- seq(-1, 1, length.out = 16)
  network$obs$off <- seq(-.2, .2, length.out = 16)
  network$obs$y <- c(1, 3, 2, 4, 2, 5, 4, 2, 5, 3, 6, 4, 2, 5, 1, 4)
  network$obs$part <- factor(rep(c("one", "two"), 8))
  network$preds$iid$x <- c(-.5, .2, .8, 0)
  network$preds$iid$off <- c(.1, -.3, .2, .4)
  network$preds$iid$part <- factor(c("one", "two", "one", "two"))
  if (identical(family, "beta")) network$obs$y <- network$obs$y / 8
  if (missing) network$obs$y[c(2, 7)] <- NA
  args <- list(formula = y ~ x + offset(off), ssn.object = network,
    tailup_type = "none", taildown_type = "none", euclid_type = "none",
    partition_factor = ~part, local = local,
    nugget_initial = nugget_initial("nugget", variance,
      known = if (estimated) NULL else "given"))
  if (is.null(family)) args$ddf <- "asymptotic" else args$family <- family
  if (random) {
    args$random <- ~part
    args$randcov_initial <- randcov_initial(part = .2, known = "given")
  }
  if (identical(family, "binomial")) args$formula <- cbind(y, 8 - y) ~ x + offset(off)
  if (!is.null(family) && !family %in% c("poisson", "binomial")) {
    args$dispersion_initial <- do.call(dispersion_initial, list(family, dispersion = 3, known = "given"))
  }
  do.call(if (is.null(family)) ssn_lm else ssn_glm, args)
}

test_that("independent Gaussian draws retain the dense exact random fixture", {
  set.seed(2)
  for (missing in c(FALSE, TRUE)) for (estimated in c(FALSE, TRUE)) {
    fit <- conditional_iid_fixture(missing = missing, estimated = estimated)
    name <- if (missing) ".missing" else "iid"
    for (samples in c(1, 1003)) {
      set.seed(2)
      dense <- draw_conditional_exact(fit, name, samples, c("newdata", "beta", "object"))
      set.seed(2)
      direct <- conditional(fit, name, samples = samples, output = "all", local = FALSE)
      expect_identical(direct, dense[names(direct)])
    }
  }
})

test_that("independent models override local controls before validating them", {
  set.seed(2)
  controls <- list(TRUE, FALSE, NA, "invalid", list(size = -1),
    list(approximation = "low-rank", ordering = "maxmin", kmeans_new = TRUE),
    list(approximation = "vecchia", size = 1, parallel = TRUE))
  for (family in list(NULL, "poisson")) {
    fit <- conditional_iid_fixture(family)
    set.seed(2)
    expected <- conditional(fit, "iid", samples = 5, local = FALSE, output = "all")
    for (control in controls) {
      set.seed(2)
      expect_identical(conditional(fit, "iid", samples = 5, local = control, output = "all"), expected)
    }
    expect_error(conditional(fit, output = "random"), "requires a model with random effects")
    expect_error(conditional(fit, "absent"), "not a valid")
    expect_error(conditional(fit, "iid", samples = 0), "positive integer")
  }
  fit <- conditional_iid_fixture()
  expected <- loocv(fit, cv_predict = TRUE, se.fit = TRUE, interval = "prediction")
  for (control in controls) {
    expect_identical(loocv(fit, cv_predict = TRUE, se.fit = TRUE,
      interval = "prediction", local = control), expected)
  }
})

test_that("independent eligibility preserves variance floors and structural exclusions", {
  set.seed(2)
  for (family in list(NULL, "poisson")) {
    fit <- conditional_iid_fixture(family)
    expect_equal(get_conditional_iid_variance(fit), .35)
    expect_null(get_conditional_iid_variance(fit, simulate_covparams = TRUE))
    for (type in c("tailup", "taildown", "euclid")) {
      spatial <- fit
      class(spatial$coefficients$params_object[[type]]) <- paste0(type, "_exponential")
      expect_false(is_conditional_iid(spatial))
      expect_null(get_conditional_iid_variance(spatial))
    }
    random <- fit
    random$random <- ~part
    expect_false(is_conditional_iid(random))
    expect_null(get_conditional_iid_variance(random))
    # Zero-variance Gaussian fits cannot reach simulation through public fitting.
    fit$coefficients$params_object$nugget[["nugget"]] <- 0
    if (is.null(family)) {
      expect_null(get_conditional_iid_variance(fit))
      expect_error(conditional(fit, "iid", samples = 1, local = list(size = -1)), "not positive")
    } else {
      expect_equal(get_conditional_iid_variance(fit), 1e-4)
      class(fit$coefficients$params_object$nugget) <- "nugget_none"
      expect_equal(get_conditional_iid_variance(fit), 1e-4)
    }
    class(fit$coefficients$params_object$nugget) <- "nugget_nugget"
    fit$coefficients$params_object$nugget[["nugget"]] <- NaN
    expect_error(get_conditional_iid_variance(fit), "Invalid independent-error variance")
  }
})

test_that("independent GLM outputs agree with dense formulas across families", {
  set.seed(2)
  for (family in c("poisson", "binomial", "Gamma", "inverse.gaussian", "nbinomial", "beta")) {
    fit <- conditional_iid_fixture(family)
    outputs <- c("newdata", "beta", "object", "latent")
    set.seed(2)
    dense <- draw_conditional_exact(fit, "iid", 1003, outputs)
    set.seed(2)
    direct <- conditional(fit, "iid", samples = 1003, output = outputs)
    expect_identical(direct$beta, dense$beta)
    expect_equal(direct, dense[names(direct)], tolerance = 1e-12)
    set.seed(2)
    link <- conditional(fit, "iid", samples = 5)
    mu <- invlink(link, family, size = 1)
    set.seed(2)
    response <- conditional(fit, "iid", samples = 5, type = "response", newdata_size = 8)
    expect_equal(response, if (family == "binomial") 8 * mu else mu)
    set.seed(2)
    conditional(fit, "iid", samples = 5)
    reference <- draw_glm_response(family, mu, as.vector(coef(fit, "dispersion")), 8)
    set.seed(2)
    expect_equal(conditional(fit, "iid", samples = 5, type = "new", newdata_size = 8), reference)
    for (output in c("beta", "object", "latent")) {
      set.seed(2)
      expect_equal(dim(conditional(fit, output = output, samples = 1)),
        c(if (output == "beta") length(coef(fit)) else fit$n, 1))
    }
  }
})

test_that("independent draws recover shared coefficient and latent uncertainty", {
  set.seed(2)
  fit <- conditional_iid_fixture("poisson")
  X <- model.matrix(fit)
  X0 <- get_newdata_model_matrix(fit, fit$ssn.object$preds$iid)$newdata_model
  B <- vcov(fit)
  v <- get_conditional_iid_variance(fit)
  eta <- fitted(fit, type = "link")
  q <- 1 / (1 / v + exp(eta))
  A <- (q / v) * X
  mapping <- rbind(diag(NCOL(X)), A, X0)
  target <- mapping %*% B %*% t(mapping) + diag(c(rep(0, NCOL(X)), q, rep(v, NROW(X0))))
  mean <- c(coef(fit), eta, X0 %*% coef(fit) + fit$ssn.object$preds$iid$off)
  draws <- conditional(fit, "iid", output = c("beta", "latent", "newdata"), samples = 12000)
  value <- do.call(rbind, draws)
  expect_lt(max(abs(rowMeans(value) - mean) / sqrt(diag(target) / NCOL(value))), 6)
  mcse <- sqrt((outer(diag(target), diag(target)) + target^2) / (NCOL(value) - 1))
  expect_lt(max(abs(cov(t(value)) - target) / mcse), 6)
})

test_that("independent simulation avoids covariance and distance preparation", {
  set.seed(2)
  fits <- list(conditional_iid_fixture(), conditional_iid_fixture("poisson"))
  blocked <- function(...) stop("unexpected dense or distance work")
  local_mocked_bindings(covmatrix = blocked, get_dist_object = blocked,
    get_dist_pred_object = blocked, get_dist_predbk_object = blocked,
    get_conditional_glm_joint = blocked, prepare_conditional_exact_prediction = blocked,
    partition_matrix = blocked)
  for (fit in fits) {
    for (n in c(4999, 5000, 5001)) {
      large <- fit
      large$n <- n
      large$ssn.object$preds$iid <- fit$ssn.object$preds$iid[rep(1:4, length.out = n), ]
      large$ssn.object$preds$iid$pid <- as.character(seq_len(n))
      expect_message(draws <- conditional(large, "iid", samples = 2), NA)
      expect_equal(dim(draws), c(n, 2))
    }
    expect_equal(dim(conditional(fit, output = "beta", samples = 2)), c(2, 2))
    if (inherits(fit, "ssn_glm")) {
      expect_equal(dim(conditional(fit, output = "latent", samples = 2)), c(16, 2))
    }
    fit$ssn.object$preds$iid$part <- NULL
    expect_error(conditional(fit, "iid", samples = 1), "part")
  }
})

test_that("independent latent draws reject invalid conditional precision", {
  set.seed(2)
  fit <- conditional_iid_fixture("beta")
  v <- get_conditional_iid_variance(fit)
  X0 <- model.matrix(fit)[FALSE, , drop = FALSE]
  curvature <- .5 / v
  local_mocked_bindings(get_D = function(...) Matrix::Diagonal(fit$n, curvature))
  expect_true(all(is.finite(get_conditional_iid_glm(fit, X0, v, 3, latent = TRUE)$latent)))
  for (curvature in c(1 / v, 2 / v, NA, Inf)) {
    expect_error(get_conditional_iid_glm(fit, X0, v, 3, latent = TRUE), "latent precision is not positive definite")
  }
})

test_that("independent covariance-parameter simulation is not cancelled by local", {
  set.seed(2)
  fit <- conditional_iid_fixture()
  local_mocked_bindings(get_satterthwaite_cached = function(...) list(vcov_theta = diag(1)),
    draw_conditional_covparams = function(...) list(cov = matrix(7, 1, 2)),
    draw_conditional_iid = function(...) stop("unexpected shortcut"))
  expect_message(value <- conditional(fit, output = "cov", samples = 2,
    simulate_covparams = TRUE, local = list(size = -1)), NA)
  expect_identical(value, matrix(7, 1, 2))
})

test_that("low-rank grouping diagnoses repeated locations without changing valid groups", {
  set.seed(2)
  fit <- conditional_iid_fixture()
  pred <- fit$ssn.object$preds$iid
  index <- c(4, 2, 1, 3)
  set.seed(2)
  reference <- split(index, stats::kmeans(sf::st_coordinates(pred[index, ]), centers = 2, iter.max = 30)$cluster)
  set.seed(2)
  expect_identical(get_local_lowrank_blocks(index, pred, 2, TRUE), reference)
  sf::st_geometry(pred) <- sf::st_geometry(pred)[rep(1, NROW(pred))]
  settings <- list(approximation = "low-rank", size_base = 1, size_new = 1,
                   reorder_base = "none", reorder_new = "none", kmeans_new = TRUE)
  expect_error(get_local_list_conditional(settings, fit, pred), "Too few distinct prediction locations")
  expect_error(get_ssn_simulate_local_lowrank(settings, NROW(pred), pred), "Too few distinct prediction locations")
  settings$kmeans_new <- FALSE
  expect_identical(unname(unlist(get_local_list_conditional(settings, fit, pred)$index$new)), 1:4)
  expect_identical(unname(unlist(get_ssn_simulate_local_lowrank(settings, NROW(pred), pred)$index$new)), 2:4)
})

test_that("local fitted independent models retain global coefficient coupling", {
  set.seed(2)
  fit <- conditional_iid_fixture("poisson", local = list(index = rep(1:2, each = 8)))
  set.seed(2)
  reference <- draw_conditional_exact(fit, "iid", 1003, c("beta", "latent", "newdata"))
  set.seed(2)
  value <- conditional(fit, "iid", output = c("beta", "latent", "newdata"), samples = 1003)
  expect_identical(value$beta, reference$beta)
  expect_equal(value, reference[names(value)], tolerance = 1e-12)
})

test_that("stored and named missing-site geometry agree with random effects", {
  set.seed(2)
  for (family in list(NULL, "poisson")) {
    fit <- conditional_iid_fixture(family, missing = TRUE, random = TRUE)
    fit$ssn.object$preds$alias <- fit$ssn.object$preds$.missing
    for (local in list(FALSE,
      list(approximation = "vecchia", size = 3, ordering = "coordinate", method_random = "random"),
      list(approximation = "low-rank", size_base = 8, reorder_base = "none", method_new = "all"))) {
      set.seed(2)
      stored <- conditional(fit, samples = 3, local = local)
      set.seed(2)
      named <- conditional(fit, "alias", samples = 3, local = local)
      expect_equal(unname(stored), unname(named), tolerance = 1e-12)
      expect_identical(rownames(stored), as.character(fit$missing_index))
    }
  }
})

test_that("zero stored GLM nugget retains the effective variance floor", {
  set.seed(2)
  fit <- conditional_iid_fixture("poisson", variance = 0)
  expect_equal(unname(coef(fit, "nugget")[["nugget"]]), 0)
  expect_equal(get_conditional_iid_variance(fit), fit$diagtol)
  outputs <- c("beta", "latent", "newdata")
  set.seed(2)
  dense <- draw_conditional_exact(fit, "iid", 1003, outputs)
  set.seed(2)
  value <- conditional(fit, "iid", samples = 1003, output = outputs)
  expect_identical(value$beta, dense$beta)
  expect_equal(value, dense[names(value)], tolerance = 1e-12)
})
