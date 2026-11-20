skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("nested slope effect conditionals retain true partitions", {
  set.seed(2)
  observed <- data.frame(g = c("a", "a", "a", "a", "b", "b", "a", "b"),
    h = c("u", "u", "v", "v", "u", "u", "u", "v"),
    part = rep(c("P", "Q"), 4), x = c(0, 2, -1, 3, 2, 1, 0, -2))
  prediction <- data.frame(g = c("a", "b", "new", "a"), h = c("u", "u", "v", "v"),
    part = c("P", "P", "Q", "Q"), x = c(1, 2, -1, 0))
  params <- c("1 | g" = .3, "1 | g:h" = .2, "x | g:h" = .15, "1 | h" = 0)
  object <- list(ssn.object = list(obs = observed), random = ~g / h + (0 + x | g:h),
                 coefficients = list(params_object = list(randcov = params)))
  all <- rbind(observed, prediction)
  same_g <- outer(all$g, all$g, "==")
  same_h <- outer(all$h, all$h, "==")
  same_part <- outer(all$part, all$part, "==")
  prior <- same_part * same_g * (.3 + same_h * (.2 + .15 * outer(all$x, all$x)))
  spatial <- (diag(.5, 8) + .03 * outer(1:8, 1:8)) * same_part[1:8, 1:8]
  for (scenario in c("exact", "reduced")) {
    index <- if (scenario == "reduced") 1:4 else 1:8
    design <- get_randcov_design(observed, params, ~part, prediction, index)
    Z <- as.matrix(design$observed)
    Zp <- as.matrix(design$newdata)
    combined <- rbind(Z, Zp)
    support <- c(index, 9:12)
    expect_equal(combined %*% diag(design$variance) %*% t(combined),
                 prior[support, support], ignore_attr = TRUE, tolerance = 1e-12)
    S <- spatial[index, index]
    effect <- prepare_conditional_random(object, design, covariance = S)
    posterior <- diag(design$variance)
    posterior[effect$active, effect$active] <- effect$covariance
    Rpo <- prior[9:12, index]
    marginal <- S + prior[index, index]
    expected_cov <- prior[9:12, 9:12] - Rpo %*% solve(marginal, t(Rpo))
    expect_equal(Zp %*% posterior %*% t(Zp), expected_cov, ignore_attr = TRUE, tolerance = 1e-10)
    residual <- cbind(seq_along(index) / 5, rev(seq_along(index)) / 3)
    rhs <- matrix(0, length(effect$active), 2)
    for (block in effect$blocks) {
      rhs[block$cols, ] <- rhs[block$cols, , drop = FALSE] +
        crossprod(block$weights, residual[block$rows, , drop = FALSE])
    }
    mean <- Zp[, effect$active, drop = FALSE] %*% effect$covariance %*% rhs
    expect_equal(mean, Rpo %*% solve(marginal, residual), ignore_attr = TRUE, tolerance = 1e-10)
    set.seed(2)
    noise <- matrix(rnorm(length(design$variance) * 2), length(design$variance), 2)
    set.seed(2)
    draws <- draw_conditional_random(effect, residual)
    inactive <- setdiff(seq_along(design$variance), effect$active)
    expect_equal(draws[inactive, , drop = FALSE],
                 sqrt(design$variance[inactive]) * noise[inactive, , drop = FALSE])
    expect_true(all(draws[design$variance == 0, , drop = FALSE] == 0))
    centered <- draws
    centered[effect$active, ] <- centered[effect$active, ] - effect$factor %*% noise[effect$active, ]
    centered[inactive, ] <- 0
    expect_equal(Zp %*% centered, mean, ignore_attr = TRUE, tolerance = 1e-10)
    if (scenario == "reduced") {
      expect_true(any(design$info$observed & !design$info$conditioned & design$variance > 0))
      expect_true(any(!design$info$observed & !design$info$conditioned & design$variance > 0))
    }
  }
})

conditional_stream_network <- function() {
  network <- mf04p
  select <- function(data, count) {
    unlist(lapply(split(seq_len(NROW(data)), data$netID), function(rows) {
      rows[unique(round(seq(1, length(rows), length.out = min(count, length(rows)))))]
    }), use.names = FALSE)
  }
  network$obs <- network$obs[rev(select(network$obs, 8)), ]
  network$preds <- list(stream_test = network$preds$pred1km[rev(select(network$preds$pred1km, 3)), ])
  network$obs$x <- seq(-1, 1, length.out = NROW(network$obs))
  network$obs$y <- rep(c(1, 3, 2, 4), length.out = NROW(network$obs))
  network$preds$stream_test$x <- seq(-.5, .5, length.out = NROW(network$preds$stream_test))
  network
}

conditional_stream_fit <- function(network, component) {
  ssn_lm(y ~ x, network, additive = "afvArea", ddf = "asymptotic",
    tailup_initial = if (component %in% c("up", "mixed"))
      tailup_initial("exponential", .4, 10000, known = "given") else tailup_initial("none"),
    taildown_initial = if (component %in% c("down", "mixed"))
      taildown_initial("exponential", .3, 10000, known = "given") else taildown_initial("none"),
    euclid_initial = if (component == "euclid")
      euclid_initial("exponential", .2, 12000, known = "given") else euclid_initial("none"),
    nugget_initial = nugget_initial("nugget", 0, known = "given"))
}

test_that("isolated stream components preserve multi-network identities", {
  set.seed(2)
  network <- conditional_stream_network()
  root <- withr::local_tempdir()
  stopifnot(file.copy(network$path, root, recursive = TRUE, copy.mode = FALSE))
  network$path <- file.path(root, basename(network$path))
  unlink(file.path(network$path, "distance"), recursive = TRUE)
  ssn_create_distmat(network, predpts = "stream_test", among_predpts = TRUE, overwrite = TRUE)
  networks <- c(as.character(network$obs$netID), as.character(network$preds$stream_test$netID))
  expect_gt(length(unique(network$obs$netID)), 1)
  expect_gt(length(unique(network$preds$stream_test$netID)), 1)
  references <- list()
  fits <- lapply(c("up", "down", "mixed"), function(component) conditional_stream_fit(network, component))
  for (k in seq_along(fits)) {
    fit <- fits[[k]]
    obs <- as.matrix(covmatrix(fit))
    cross <- as.matrix(covmatrix(fit, "stream_test", cov_type = "obs.pred"))
    pred <- as.matrix(covmatrix(fit, "stream_test", cov_type = "pred.pred"))
    references[[k]] <- rbind(cbind(obs, cross), cbind(t(cross), pred))
    expect_true(all(references[[k]][outer(networks, networks, "!=")] == 0))
  }
  unconnected <- references[[1]] == 0 & references[[2]] > 0 & outer(networks, networks, "==")
})
