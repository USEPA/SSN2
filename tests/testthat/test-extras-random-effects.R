skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# randcov slope
test_that("point prediction se.fit is heteroscedastic for a random slope and matches a manual reference (Bug A)", {
  fit <- ssn_lm(Summer_mn ~ 1, mf04p,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    random = ~ (ELEV_DEM | as.factor(netID))
  )

  preds <- fit$ssn.object$preds$CapeHorn
  pred_pt <- predict(fit, "CapeHorn", se.fit = TRUE)

  # heteroscedasticity: se.fit must vary across rows now that it reflects
  # each row's own ELEV_DEM value, not a single constant
  expect_true(length(unique(round(pred_pt$se.fit, 8))) > 1)

  # manual reference at the min- and max-ELEV_DEM prediction rows, replicating
  # get_pred()'s kriging math directly from the fitted object's own pieces
  randcov <- fit$coefficients$params_object$randcov
  randcov_intercept <- as.numeric(randcov["1 | as.factor(netID)"])
  randcov_slope <- as.numeric(randcov["ELEV_DEM | as.factor(netID)"])
  nugget <- fit$coefficients$params_object$nugget[["nugget"]]
  tailup_de <- fit$coefficients$params_object$tailup[["de"]]

  cov_vector <- covmatrix(fit, "CapeHorn")
  cov_matrix_val <- covmatrix(fit)
  cov_lowchol <- t(chol(cov_matrix_val))
  Xmat <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  cov_betahat <- vcov(fit)

  manual_se <- function(i) {
    c0 <- as.numeric(cov_vector[i, ])
    SqrtSigInv_X <- forwardsolve(cov_lowchol, Xmat)
    SqrtSigInv_c0 <- forwardsolve(cov_lowchol, c0)
    x0 <- matrix(1, nrow = 1, ncol = 1)
    H <- x0 - crossprod(SqrtSigInv_c0, SqrtSigInv_X)
    total_var <- nugget + tailup_de + randcov_intercept + randcov_slope * preds$ELEV_DEM[i]^2
    var <- as.numeric(total_var - crossprod(SqrtSigInv_c0, SqrtSigInv_c0) + H %*% tcrossprod(cov_betahat, H))
    sqrt(var)
  }

  i_min <- which.min(preds$ELEV_DEM)
  i_max <- which.max(preds$ELEV_DEM)
  expect_equal(pred_pt$se.fit[[i_min]], manual_se(i_min), tolerance = 1e-8)
  expect_equal(pred_pt$se.fit[[i_max]], manual_se(i_max), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(pred_pt$se.fit[[i_min]], pred_pt$se.fit[[i_max]])))
})

test_that("random-intercept-only prediction variance is unchanged (invariant that keeps existing coverage byte-identical)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    random = ~ as.factor(netID)
  )

  preds <- fit$ssn.object$preds$CapeHorn
  row_a <- preds[1, , drop = FALSE]
  row_b <- preds[which.max(preds$ELEV_DEM), , drop = FALSE]
  randcov_contribution_a <- randcov_newvar(fit$coefficients$params_object$randcov, row_a)
  randcov_contribution_b <- randcov_newvar(fit$coefficients$params_object$randcov, row_b)
  expect_equal(randcov_contribution_a, as.numeric(fit$coefficients$params_object$randcov))
  expect_equal(randcov_contribution_a, randcov_contribution_b)

  spatial_nugget_var <- get_spatial_nugget_var(fit$coefficients$params_object, fit$diagtol)
  cov_matrix_val <- covmatrix(fit)
  expect_equal(spatial_nugget_var + randcov_contribution_a, cov_matrix_val[1, 1])
})

test_that("an NA in a random-slope covariate in newdata errors instead of silently corrupting predictions (Bug B)", {
  fit <- ssn_lm(Summer_mn ~ 1, mf04p,
    tailup_type = "none", taildown_type = "none", euclid_type = "none",
    additive = "afvArea", random = ~ (ELEV_DEM | as.factor(netID))
  )

  broken <- fit
  preds <- broken$ssn.object$preds$CapeHorn
  preds$ELEV_DEM[2] <- NA
  broken$ssn.object$preds$CapeHorn <- preds

  expect_error(
    predict(broken, "CapeHorn", se.fit = TRUE),
    "Cannot have NA values in predictors"
  )
})

test_that("block prediction for a random-slope model still runs correctly (predict_block()'s own pred.pred-based total_var is untouched by this fix)", {
  fit <- ssn_lm(Summer_mn ~ 1, mf04p,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    random = ~ (ELEV_DEM | as.factor(netID))
  )

  block_pred <- predict(fit, "CapeHorn", block = TRUE, se.fit = TRUE)
  expect_true(is.finite(block_pred$fit))
  expect_true(is.finite(block_pred$se.fit) && block_pred$se.fit > 0)

  expect_true(is.finite(mean(covmatrix(fit, "CapeHorn", cov_type = "pred.pred"))))
})

test_that("get_pred_glm() uses randcov_newvar() too, via the same shared helper as the LM path (Bug A, GLM path)", {
  randcov_params <- c("1 | grp" = 2, "x | grp" = 3)
  row_a <- data.frame(grp = "a", x = 2)
  row_b <- data.frame(grp = "a", x = 5)
  expect_equal(randcov_newvar(randcov_params, row_a), 2 + 3 * 2^2)
  expect_equal(randcov_newvar(randcov_params, row_b), 2 + 3 * 5^2)
  expect_false(isTRUE(all.equal(
    randcov_newvar(randcov_params, row_a),
    randcov_newvar(randcov_params, row_b)
  )))

  s <- mf04p
  s$obs$count_response <- rpois(nrow(s$obs), lambda = 5)

  fit <- ssn_glm(count_response ~ 1, s, family = "poisson",
    tailup_type = "exponential", additive = "afvArea",
    random = ~ (ELEV_DEM | as.factor(netID))
  )

  pred_pt <- predict(fit, "CapeHorn", se.fit = TRUE)
  expect_true(all(is.finite(pred_pt$se.fit)) && all(pred_pt$se.fit > 0))

  expect_true(fit$diagtol > 0)
  spatial_nugget_var <- get_spatial_nugget_var(fit$coefficients$params_object, fit$diagtol)
  obs_row_1 <- fit$ssn.object$obs[1, , drop = FALSE]
  randcov_contribution <- randcov_newvar(fit$coefficients$params_object$randcov, obs_row_1)
  cov_matrix_val <- covmatrix(fit)
  expect_equal(spatial_nugget_var + randcov_contribution, cov_matrix_val[1, 1])
})

# randcov xlev
test_that("block prediction with a random-effect level absent from the prediction subset no longer crashes and computes a correct random design matrix", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    random = ~ as.factor(netID)
  )
  netgeom_pred <- ssn_get_netgeom(fit$ssn.object$preds$CapeHorn, reformat = TRUE)
  expect_true(all(netgeom_pred$NetworkID == 2))
  netgeom_obs <- ssn_get_netgeom(fit$ssn.object$obs, reformat = TRUE)
  expect_true(length(unique(netgeom_obs$NetworkID)) > 1)

  # the fix captures fitted random-effect levels on the object
  expect_named(fit$random_xlev, get_randcov_names(fit$random))
  expect_true(length(fit$random_xlev[[1]][["as.factor(netID)"]]) > 1)

  # point prediction was already correct before this fix and must remain so
  point_pred <- predict(fit, "CapeHorn")
  expect_true(all(is.finite(point_pred)))

  broken_fit <- fit
  broken_fit$random_xlev <- NULL
  expect_error(
    predict(broken_fit, "CapeHorn", block = TRUE, se.fit = TRUE),
    "contrasts can be applied only to factors with 2 or more levels"
  )

  block_pred <- predict(fit, "CapeHorn", block = TRUE, se.fit = TRUE)
  expect_true(is.finite(block_pred$fit))
  expect_true(is.finite(block_pred$se.fit) && block_pred$se.fit > 0)

  randcov_Zs <- get_randcov_Zs(
    fit$ssn.object$preds$CapeHorn, get_randcov_names(fit$random),
    xlev_list = fit$random_xlev
  )
  Z <- as.matrix(randcov_Zs[[1]]$Z)
  expect_equal(ncol(Z), 1)
  expect_true(all(Z == 1))
})

test_that("block prediction is unaffected when every training random-effect level is present in the prediction subset", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    random = ~ as.factor(netID)
  )
  # pred1km (packaged fixture) spans the same networks as the observed data
  netgeom_obs <- ssn_get_netgeom(fit$ssn.object$obs, reformat = TRUE)
  netgeom_pred <- ssn_get_netgeom(fit$ssn.object$preds$pred1km, reformat = TRUE)
  expect_true(setequal(unique(netgeom_obs$NetworkID), unique(netgeom_pred$NetworkID)))

  block_pred <- predict(fit, "pred1km", block = TRUE, se.fit = TRUE)
  expect_true(is.finite(block_pred$fit))
  expect_true(is.finite(block_pred$se.fit) && block_pred$se.fit > 0)
})

# random covariance
ssn_create_bigdist(mf04p, predpts = "CapeHorn", overwrite = TRUE, no_cores = 1, verbose = FALSE)

test_that("multi-group theoretical variance with a known random intercept matches an independent sandwich", {
  set.seed(2)
  s <- mf04p
  s$obs <- s$obs[sample(nrow(s$obs)), ]
  s$obs$audit_group <- factor(rep(1:4, length.out = nrow(s$obs)))

  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  td <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  ng <- nugget_initial("nugget", nugget = 0.5, known = "nugget")
  rc <- randcov_initial(audit_group = 1, known = "given")
  group_index <- rep(c(2, 7, 10), length.out = nrow(s$obs))

  dense <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea",
    random = ~audit_group, randcov_initial = rc
  )
  grouped <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea",
    local = list(index = group_index, var_adjust = "theoretical"),
    random = ~audit_group, randcov_initial = rc
  )

  full <- covmatrix(dense)
  B <- matrix(0, nrow(full), ncol(full))
  for (ids in split(seq_len(nrow(full)), group_index)) B[ids, ids] <- solve(full[ids, ids])
  XX <- model.matrix(dense)
  AA <- solve(crossprod(XX, B %*% XX)) %*% t(XX) %*% B
  sandwich <- AA %*% full %*% t(AA)

  expect_equal(as.matrix(vcov(grouped)), sandwich, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("multi-group theoretical variance with a random intercept agrees between serial and parallel", {
  set.seed(2)
  s <- mf04p
  s$obs$audit_group <- factor(rep(1:4, length.out = nrow(s$obs)))

  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  td <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  ng <- nugget_initial("nugget", nugget = 0.5, known = "nugget")
  rc <- randcov_initial(audit_group = 1, known = "given")
  group_index <- rep(c(2, 7, 10), length.out = nrow(s$obs))

  fit_serial <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea",
    local = list(index = group_index, var_adjust = "theoretical"),
    random = ~audit_group, randcov_initial = rc
  )
  fit_parallel <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea",
    local = list(index = group_index, var_adjust = "theoretical", parallel = TRUE, ncores = 2),
    random = ~audit_group, randcov_initial = rc
  )
  expect_equal(vcov(fit_serial), vcov(fit_parallel), tolerance = 1e-10)
})

test_that("prediction/observation cross covariance includes the random-effect term when every spatial component is none", {
  set.seed(2)
  s <- mf04p
  n_obs <- nrow(s$obs)
  n_pred <- nrow(s$preds$CapeHorn)
  s$obs$audit_group <- factor(rep(1:3, length.out = n_obs))
  s$preds$CapeHorn$audit_group <- factor(rep(1:3, length.out = n_pred))

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_type = "none", taildown_type = "none", euclid_type = "none", nugget_type = "nugget",
    random = ~audit_group
  )
  cv <- as.matrix(covmatrix(fit, "CapeHorn"))

  G <- as.numeric(coef(fit, type = "randcov"))
  Z_obs <- model.matrix(~ audit_group - 1, data = s$obs)
  Z_pred <- model.matrix(~ audit_group - 1, data = s$preds$CapeHorn)
  expected <- G * (Z_pred %*% t(Z_obs))

  expect_equal(cv, expected, tolerance = 1e-8, ignore_attr = TRUE)
  # a real regression guard, not just a numerical-tolerance check: the old
  # code returned exactly zero everywhere here
  expect_true(any(cv != 0))
})

test_that("prediction/observation cross covariance still applies the partition mask after the random-effect fix", {
  set.seed(2)
  s <- mf04p
  n_obs <- nrow(s$obs)
  n_pred <- nrow(s$preds$CapeHorn)
  s$obs$audit_group <- factor(rep(1:3, length.out = n_obs))
  s$preds$CapeHorn$audit_group <- factor(rep(1:3, length.out = n_pred))
  s$obs$audit_partition <- factor(rep(c("a", "b"), length.out = n_obs))
  s$preds$CapeHorn$audit_partition <- factor(rep(c("a", "b"), length.out = n_pred))

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_type = "none", taildown_type = "none", euclid_type = "none", nugget_type = "nugget",
    random = ~audit_group, partition_factor = ~audit_partition
  )
  cv <- as.matrix(covmatrix(fit, "CapeHorn"))

  diff_partition <- outer(as.character(s$preds$CapeHorn$audit_partition), as.character(s$obs$audit_partition), "!=")
  same_group <- outer(as.character(s$preds$CapeHorn$audit_group), as.character(s$obs$audit_group), "==")

  expect_true(all(cv[diff_partition] == 0))
  expect_true(any(cv[!diff_partition & same_group] > 0))
})

# random effect grouping via interaction/nesting
test_that("get_randcov_Z() correctly crosses interaction/nested grouping variables instead of silently dropping all but the first (spmodel bug, fixed here via interaction() + fac2sparse())", {
  set.seed(2)
  data <- data.frame(
    g1 = factor(sample(c("a", "b", "c"), 30, replace = TRUE)),
    g2 = factor(sample(c("x", "y"), 30, replace = TRUE))
  )

  result <- get_randcov_Z("1 | g1:g2", data)

  # reference built the old (pre-optimization) way: a dense model.matrix() on
  # the interaction term, with unobserved combinations dropped identically to
  # get_randcov_Z()
  mm <- model.matrix(~ g1:g2 - 1, data)
  mm <- mm[, colSums(abs(mm)) > 0, drop = FALSE]

  expect_equal(ncol(result$Z), ncol(mm))
  expect_equal(as.matrix(result$ZZt), tcrossprod(mm, mm), ignore_attr = TRUE)

  # a real regression guard: with the Z_frame[[1]]-only bug, ncol(result$Z)
  # would equal nlevels(g1) (g2 silently dropped), not the number of observed
  # g1:g2 combinations
  expect_true(ncol(result$Z) > nlevels(data$g1))

  # one-hot: every row belongs to exactly one group
  expect_true(all(Matrix::rowSums(result$Z) == 1))
})

test_that("ssn_lm() fits a nested random effect (random = ~ g1/g2) and its random covariance structure matches an independently-built reference", {
  set.seed(2)
  s <- mf04p
  n_obs <- nrow(s$obs)
  s$obs$g1 <- factor(rep(c("n1", "n2"), length.out = n_obs))
  s$obs$g2 <- factor(rep(c("s1", "s2", "s3"), length.out = n_obs))

  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_type = "none", taildown_type = "none", euclid_type = "none", nugget_type = "nugget",
    random = ~ g1 / g2
  )

  expect_named(fit$coefficients$params_object$randcov, c("1 | g1", "1 | g1:g2"))

  cv <- as.matrix(covmatrix(fit))
  G <- coef(fit, type = "randcov")
  Z_g1 <- model.matrix(~g1 - 1, data = s$obs)
  Z_g1g2 <- model.matrix(~ g1:g2 - 1, data = s$obs)
  Z_g1g2 <- Z_g1g2[, colSums(abs(Z_g1g2)) > 0, drop = FALSE]
  nugget <- fit$coefficients$params_object$nugget[["nugget"]]
  expected <- as.numeric(G["1 | g1"]) * tcrossprod(Z_g1, Z_g1) +
    as.numeric(G["1 | g1:g2"]) * tcrossprod(Z_g1g2, Z_g1g2) +
    diag(nugget, n_obs)

  expect_equal(cv, expected, tolerance = 1e-8, ignore_attr = TRUE)
  # a real regression guard: with the g2-dropping bug, the "1 | g1:g2" term
  # collapses onto "1 | g1" and this comparison fails
  expect_true(any(Z_g1g2 %*% t(Z_g1g2) != Z_g1 %*% t(Z_g1)))
})

test_that("theoretical variance with no random effects retains previously observed agreement (unaffected by F5)", {
  set.seed(2)
  s <- mf04p
  s$obs <- s$obs[sample(nrow(s$obs)), ]
  group_index <- rep(c(2, 7, 10), length.out = nrow(s$obs))

  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  td <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  ng <- nugget_initial("nugget", nugget = 0.5, known = "nugget")

  dense <- ssn_lm(Summer_mn ~ ELEV_DEM, s, tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea")
  grouped <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng, additive = "afvArea",
    local = list(index = group_index, var_adjust = "theoretical")
  )

  full <- covmatrix(dense)
  B <- matrix(0, nrow(full), ncol(full))
  for (ids in split(seq_len(nrow(full)), group_index)) B[ids, ids] <- solve(full[ids, ids])
  XX <- model.matrix(dense)
  AA <- solve(crossprod(XX, B %*% XX)) %*% t(XX) %*% B
  sandwich <- AA %*% full %*% t(AA)

  expect_equal(as.matrix(vcov(grouped)), sandwich, tolerance = 1e-6, ignore_attr = TRUE)
})


test_that("get_randcov_context()'s reusable pieces are built once per prediction operation, not once per chunk", {
  # stream covariance must be active (not "none") for point prediction's
  # covariance construction to take the chunked (.bmat) path at all -- a
  # pure nugget/random-effect fit resolves to the unchunked dense path
  # instead, which never calls get_randcov_context() per chunk in the
  # first place
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, random = ~ as.factor(netID),
    tailup_type = "exponential", additive = "afvArea"
  )

  call_count <- 0
  real_context <- get_randcov_context
  counting_context <- function(...) {
    call_count <<- call_count + 1
    real_context(...)
  }
  testthat::local_mocked_bindings(get_randcov_context = counting_context, .package = "SSN2")

  # a tiny chunk_size forces several chunks over CapeHorn's prediction rows
  predict(fit, "CapeHorn", local = list(method = "all", chunk_size = 3))
  expect_equal(call_count, 1)
})

test_that("block prediction's diagonal-sum context is built once, not once per get_block_diagonal_sum() call", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, random = ~ as.factor(netID),
    tailup_type = "exponential", additive = "afvArea"
  )

  call_count <- 0
  real_context <- get_randcov_context
  counting_context <- function(...) {
    call_count <<- call_count + 1
    real_context(...)
  }
  testthat::local_mocked_bindings(get_randcov_context = counting_context, .package = "SSN2")

  # method_new = "subset" with size_new below the grid size makes
  # get_block_diagonal_sum() run twice (once for the node subset, once for
  # the full grid); se.fit = TRUE makes it run at all
  predict(fit, "CapeHorn",
    block = TRUE, se.fit = TRUE,
    local = list(method = "all", method_new = "subset", size_new = 5L, ordering = "pid", chunk_size = 1000L, parallel = FALSE)
  )
  expect_equal(call_count, 1)
})

test_that("randcov context reuse produces identical results for a random-slope term, matching fresh derivation across row-chunks", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "none", taildown_type = "none", euclid_type = "none",
    random = ~ (ELEV_DEM | as.factor(netID))
  )
  randcov_params <- fit$coefficients$params_object$randcov

  data <- mf04p$obs
  newdata_all <- mf04p$preds$CapeHorn

  # reference: no context, computed directly on the full newdata set
  reference <- randcov_vector(randcov_params, data, newdata_all)

  # cached: one context built once, reused across row-chunks
  context <- get_randcov_context(randcov_params, data, newdata_all)
  chunks <- split(seq_len(nrow(newdata_all)), ceiling(seq_len(nrow(newdata_all)) / 4))
  chunked <- do.call(rbind, lapply(chunks, function(rows) {
    randcov_vector(randcov_params, data, newdata_all[rows, , drop = FALSE], context = context)
  }))

  expect_equal(as.matrix(chunked), as.matrix(reference), ignore_attr = TRUE)

  # a genuine slope term: entries are not all equal to the variance
  # component alone, the way they would be for an intercept-only term
  expect_true(length(unique(as.numeric(reference))) > 2)
})

test_that("point prediction's parallel dispatch matches serial output for a random-slope model", {
  # regression test: get_pred()'s call to randcov_newvar() passes a context
  # argument that only the current development source's randcov_newvar()
  # accepts. A worker process resolves that call through the SSN2 namespace
  # reconstructed on the worker, not through clusterExport()-ed bindings in
  # its global environment, so run_pred_dispatch() must reload SSN2's
  # current source (pkgload::load_all()) on every worker before dispatching,
  # the same way get_block_chunk_apply() already does -- omitting that step
  # fails with "unused argument (context = randcov_context)" as soon as any
  # random-effect term reaches a worker, whether or not it is actually used.
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    random = ~ (ELEV_DEM | as.factor(netID))
  )

  serial <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "all", parallel = FALSE))
  parallel_res <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "all", parallel = TRUE, ncores = 2))
  expect_equal(parallel_res, serial)
})

test_that("get_randcov_context()'s formulas do not retain dense training-side temporaries when serialized", {
  # Serialization exposes retained environments that object.size() omits.
  serialized_bytes <- function(n) {
    d <- data.frame(g = factor(rep(seq_len(n / 2), each = 2)), x = seq_len(n))
    params <- setNames(1, "x | g")
    context <- get_randcov_context(params, d, d[1:10, ])
    length(serialize(context, NULL))
  }

  small <- serialized_bytes(250)
  large <- serialized_bytes(2000)

  expect_lt(large, small * 20)
  expect_lt(large, 2e5)

  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 3)), x = 1:9)
  params <- setNames(1, "x | as.factor(g)")
  context <- get_randcov_context(params, d, d)
  restored <- unserialize(serialize(context, NULL))
  for (term in list(context[[1]], restored[[1]])) {
    for (name in c("reform_bar1", "reform_bar2")) {
      expect_identical(environment(term[[name]]), asNamespace("SSN2"))
      expect_false(exists("Z_index_data_mx", envir = environment(term[[name]]), inherits = FALSE))
      expect_false(exists("Z_index_data_split", envir = environment(term[[name]]), inherits = FALSE))
    }
  }
  with_context <- get_randcov_vectors("x | as.factor(g)", params, d, d, context = context)
  without_context <- get_randcov_vectors("x | as.factor(g)", params, d, d, context = NULL)
  expect_equal(with_context, without_context)
  expect_equal(randcov_newvar(params, d[1, , drop = FALSE], context = context), params[[1]] * d$x[1]^2)
})

local_randcov_function_masks <- function(.local_envir = parent.frame()) {
  bindings <- list(
    log = function(x, base = exp(1)) base::log(x, base = base) + 0.001,
    as.factor = function(...) stop("global grouping function used"),
    pnorm = function(...) stop("global transformation function used")
  )
  existing <- intersect(names(bindings), ls(envir = globalenv(), all.names = TRUE))
  saved <- mget(existing, envir = globalenv(), inherits = FALSE)
  withr::defer({
    rm(list = names(bindings), envir = globalenv())
    for (name in names(saved)) assign(name, saved[[name]], envir = globalenv())
  }, envir = .local_envir)
  for (name in names(bindings)) assign(name, bindings[[name]], envir = globalenv())
  invisible(NULL)
}

test_that("random-effect formula lookup agrees with fitting despite global function masks", {
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 3)), x = 1:9)
  local_randcov_function_masks()
  for (name in c("1 | as.factor(g)", "log(x) | as.factor(g)", "pnorm(x) | as.factor(g)")) {
    params <- setNames(1, name)
    reference <- as.matrix(get_randcov_Z(name, d)$ZZt)
    context <- get_randcov_context(params, d, d)
    for (ctx in list(NULL, context, unserialize(serialize(context, NULL)))) {
      actual <- as.matrix(randcov_vector(params, d, d, context = ctx))
      expect_equal(unname(actual), unname(reference), tolerance = 1e-12)
      variances <- vapply(seq_len(nrow(d)), function(i) {
        randcov_newvar(params, d[i, , drop = FALSE], context = ctx)
      }, numeric(1))
      expect_equal(variances, unname(diag(reference)), tolerance = 1e-12)
    }
  }
})

test_that("transformed random-slope predictions are stable under global function masks", {
  fits <- list(
    ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
      tailup_type = "exponential", additive = "afvArea",
      random = ~ (log(ELEV_DEM) | as.factor(netID))
    ),
    ssn_glm(Summer_mn ~ ELEV_DEM, mf04p, family = "Gamma",
      tailup_type = "exponential", additive = "afvArea",
      random = ~ (log(ELEV_DEM) | as.factor(netID))
    )
  )
  settings <- list(list(method = "all"), list(method = "covariance", size = 20))
  baseline <- lapply(fits, function(fit) lapply(settings, function(local) {
    predict(fit, "CapeHorn", se.fit = TRUE, local = local)
  }))
  local_randcov_function_masks()
  for (i in seq_along(fits)) {
    for (j in seq_along(settings)) {
      local <- settings[[j]]
      serial <- predict(fits[[i]], "CapeHorn", se.fit = TRUE, local = local)
      parallel_result <- predict(fits[[i]], "CapeHorn", se.fit = TRUE,
        local = c(local, list(parallel = TRUE, ncores = 2))
      )
      expect_true(all(is.finite(serial$se.fit)))
      expect_true(all(is.finite(parallel_result$se.fit)))
      expect_equal(serial, baseline[[i]][[j]])
      expect_equal(parallel_result, serial)
    }
  }
})

test_that("model_matrix_group_labels() matches an independently computed group label for single and crossed-variable formulas", {
  d <- data.frame(g = factor(c("a", "a", "b", "b", "b")), h = factor(c("x", "y", "x", "y", "x")), x = 1:5)

  independent_labels <- function(reform, data) {
    mf <- model.frame(reform, data)
    mx <- model.matrix(reform, mf)
    vapply(seq_len(nrow(mx)), function(i) colnames(mx)[which(as.logical(mx[i, ]))], character(1))
  }

  # a single combined grouping term in each case -- an interaction (g:h), not
  # an additive multi-term formula (g + h), matching how this helper is
  # actually used elsewhere (one grouping/partition formula per call); an
  # additive formula does not one-hot encode to a single column per row once
  # more than one factor is involved, so it is not a supported input here
  for (reform in list(reformulate("g", intercept = FALSE), reformulate("g:h", intercept = FALSE))) {
    expect_equal(model_matrix_group_labels(reform, d), independent_labels(reform, d))
  }

  # A missing grouping value cannot identify a single indicator column.
  d_na <- d
  d_na$g[2] <- NA
  expect_error(model_matrix_group_labels(reformulate("g", intercept = FALSE), d_na, na_pass = TRUE))

  # xlev extends the label vocabulary without changing labels already present
  xlev <- list(g = c("a", "b", "c"))
  expect_equal(
    model_matrix_group_labels(reformulate("g", intercept = FALSE), d, xlev = xlev),
    independent_labels(reformulate("g", intercept = FALSE), d)
  )
})

test_that("get_partition_context()'s reusable pieces are built once per point-prediction operation, not once per chunk", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    partition_factor = ~ as.factor(netID)
  )

  call_count <- 0
  real_context <- get_partition_context
  counting_context <- function(...) {
    call_count <<- call_count + 1
    real_context(...)
  }
  testthat::local_mocked_bindings(get_partition_context = counting_context, .package = "SSN2")

  # a tiny chunk_size forces several chunks over CapeHorn's prediction rows
  predict(fit, "CapeHorn", local = list(method = "all", chunk_size = 3))
  expect_equal(call_count, 1)
})

test_that("point prediction without a partition factor never builds a partition context", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential", additive = "afvArea")

  call_count <- 0
  real_context <- get_partition_context
  counting_context <- function(...) {
    call_count <<- call_count + 1
    real_context(...)
  }
  testthat::local_mocked_bindings(get_partition_context = counting_context, .package = "SSN2")

  predict(fit, "CapeHorn", local = list(method = "all", chunk_size = 3))
  expect_equal(call_count, 0)
})

test_that("partition context reuse produces identical predictions across chunks, matching a single unchunked call", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    partition_factor = ~ as.factor(netID)
  )

  unchunked <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "all", chunk_size = 1000L))
  chunked <- predict(fit, "CapeHorn", se.fit = TRUE, local = list(method = "all", chunk_size = 3L))
  expect_equal(chunked, unchunked)

  # partition_vector() itself: a precomputed context built from the full
  # prediction set matches a fresh, uncontexted call for a row subset --
  # the same equivalence get_randcov_context()'s own test establishes
  data <- mf04p$obs
  newdata_all <- mf04p$preds$CapeHorn
  reference <- partition_vector(fit$partition_factor, data, newdata_all)
  context <- get_partition_context(fit$partition_factor, data, newdata_all)
  chunks <- split(seq_len(nrow(newdata_all)), ceiling(seq_len(nrow(newdata_all)) / 4))
  chunked_partition <- do.call(rbind, lapply(chunks, function(rows) {
    partition_vector(
      fit$partition_factor, data, newdata_all[rows, , drop = FALSE],
      reform_bar2 = context$reform_bar2, partition_index_data = context$partition_index_data
    )
  }))
  expect_equal(as.matrix(chunked_partition), as.matrix(reference), ignore_attr = TRUE)
})

test_that("get_partition_context()'s level_index_map is built once and actually reused by partition_vector(), not silently rebuilt", {
  # regression test: an earlier version of get_partition_context() cached
  # reform_bar2_vals/reform_bar2_xlev but not the group-label-to-training-row
  # lookup map, so partition_vector() unconditionally rebuilt it
  # (split(seq_along(group_label), group_label)) on every call regardless of
  # what was supplied -- this failed to detect that, since every other cached
  # piece still produced a correct result. Confirming the map field exists
  # would have caught the missing piece directly; deliberately corrupting a
  # supplied map and confirming the result changes confirms partition_vector()
  # actually consumes it rather than quietly recomputing its own copy.
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea",
    partition_factor = ~ as.factor(netID)
  )
  data <- mf04p$obs
  newdata_all <- mf04p$preds$CapeHorn

  context <- get_partition_context(fit$partition_factor, data, newdata_all)
  expect_false(is.null(context$partition_index_data$level_index_map))

  reference <- partition_vector(fit$partition_factor, data, newdata_all)
  reused <- partition_vector(
    fit$partition_factor, data, newdata_all,
    reform_bar2 = context$reform_bar2, partition_index_data = context$partition_index_data
  )
  expect_equal(as.matrix(reused), as.matrix(reference), ignore_attr = TRUE)

  corrupted_index_data <- context$partition_index_data
  corrupted_index_data$level_index_map <- lapply(corrupted_index_data$level_index_map, function(x) integer(0))
  corrupted_result <- partition_vector(
    fit$partition_factor, data, newdata_all,
    reform_bar2 = context$reform_bar2, partition_index_data = corrupted_index_data
  )
  expect_true(all(as.matrix(corrupted_result) == 0))
  expect_false(all(as.matrix(reference) == 0))

  # uncached fallback preserved: partition_index_data without a
  # level_index_map field still produces the correct result
  no_map_index_data <- context$partition_index_data
  no_map_index_data$level_index_map <- NULL
  fallback_result <- partition_vector(
    fit$partition_factor, data, newdata_all,
    reform_bar2 = context$reform_bar2, partition_index_data = no_map_index_data
  )
  expect_equal(as.matrix(fallback_result), as.matrix(reference), ignore_attr = TRUE)
})

test_that("get_partition_context()'s formulas do not retain dense training-side temporaries when serialized", {
  serialized_bytes <- function(n) {
    d <- data.frame(g = factor(rep(seq_len(n / 2), each = 2)), x = seq_len(n))
    context <- get_partition_context(~g, d, d[1:10, ])
    length(serialize(context, NULL))
  }

  small <- serialized_bytes(250)
  large <- serialized_bytes(2000)

  # a retained dense n x (n/2) group-indicator matrix would grow roughly with
  # the square of n; a cached level_index_map (one integer per row, spread
  # across groups) grows close to linearly -- well under quadratic
  expect_lt(large, small * 20)

  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 3)), x = 1:9)
  context <- get_partition_context(~g, d, d)
  expect_identical(environment(context$reform_bar2), asNamespace("SSN2"))
  expect_false(exists("p_index_data_mf", envir = environment(context$reform_bar2), inherits = FALSE))
})

# Sparse partition indicator
partition_matrix_reference <- function(data, formula) {
  mf <- model.frame(reformulate(labels(terms(formula)), intercept = FALSE), data)
  Reduce(`&`, lapply(mf, function(x) {
    outer(as.character(x), as.character(x), `==`)
  })) * 1
}

test_that("compound partition labels cannot merge distinct groups", {
  d <- data.frame(
    g = factor(c("a", "a.b", "a", "a.b", "a"), levels = c("a", "a.b", "unused")),
    h = factor(c("b.c", "c", "b.c", "b.c", "c"))
  )
  cases <- list(
    d,
    transform(d, g = as.character(g), h = as.character(h)),
    transform(d, g = ordered(g), h = ordered(h)),
    transform(d, k = factor(c("x", "x", "x", "y", "x")))
  )
  for (data in cases) {
    rownames(data) <- paste0("obs", c(9, 2, 15, 4, 11))
    form <- if ("k" %in% names(data)) ~g:h:k else ~g:h
    expected <- partition_matrix_reference(data, form)
    dimnames(expected) <- list(rownames(data), rownames(data))
    design <- model.matrix(reformulate(labels(terms(form)), intercept = FALSE), data)
    actual <- partition_matrix(form, data)

    expect_s4_class(actual, "sparseMatrix")
    expect_equal(as.matrix(actual), expected)
    expect_equal(as.matrix(actual), tcrossprod(design))
    expect_equal(as.numeric(actual[1, 2]), 0)
    expect_equal(as.numeric(actual[1, 3]), 1)

    rows <- c(5, 2, 4, 1, 3)
    expect_equal(as.matrix(partition_matrix(form, data[rows, ])), expected[rows, rows])
  }
})

test_that("compound partition groups preserve zero cross-group covariance in a public fit", {
  ssn <- mf04p
  n <- nrow(ssn$obs)
  ssn$obs$g <- factor(rep(c("a", "a.b"), length.out = n))
  ssn$obs$h <- factor(rep(c("b.c", "c"), length.out = n))
  euclid <- euclid_initial("exponential", de = 2, range = 10000, known = "given")
  nugget <- nugget_initial("nugget", nugget = 1, known = "given")
  unpartitioned <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn,
    euclid_initial = euclid, nugget_initial = nugget, ddf = "asymptotic"
  )
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn,
    euclid_initial = euclid, nugget_initial = nugget,
    partition_factor = ~g:h, ddf = "asymptotic"
  )
  mask <- partition_matrix_reference(ssn$obs, ~g:h)
  expected <- as.matrix(covmatrix(unpartitioned)) * mask
  actual <- as.matrix(covmatrix(fit))
  expect_equal(actual, expected)
  expect_true(all(actual[mask == 0] == 0))
  expect_gt(sum(actual[mask == 1 & row(actual) != col(actual)]), 0)

  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  coefficients <- solve(crossprod(X, solve(expected, X)), crossprod(X, solve(expected, y)))
  expect_equal(unname(coef(fit)), as.numeric(coefficients), tolerance = 1e-8)
})

test_that("partition_matrix() built directly as sparse matches an independent reference", {
  d <- data.frame(g = factor(c("a", "a", "b", "c", "c", "c"), levels = c("a", "b", "c")))
  row.names(d) <- c("obs1", "obs2", "obs3", "obs4", "obs5", "obs6")

  pm <- SSN2:::partition_matrix(~g, d)
  expected <- partition_matrix_reference(d, ~g)
  dimnames(expected) <- list(row.names(d), row.names(d))
  expect_equal(as.matrix(pm), expected)
  expect_true(methods::is(pm, "sparseMatrix"))

  # reordered rows preserve row identity, not positional order
  d_reordered <- d[c(6, 1, 4, 2, 5, 3), , drop = FALSE]
  pm_reordered <- SSN2:::partition_matrix(~g, d_reordered)
  expected_reordered <- partition_matrix_reference(d_reordered, ~g)
  dimnames(expected_reordered) <- list(row.names(d_reordered), row.names(d_reordered))
  expect_equal(as.matrix(pm_reordered), expected_reordered)

  # unused levels contribute all-zero columns internally but do not change
  # the returned group-membership matrix
  d_unused <- data.frame(g = factor(c("a", "a", "b"), levels = c("a", "b", "c")))
  pm_unused <- SSN2:::partition_matrix(~g, d_unused)
  expected_unused <- partition_matrix_reference(d_unused, ~g)
  dimnames(expected_unused) <- list(row.names(d_unused), row.names(d_unused))
  expect_equal(as.matrix(pm_unused), expected_unused)

  # a single row always has exactly one distinct value present, so it takes
  # the pre-existing single-group branch (unchanged by this construction,
  # dimnames always NULL) rather than the fac2sparse() branch -- true both
  # before and after this change
  d_single <- d[1, , drop = FALSE]
  pm_single_row <- SSN2:::partition_matrix(~g, d_single)
  expect_null(dimnames(pm_single_row)[[1]])
  expect_equal(as.matrix(pm_single_row), matrix(1, 1, 1), ignore_attr = TRUE)

  # crossed/interaction grouping terms
  d_crossed <- data.frame(g = factor(c("a", "a", "b", "c")), h = factor(c("x", "y", "x", "x")))
  pm_crossed <- SSN2:::partition_matrix(~ g:h, d_crossed)
  expected_crossed <- partition_matrix_reference(d_crossed, ~ g:h)
  dimnames(expected_crossed) <- list(row.names(d_crossed), row.names(d_crossed))
  expect_equal(as.matrix(pm_crossed), expected_crossed)

  # transformed grouping term
  d_transform <- data.frame(g = factor(c("a", "a", "b", "c")))
  pm_transform <- SSN2:::partition_matrix(reformulate("droplevels(g)"), d_transform)
  expected_transform <- partition_matrix_reference(d_transform, reformulate("droplevels(g)"))
  dimnames(expected_transform) <- list(row.names(d_transform), row.names(d_transform))
  expect_equal(as.matrix(pm_transform), expected_transform)

  # single-group case is unchanged: NULL dimnames, an all-ones matrix
  d_onegroup <- data.frame(g = factor(rep("only", 4)))
  pm_onegroup <- SSN2:::partition_matrix(~g, d_onegroup)
  expect_null(dimnames(pm_onegroup)[[1]])
  expect_null(dimnames(pm_onegroup)[[2]])
  expect_true(all(as.matrix(pm_onegroup) == 1))

  # rows with no explicit row.names() still get default (character-integer)
  # dimnames, matching the previous model.matrix()-based construction
  d_default_names <- data.frame(g = factor(c("a", "a", "b")))
  pm_default_names <- SSN2:::partition_matrix(~g, d_default_names)
  expect_equal(dimnames(pm_default_names), list(c("1", "2", "3"), c("1", "2", "3")))
})

test_that("partition_matrix() no longer builds a dense model.matrix() for multi-group data", {
  d <- data.frame(g = factor(c("a", "a", "b", "c", "c", "c")))
  called <- FALSE
  testthat::local_mocked_bindings(
    "model.matrix" = function(...) {
      called <<- TRUE
      stop("model.matrix() should not be called by partition_matrix()")
    },
    .package = "SSN2"
  )
  pm <- SSN2:::partition_matrix(~g, d)
  expect_false(called)
  expect_true(any(as.matrix(pm) == 0)) # sanity: this fixture actually has multiple groups
})

test_that("partition_matrix() results feed correctly into fitting, covmatrix(), and Torgegram with a random effect present", {
  ssn <- mf04p
  ssn$obs$block_group <- factor(rep(c("a", "b"), length.out = NROW(ssn$obs)))
  ssn$obs$block_partition <- factor(rep(c("in", "out"), length.out = NROW(ssn$obs)))
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    random = ~block_group, partition_factor = ~block_partition
  )
  expect_true(all(is.finite(coef(fit))))
  cm <- covmatrix(fit)
  groups <- ssn$obs$block_partition[fit$observed_index]
  expect_true(all(cm[outer(groups, groups, `!=`)] == 0))
  tg <- Torgegram(Summer_mn ~ ELEV_DEM, ssn, partition_factor = ~block_partition)
  expect_s3_class(tg, "Torgegram")
})

# Ordered partition factors
test_that("ordered partition factors are accepted and behave identically to an equivalent unordered factor", {
  ssn <- mf04p
  ssn$obs$level_group <- factor(rep(c("lo", "mid", "hi"), length.out = NROW(ssn$obs)),
    levels = c("lo", "mid", "hi"), ordered = TRUE
  )
  ssn$preds$CapeHorn$level_group <- factor(rep(c("lo", "mid", "hi"), length.out = NROW(ssn$preds$CapeHorn)),
    levels = c("lo", "mid", "hi"), ordered = TRUE
  )
  ssn_unordered <- ssn
  ssn_unordered$obs$level_group <- factor(as.character(ssn$obs$level_group), levels = c("lo", "mid", "hi"))
  ssn_unordered$preds$CapeHorn$level_group <- factor(
    as.character(ssn$preds$CapeHorn$level_group),
    levels = c("lo", "mid", "hi")
  )

  fit_ordered <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    partition_factor = ~level_group
  )
  fit_unordered <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn_unordered,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    partition_factor = ~level_group
  )
  expect_equal(coef(fit_ordered), coef(fit_unordered))
  expect_equal(as.matrix(covmatrix(fit_ordered)), as.matrix(covmatrix(fit_unordered)))

  pred_ordered <- predict(fit_ordered, "CapeHorn", se.fit = TRUE)
  pred_unordered <- predict(fit_unordered, "CapeHorn", se.fit = TRUE)
  expect_equal(pred_ordered, pred_unordered)

  # GLM exact and grouped/local paths
  fit_glm_ordered <- ssn_glm(Summer_mn ~ ELEV_DEM, ssn,
    family = "Gamma", tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    partition_factor = ~level_group
  )
  fit_glm_unordered <- ssn_glm(Summer_mn ~ ELEV_DEM, ssn_unordered,
    family = "Gamma", tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    partition_factor = ~level_group
  )
  expect_equal(coef(fit_glm_ordered), coef(fit_glm_unordered))

  ssn_create_bigdist(ssn, predpts = "CapeHorn", overwrite = TRUE, no_cores = 1, verbose = FALSE)
  ssn_create_bigdist(ssn_unordered, predpts = "CapeHorn", overwrite = TRUE, no_cores = 1, verbose = FALSE)
  fit_local_ordered <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    partition_factor = ~level_group, local = list(size = 30)
  )
  fit_local_unordered <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn_unordered,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    partition_factor = ~level_group, local = list(size = 30)
  )
  expect_equal(coef(fit_local_ordered), coef(fit_local_unordered))
})

test_that("ordered partition factors preserve unused and prediction-only levels", {
  ssn <- mf04p
  ssn$obs$level_group <- factor(rep(c("lo", "hi"), length.out = NROW(ssn$obs)),
    levels = c("lo", "mid", "hi"), ordered = TRUE
  )
  ssn$preds$CapeHorn$level_group <- factor(rep(c("lo", "mid", "hi"), length.out = NROW(ssn$preds$CapeHorn)),
    levels = c("lo", "mid", "hi"), ordered = TRUE
  )
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, ssn,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    partition_factor = ~level_group
  )
  expect_true(all(is.finite(coef(fit))))
  pred <- predict(fit, "CapeHorn", se.fit = TRUE)
  expect_true(all(is.finite(pred$fit)))
})

test_that("non-categorical and multi-variable partition factors are still rejected when ordered", {
  d <- data.frame(y = c(1, 2, 3))
  expect_error(SSN2:::coerce_partition_factor(~y, d), "categorical or factor")

  d2 <- data.frame(
    g = factor(c("a", "b"), ordered = TRUE),
    h = factor(c("x", "y"), ordered = TRUE)
  )
  expect_error(SSN2:::coerce_partition_factor(~ g + h, d2), "Only one variable")
})
