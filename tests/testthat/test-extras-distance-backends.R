skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("binary-ID matching preserves native prefix and connectivity conventions", {
  ids <- c("1", "111", "101", "011", "11", "")
  expect_identical(get_binary_id_match(ids, "11"), c(-1L, -2L, 1L, 0L, -2L, 0L))
  expect_identical(get.rid.fc(ids, "11"), data.frame(
    fc = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
    binaryID = c("1", "11", "1", "", "11", "")))
  expect_identical(get_binary_id_match(character(), "1"), integer())
  expect_identical(get_binary_id_match(ids, ""), integer(length(ids)))

  prefix <- paste(rep("1", 4096), collapse = "")
  expect_identical(get_binary_id_match(c(prefix, paste0(prefix, "0"),
    paste0(prefix, "11")), paste0(prefix, "1")), c(-4096L, 4096L, -4097L))
})

test_that("binary-ID matching agrees with character-wise comparison", {
  ids <- c("", unlist(lapply(1:5, function(n) {
    apply(expand.grid(rep(list(c("0", "1")), n)), 1, paste0, collapse = "")
  }), use.names = FALSE))
  expected <- vapply(ids, function(reference) {
    vapply(ids, function(id) {
      size <- min(nchar(id), nchar(reference))
      a <- strsplit(id, "", fixed = TRUE)[[1]][seq_len(size)]
      b <- strsplit(reference, "", fixed = TRUE)[[1]][seq_len(size)]
      mismatch <- which(a != b)
      if (length(mismatch)) mismatch[1] - 1L else -size
    }, integer(1), USE.NAMES = FALSE)
  }, integer(length(ids)), USE.NAMES = FALSE)
  actual <- vapply(ids, function(reference) get_binary_id_match(ids, reference),
                   integer(length(ids)), USE.NAMES = FALSE)
  expect_identical(actual, expected)
})

test_that("binary-ID matching rejects malformed helper inputs", {
  expect_error(get_binary_id_match(1, "1"), "invalid arguments")
  expect_error(get_binary_id_match("1", c("1", "11")), "length one")
  expect_error(get_binary_id_match(NA_character_, "1"), "missing binary ID")
  expect_error(get_binary_id_match("1", NA_character_), "missing binary ID")
})

# distance blocks
ssn_create_bigdist(mf04p, predpts = "CapeHorn", overwrite = TRUE, no_cores = 1, verbose = FALSE)

test_that("cross-block distance reader honors arbitrary requested pid order", {
  initial <- get_initial_object("exponential", "exponential", "exponential", "nugget", NULL, NULL, NULL, NULL)
  dist0 <- get_dist_object(mf04p, initial, "afvArea", FALSE)
  ng <- ssn_get_netgeom(mf04p$obs)
  id1 <- order(ng$NetworkID, ng$pid)[seq_len(9)]
  id2 <- order(ng$NetworkID, ng$pid)[10:19]

  canonical <- get_distjunc_matlist_bigdata_cross(ng$NetworkID[id1], ng$pid[id1], ng$NetworkID[id2], ng$pid[id2], mf04p)
  expect_equal(canonical, as.matrix(dist0$distjunc_mat[id1, id2]), ignore_attr = TRUE)

  reversed <- get_distjunc_matlist_bigdata_cross(ng$NetworkID[rev(id1)], ng$pid[rev(id1)], ng$NetworkID[rev(id2)], ng$pid[rev(id2)], mf04p)
  expect_equal(reversed, as.matrix(dist0$distjunc_mat[rev(id1), rev(id2)]), ignore_attr = TRUE)

  # a genuinely arbitrary (non-monotonic) permutation, not just a full reversal
  set.seed(2)
  perm1 <- sample(id1)
  perm2 <- sample(id2)
  arbitrary <- get_distjunc_matlist_bigdata_cross(ng$NetworkID[perm1], ng$pid[perm1], ng$NetworkID[perm2], ng$pid[perm2], mf04p)
  expect_equal(arbitrary, as.matrix(dist0$distjunc_mat[perm1, perm2]), ignore_attr = TRUE)
})

test_that("cross-block distance reader errors on a pid absent from its claimed network's file", {
  ng <- ssn_get_netgeom(mf04p$obs)
  id1 <- order(ng$NetworkID, ng$pid)[seq_len(5)]
  bogus_pid <- c(ng$pid[id1], "999999")
  bogus_network <- c(ng$NetworkID[id1], ng$NetworkID[id1[1]])
  expect_error(
    get_distjunc_matlist_bigdata_cross(bogus_network, bogus_pid, ng$NetworkID[id1], ng$pid[id1], mf04p),
    "Unable to locate stored distance information"
  )
})

test_that("cross-block distance reader still returns legitimate cross-network zeros without erroring", {
  ng <- ssn_get_netgeom(mf04p$obs)
  id_net1 <- which(ng$NetworkID == 1)[seq_len(3)]
  id_net2 <- which(ng$NetworkID == 2)[seq_len(3)]
  cross <- get_distjunc_matlist_bigdata_cross(ng$NetworkID[id_net1], ng$pid[id_net1], ng$NetworkID[id_net2], ng$pid[id_net2], mf04p)
  expect_true(all(cross == 0))
})

test_that("prediction distance orientation is correct even when nobs == npreds", {
  # a separate on-disk copy, distinct from the shared mf04p fixture's
  # tempdir()/MiddleFork04.ssn -- this test subsets and overwrites distance
  # files in place, which would corrupt mf04p's own copy if reused directly
  iso_path <- file.path(tempdir(), "MiddleFork04_equalcount.ssn")
  if (!dir.exists(iso_path)) {
    dir.create(iso_path, recursive = TRUE, showWarnings = FALSE)
    src <- system.file("lsndata/MiddleFork04.ssn", package = "SSN2")
    file.copy(list.files(src, full.names = TRUE), iso_path, recursive = TRUE, overwrite = FALSE, copy.mode = FALSE)
  }
  square <- ssn_import(iso_path, predpts = "CapeHorn", overwrite = TRUE, verbose = FALSE)
  ng <- ssn_get_netgeom(square$obs)
  net2_ids <- which(ng$NetworkID == 2)
  square$obs <- square$obs[net2_ids[c(3, 7, 1, 12, 2, 10, 5, 11, 4, 8, 6, 9)], ]
  square$preds$CapeHorn <- square$preds$CapeHorn[12:1, ]
  ssn_create_distmat(square, predpts = "CapeHorn", among_predpts = TRUE, overwrite = TRUE)

  tu <- tailup_initial("exponential", de = 2, range = 10000, known = c("de", "range"))
  td <- taildown_initial("exponential", de = 1, range = 15000, known = c("de", "range"))
  ng_init <- nugget_initial("nugget", nugget = 0.5, known = "nugget")
  square_model <- ssn_lm(Summer_mn ~ ELEV_DEM, square,
    tailup_initial = tu, taildown_initial = td, nugget_initial = ng_init,
    additive = "afvArea"
  )

  initial <- get_initial_object("exponential", "exponential", "exponential", "nugget", NULL, NULL, NULL, NULL)
  square_dist <- get_dist_pred_object(square_model, "CapeHorn", initial)
  reference <- get_distjunc_pred_matlist(square_model$ssn.object, "CapeHorn", square_dist)
  raw_a <- as.matrix(Matrix::bdiag(reference$distjunca))
  raw_b <- as.matrix(Matrix::bdiag(reference$distjuncb))
  expected_a <- t(raw_a[square_dist$inv_dist_order, square_dist$inv_dist_order_pred, drop = FALSE])
  expected_b <- raw_b[square_dist$inv_dist_order_pred, square_dist$inv_dist_order, drop = FALSE]

  expect_equal(as.matrix(square_dist$distjunca_pred_mat), expected_a, ignore_attr = TRUE)
  # a real regression guard, not just a numerical-tolerance check: the old
  # code silently applied the obs-by-pred index order to this pred-by-obs
  # field whenever nobs == npreds, scrambling entries rather than erroring
  expect_equal(t(as.matrix(square_dist$distjuncb_pred_mat)), expected_b, ignore_attr = TRUE)
  expect_equal(as.matrix(square_dist$hydro_pred_mat), expected_a + expected_b, ignore_attr = TRUE)

  preds <- predict(square_model, "CapeHorn")
  expect_length(preds, 12)
  expect_true(all(is.finite(preds)))
})

# backend routing
iso_root <- file.path(tempdir(), "backend_routing_iso")
iso_path <- file.path(iso_root, "MiddleFork04.ssn")
dir.create(iso_path, recursive = TRUE, showWarnings = FALSE)
src <- system.file("lsndata/MiddleFork04.ssn", package = "SSN2")
file.copy(list.files(src, full.names = TRUE), iso_path, recursive = TRUE, overwrite = FALSE, copy.mode = FALSE)

br_ssn <- ssn_import(iso_path, predpts = c("pred1km", "CapeHorn"), overwrite = TRUE)
ssn_create_distmat(br_ssn, predpts = c("pred1km", "CapeHorn"), overwrite = TRUE, among_predpts = TRUE)
ssn_create_bigdist(br_ssn, predpts = "CapeHorn", overwrite = TRUE, among_predpts = TRUE, verbose = FALSE)
ssn_create_bigdist(br_ssn, predpts = "pred1km", overwrite = TRUE, among_predpts = TRUE, verbose = FALSE)

br_rdata_files <- function(dirs) {
  unlist(lapply(dirs, function(d) list.files(file.path(iso_path, "distance", d), pattern = "\\.RData$", full.names = TRUE)))
}

hide_rdata <- function(dirs) {
  files <- br_rdata_files(dirs)
  hidden <- paste0(files, ".hidden")
  file.rename(files, hidden)
  hidden
}

restore_rdata <- function(hidden) {
  file.rename(hidden, sub("\\.hidden$", "", hidden))
}

br_bmat_files <- function(dirs) {
  unlist(lapply(dirs, function(d) list.files(file.path(iso_path, "distance", d), pattern = "\\.bmat$", full.names = TRUE)))
}

hide_bmat <- function(dirs) {
  files <- br_bmat_files(dirs)
  hidden <- paste0(files, ".hidden")
  file.rename(files, hidden)
  hidden
}

restore_bmat <- function(hidden) {
  file.rename(hidden, sub("\\.hidden$", "", hidden))
}

test_that("all four covmatrix() modes match between dense and filematrix-only backends, including named transpose equality", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  dense_obsobs <- covmatrix(fit)
  dense_obspred <- covmatrix(fit, "CapeHorn", cov_type = "obs.pred")
  dense_predobs <- covmatrix(fit, "CapeHorn", cov_type = "pred.obs")
  dense_predpred <- covmatrix(fit, "CapeHorn", cov_type = "pred.pred")
  expect_equal(dense_obspred, t(dense_predobs))

  hidden <- hide_rdata(c("obs", "CapeHorn"))
  on.exit(restore_rdata(hidden), add = TRUE)

  expect_equal(covmatrix(fit), dense_obsobs)
  bigdata_obspred <- covmatrix(fit, "CapeHorn", cov_type = "obs.pred")
  bigdata_predobs <- covmatrix(fit, "CapeHorn", cov_type = "pred.obs")
  expect_equal(bigdata_obspred, dense_obspred)
  expect_equal(bigdata_predobs, dense_predobs)
  expect_equal(covmatrix(fit, "CapeHorn", cov_type = "pred.pred"), dense_predpred)
  expect_equal(bigdata_obspred, t(bigdata_predobs))
})

test_that("ssn_lm(local = TRUE) model fitting with the .bmat-only backend matches the dense reference", {
  set.seed(2)
  dense_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000),
    local = list(size = 20)
  )

  hidden <- hide_rdata("obs")
  on.exit(restore_rdata(hidden), add = TRUE)

  set.seed(2)
  bigdata_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000),
    local = list(size = 20)
  )
  expect_equal(coef(bigdata_fit), coef(dense_fit), tolerance = 1e-8)
  expect_equal(
    unname(unlist(coef(bigdata_fit, type = "tailup"))),
    unname(unlist(coef(dense_fit, type = "tailup"))),
    tolerance = 1e-6
  )
  expect_equal(logLik(bigdata_fit), logLik(dense_fit), tolerance = 1e-6)
})

test_that("ssn_lm(local = TRUE) model fitting succeeds with only .RData distance matrices (no .bmat at all)", {
  # regression test: get_dist_object_bigdata() (the distance-object builder
  # local model fitting uses) previously hardcoded backend = "bigdata" with
  # no fallback, so local model fitting hard-errored on any dataset that
  # never ran ssn_create_bigdist() -- even though ssn_create_distmat()'s
  # .RData files were sufficient
  set.seed(2)
  dense_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000),
    local = list(size = 20)
  )

  hidden <- hide_bmat("obs")
  on.exit(restore_bmat(hidden), add = TRUE)

  set.seed(2)
  rdata_only_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000),
    local = list(size = 20)
  )
  expect_equal(coef(rdata_only_fit), coef(dense_fit), tolerance = 1e-8)
})

test_that("predict(local = 'covariance') with the .bmat-only backend matches the dense reference", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  size <- max(1, fit$n - 10)
  dense_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)

  hidden <- hide_rdata("CapeHorn")
  on.exit(restore_rdata(hidden), add = TRUE)

  bigdata_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)
  expect_equal(bigdata_pred$fit, dense_pred$fit, tolerance = 1e-8)
  expect_equal(bigdata_pred$se.fit, dense_pred$se.fit, tolerance = 1e-8)
})

test_that("predict_block(local = 'covariance') with the .bmat-only backend matches the dense reference", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  size <- max(1, fit$n - 10)
  dense_pred <- predict(fit, "CapeHorn", block = TRUE, local = list(method = "covariance", size = size), se.fit = TRUE)

  hidden <- hide_rdata("CapeHorn")
  on.exit(restore_rdata(hidden), add = TRUE)

  bigdata_pred <- predict(fit, "CapeHorn", block = TRUE, local = list(method = "covariance", size = size), se.fit = TRUE)
  expect_equal(bigdata_pred$fit, dense_pred$fit, tolerance = 1e-8)
  expect_equal(bigdata_pred$se.fit, dense_pred$se.fit, tolerance = 1e-8)
})

test_that("ssn_glm predict(local = 'covariance') with the .bmat-only backend matches the dense reference", {
  s <- br_ssn
  s$obs$y_pois <- round(s$obs$Summer_mn)
  fit <- ssn_glm(y_pois ~ ELEV_DEM, s,
    family = "poisson", tailup_type = "exponential", taildown_type = "exponential",
    euclid_type = "exponential", nugget_type = "nugget", additive = "afvArea"
  )
  size <- max(1, fit$n - 10)
  dense_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)

  hidden <- hide_rdata("CapeHorn")
  on.exit(restore_rdata(hidden), add = TRUE)

  bigdata_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)
  expect_equal(bigdata_pred$fit, dense_pred$fit, tolerance = 1e-6)
  expect_equal(bigdata_pred$se.fit, dense_pred$se.fit, tolerance = 1e-6)
})

test_that("predicting originally-missing responses ('.missing') with the .bmat-only backend matches the dense reference", {
  s <- br_ssn
  s$obs$Summer_mn[c(2, 9, 17)] <- NA
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  size <- max(1, fit$n - 5)
  dense_pred <- predict(fit, ".missing", local = list(method = "covariance", size = size), se.fit = TRUE)

  hidden <- hide_rdata("obs")
  on.exit(restore_rdata(hidden), add = TRUE)

  bigdata_pred <- predict(fit, ".missing", local = list(method = "covariance", size = size), se.fit = TRUE)
  expect_equal(bigdata_pred$fit, dense_pred$fit, tolerance = 1e-8)
  expect_equal(bigdata_pred$se.fit, dense_pred$se.fit, tolerance = 1e-8)
})

test_that("one-site prediction with the .bmat-only backend matches the dense reference", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "none", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  fit$ssn.object$preds$CapeHorn <- fit$ssn.object$preds$CapeHorn[1, , drop = FALSE]
  size <- max(1, fit$n - 10)

  dense_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)

  hidden <- hide_rdata("CapeHorn")
  on.exit(restore_rdata(hidden), add = TRUE)

  bigdata_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)
  expect_equal(bigdata_pred$fit, dense_pred$fit, tolerance = 1e-8)
  expect_equal(bigdata_pred$se.fit, dense_pred$se.fit, tolerance = 1e-8)
})

test_that("a partition factor combined with local prediction matches the dense reference on the .bmat-only backend", {
  s <- br_ssn
  s$obs$part_group <- factor(rep(1:3, length.out = nrow(s$obs)))
  s$preds$CapeHorn$part_group <- factor(rep(1:3, length.out = nrow(s$preds$CapeHorn)), levels = levels(s$obs$part_group))
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    partition_factor = ~part_group
  )
  size <- max(1, fit$n - 10)
  dense_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)

  hidden <- hide_rdata("CapeHorn")
  on.exit(restore_rdata(hidden), add = TRUE)

  bigdata_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)
  expect_equal(bigdata_pred$fit, dense_pred$fit, tolerance = 1e-8)
  expect_equal(bigdata_pred$se.fit, dense_pred$se.fit, tolerance = 1e-8)
})

test_that("a random intercept combined with local prediction matches the dense reference on the .bmat-only backend", {
  s <- br_ssn
  s$obs$rand_group <- factor(rep(1:3, length.out = nrow(s$obs)))
  s$preds$CapeHorn$rand_group <- factor(rep(1:3, length.out = nrow(s$preds$CapeHorn)), levels = levels(s$obs$rand_group))
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, s,
    tailup_type = "exponential", nugget_type = "nugget", additive = "afvArea",
    random = ~rand_group
  )
  size <- max(1, fit$n - 10)
  dense_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)

  hidden <- hide_rdata("CapeHorn")
  on.exit(restore_rdata(hidden), add = TRUE)

  bigdata_pred <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)
  expect_equal(bigdata_pred$fit, dense_pred$fit, tolerance = 1e-8)
  expect_equal(bigdata_pred$se.fit, dense_pred$se.fit, tolerance = 1e-8)
})

test_that("an informative error is raised when neither .RData nor .bmat distance files exist for a prediction dataset (item 2)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  all_files <- list.files(file.path(iso_path, "distance", "CapeHorn"), pattern = "\\.(RData|bmat)$", full.names = TRUE)
  hidden <- paste0(all_files, ".hidden")
  file.rename(all_files, hidden)
  on.exit(file.rename(hidden, all_files), add = TRUE)

  expect_error(covmatrix(fit, "CapeHorn"), "CapeHorn")
})

test_that("an informative error is raised when among-predpts distances are missing for pred.pred (item 2)", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  square_files <- list.files(file.path(iso_path, "distance", "CapeHorn"),
    pattern = "^dist\\.net[0-9]+\\.(RData|bmat)$", full.names = TRUE
  )
  hidden <- paste0(square_files, ".hidden")
  file.rename(square_files, hidden)
  on.exit(file.rename(hidden, square_files), add = TRUE)

  expect_error(covmatrix(fit, "CapeHorn", cov_type = "pred.pred"), "CapeHorn")
})

test_that("the bounded-block accumulator processes multiple network chunks and matches colMeans(covmatrix()) exactly (item 3)", {
  # pred1km (unlike CapeHorn, which lands entirely on one network) spans both
  # networks in this fixture, so this is the case that genuinely exercises
  # more than one chunk in the accumulator.
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  netgeom_pred <- ssn_get_netgeom(fit$ssn.object$preds$pred1km, reformat = TRUE)
  n_networks <- length(unique(netgeom_pred$NetworkID))
  expect_gt(n_networks, 1) # confirms this fixture genuinely exercises >1 chunk

  dense_means <- colMeans(covmatrix(fit, "pred1km"))

  hidden <- hide_rdata("pred1km")
  on.exit(restore_rdata(hidden), add = TRUE)

  initial_object_val <- get_initial_object_from_coef(fit)
  bigdata_means <- get_local_cov_means_bigdata(fit, "pred1km", initial_object_val, fit$coefficients$params_object)
  expect_equal(bigdata_means, dense_means, tolerance = 1e-9)

  # every network chunk is strictly smaller than the full prediction count,
  # confirming no single chunk read reaches the full m x n size
  chunk_sizes <- table(netgeom_pred$NetworkID)
  expect_true(all(chunk_sizes < nrow(fit$ssn.object$preds$pred1km)))
})

test_that("existing dense (.RData) behavior is unchanged when both backends are available", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  size <- max(1, fit$n - 10)
  pred_val <- predict(fit, "CapeHorn", local = list(method = "covariance", size = size), se.fit = TRUE)
  expect_true(all(is.finite(pred_val$fit)))
  expect_true(all(is.finite(pred_val$se.fit)))
  expect_equal(select_pred_dist_backend(fit$ssn.object, "CapeHorn", FALSE, FALSE), "dense")
  expect_equal(select_square_dist_backend(fit$ssn.object, "obs", FALSE, FALSE), "dense")
})

test_that("select_*_dist_backend() honor prefer = \"bigdata\" when both backends are available", {
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  expect_equal(select_square_dist_backend(fit$ssn.object, "obs", FALSE, FALSE, prefer = "bigdata"), "bigdata")
  expect_equal(select_pred_dist_backend(fit$ssn.object, "CapeHorn", FALSE, FALSE, prefer = "bigdata"), "bigdata")
  # prefer = "dense" (the default, used by the exact/non-local paths) is unaffected
  expect_equal(select_square_dist_backend(fit$ssn.object, "obs", FALSE, FALSE), "dense")
  expect_equal(select_pred_dist_backend(fit$ssn.object, "CapeHorn", FALSE, FALSE), "dense")
})

test_that("block prediction's own backend resolver now prefers .bmat, matching decorrelation's", {
  # Indexed reads avoid whole-file deserialization for each requested chunk.
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  expect_equal(get_block_pred_backend(fit, "CapeHorn"), "bigdata")
  expect_equal(get_block_pred_backend(fit, "CapeHorn", square = TRUE), "bigdata")
  expect_equal(get_block_pred_backend(fit, "CapeHorn", prefer = "dense"), "dense")
  expect_equal(get_decorrelate_observed_backend(fit), "bigdata")
})

test_that("ssn_lm(local = TRUE) reads the big-data backend (not .RData) when both exist", {
  # corrupting (not hiding) the .RData bytes means a fit that still reads
  # .RData would fail with a deserialization error; a successful,
  # numerically-matching fit is direct proof the .bmat backend was read
  # instead, confirming local model fitting now prefers .bmat over .RData
  # whenever both are present (previously the reverse)
  set.seed(2)
  dense_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000),
    local = list(size = 20)
  )

  rdata_files <- br_rdata_files("obs")
  original_bytes <- lapply(rdata_files, function(f) readBin(f, what = "raw", n = file.info(f)$size))
  on.exit(mapply(function(bytes, f) writeBin(bytes, f), original_bytes, rdata_files), add = TRUE)
  invisible(lapply(rdata_files, function(f) writeBin(as.raw(0), f)))

  set.seed(2)
  bigdata_fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000),
    local = list(size = 20)
  )
  expect_equal(coef(bigdata_fit), coef(dense_fit), tolerance = 1e-8)
})

test_that("predict()'s point-level path prefers .bmat and only chunks then (T25-02/T25-03)", {
  # Dense point prediction stays unchunked to avoid repeated whole-file reads.
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  bigdata_ref <- get_point_pred_cov_vector_list(fit, "CapeHorn", chunk_size = 1000)
  bigdata_tiny <- get_point_pred_cov_vector_list(fit, "CapeHorn", chunk_size = 2)
  expect_equal(bigdata_tiny, bigdata_ref)

  hidden <- hide_bmat(c("obs", "CapeHorn"))
  on.exit(restore_bmat(hidden), add = TRUE)

  expect_equal(get_block_pred_backend(fit, "CapeHorn", prefer = "bigdata"), "dense")
  dense_ref <- get_point_pred_cov_vector_list(fit, "CapeHorn", chunk_size = 1000)
  dense_tiny <- get_point_pred_cov_vector_list(fit, "CapeHorn", chunk_size = 2)
  expect_equal(dense_tiny, dense_ref)
  expect_equal(dense_ref, bigdata_ref)
})

test_that("predict()'s point-level path actually reads every chunk with the backend it selected, not just the outer choice (backend propagation)", {
  # regression coverage for a defect where get_point_pred_cov_vector_list()
  # resolved prefer = "bigdata" to decide *whether* to chunk, but each
  # chunk's own get_block_obs_covariance() call re-resolved the backend with
  # get_block_pred_backend()'s un-preferenced default (prefer = "dense"):
  # with both formats present, that means chunking was correctly activated,
  # but every chunk actually deserialized the dense .RData file anyway --
  # defeating the entire point of preferring .bmat. Asserting only on
  # get_block_pred_backend()'s outer return value (as the test above does)
  # cannot catch this: it needs to inspect what backend each chunk read was
  # actually issued with.
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, br_ssn,
    tailup_type = "exponential", taildown_type = "exponential", euclid_type = "exponential",
    nugget_type = "nugget", additive = "afvArea", control = list(maxit = 2000)
  )
  reference <- get_point_pred_cov_vector_list(fit, "CapeHorn", chunk_size = 1000)

  real_reader <- get_dist_pred_object
  backends_used <- character()
  reader <- function(object, newdata_name, initial_object, backend = NULL) {
    backends_used[[length(backends_used) + 1L]] <<- backend
    real_reader(object, newdata_name, initial_object, backend)
  }
  testthat::local_mocked_bindings(get_dist_pred_object = reader, .package = "SSN2")

  instrumented <- get_point_pred_cov_vector_list(fit, "CapeHorn", chunk_size = 2)
  expect_true(length(backends_used) > 1) # confirms chunking actually occurred
  expect_true(all(backends_used == "bigdata"))
  expect_equal(instrumented, reference)
})
