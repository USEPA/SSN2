test_that("conditional simulation and K-fold validation dispatch on a stream model", {
  network <- mf04p
  fit <- ssn_lm(Summer_mn ~ ELEV_DEM, network, additive = "afvArea",
    tailup_initial = tailup_initial("exponential", 1, 10000, known = "given"),
    nugget_initial = nugget_initial("nugget", 0.5, known = "given"),
    ddf = "asymptotic")
  set.seed(2)
  draws <- conditional(fit, "CapeHorn", samples = 2)
  expect_equal(dim(draws), c(nrow(network$preds$CapeHorn), 2L))
  expect_true(all(is.finite(draws)))
  cv <- kcv(fit, folds_index = rep(1:3, length.out = nobs(fit)), cv_predict = TRUE)
  expect_length(cv$cv_predict, nobs(fit))
  expect_equal(cv$stats$MSPE, mean((network$obs$Summer_mn - cv$cv_predict)^2))
})

test_that("decorrelation grids and transformations support basic use", {
  grid <- ssn_decorrelate_grid(Summer_mn ~ ELEV_DEM, mf04p,
    tailup_type = "exponential", additive = "afvArea", dense_grid = FALSE)
  expect_s3_class(grid, "ssn_decorrelate_grid")
  expect_true(any(tidy(grid)$tailup_type == "no transformation"))
  transformed <- ssn_decorrelate_data(Summer_mn ~ ELEV_DEM, mf04p,
    nugget_params = nugget_params("nugget", 1),
    ordering = "none", local = FALSE)
  expect_equal(as.numeric(transformed$tX), as.numeric(transformed$X))
  expect_equal(transformed$ty, transformed$y)
  prediction <- ssn_decorrelate_newdata(transformed, "CapeHorn")
  # tX_newdata is a computed transform, not a raw model matrix, so (unlike
  # X_newdata) it carries no "assign" attribute
  expect_equal(unname(prediction$tX_newdata), unname(prediction$X_newdata), ignore_attr = "assign")
  values <- seq_len(nrow(prediction$X_newdata))
  expect_equal(as.numeric(ssn_recorrelate_newdata(prediction, values)), as.numeric(values))
})
