test_that("stream wave constructors retain two-parameter specifications", {
  for (component in c("tailup", "taildown")) {
    params <- get(paste0(component, "_params"))("wave", de = 2, range = 10)
    initial <- get(paste0(component, "_initial"))("wave", de = 2, range = 10, known = "given")
    expect_s3_class(params, paste0(component, "_wave"))
    expect_identical(names(params), c("de", "range"))
    expect_equal(initial$initial, c(de = 2, range = 10))
    expect_true(all(initial$is_known))
    expect_equal(get_cov_initial_range(initial, c(1, 20)), c(1, 20))
    expect_warning(value <- get_effective_range(params),
      '"wave" effective range is not well defined. A potential candidate value is used.', fixed = TRUE)
    expect_equal(value, 200)
  }
})

test_that("stream wave limits and negative values differ from Euclidean wave", {
  skip_if_not_installed("pracma")
  h <- c(0, 0.5, 1, 2, 3.47)
  expected <- c(1, 0.804075, 0.557199, 0.060819, -0.262327)
  expect_equal(stream_wave_connected(h, 1), expected, tolerance = 1e-6)
  expect_equal(stream_wave_unconnected(rep(0, length(h)), h, 1),
               stream_wave_connected(h, 1), tolerance = 1e-12)
  expect_equal(stream_wave_unconnected(1, 1, 1), 0.428736, tolerance = 1e-6)
  expect_equal((1 - stream_wave_connected(1e-7, 1)) / 1e-7, 1 / pi, tolerance = 1e-6)
  expect_equal(stream_wave_connected(numeric(), 1), numeric())
  expect_error(stream_wave_connected(-1, 1), "nonnegative")
  expect_error(stream_wave_connected(1, 0), "finite and positive")
  expect_error(stream_wave_unconnected(2, 1, 1), "a <= b", fixed = TRUE)
  expect_error(stream_wave_unconnected(1, c(1, 2), 1), "equal lengths")
  euclidean <- cov_matrix(euclid_params("wave", de = 1, range = 1),
                         list(euclid_mat = h), anisotropy = FALSE)
  expect_equal(as.numeric(euclidean), c(1, sin(h[-1]) / h[-1]))
  expect_gt(abs(euclidean[2] - expected[2]), 0.1)
})

test_that("stream wave dispatch preserves masks, dimensions, and names", {
  skip_if_not_installed("pracma")
  a <- matrix(c(0, 0, 1, NA, 0, 1), 2,
              dimnames = list(c("7", "2"), c("11", "4", "19")))
  b <- matrix(c(0, 3.47, 1, NA, 0.5, 2), 2, dimnames = dimnames(a))
  mask <- matrix(c(1, 1, 1, 0, 1, 1), 2, dimnames = dimnames(a))
  w <- mask * (a == 0)
  w[mask == 0] <- 0
  w[2, 1] <- 0.5
  for (sparse in c(FALSE, TRUE)) {
    observed <- list(hydro_mat = a + b, a_mat = a, b_mat = b, mask_mat = mask, w_mat = w)
    if (sparse) observed <- lapply(observed, Matrix::Matrix, sparse = TRUE)
    predicted <- observed
    names(predicted) <- sub("_mat$", "_pred_mat", names(predicted))
    for (component in c("tailup", "taildown")) {
      params <- get(paste0(component, "_params"))("wave", de = 2, range = 1)
      square <- as.matrix(cov_matrix(params, observed))
      cross <- as.matrix(cov_vector(params, predicted))
      expect_equal(cross, square)
      expect_identical(dimnames(square), dimnames(a))
      expect_identical(dim(square), c(2L, 3L))
      expect_true(all(is.finite(square)))
      expect_equal(square[2, 2], 0)
      expect_equal(square[1, 1], 2)
      expect_lt(square[2, 1], 0)
      expect_equal(square[1, 2], if (component == "tailup") 0 else
        2 * stream_wave_unconnected(1, 1, 1))
      for (rows in list(1L, integer())) {
        block <- lapply(predicted, function(x) x[rows, 2:3, drop = FALSE])
        expect_equal(as.matrix(cov_vector(params, block)), square[rows, 2:3, drop = FALSE])
      }
    }
  }
})

test_that("only stream wave evaluation requires pracma", {
  dependency <- stream_wave_dependency
  environment(dependency) <- new.env(parent = environment(dependency))
  environment(dependency)$requireNamespace <- function(...) FALSE
  local_mocked_bindings(stream_wave_dependency = dependency)
  params <- tailup_params("wave", de = 1, range = 2)
  expect_error(cov_matrix(params, list(hydro_mat = 0, mask_mat = 1, w_mat = 1)),
               "install.packages('pracma')", fixed = TRUE)
  params <- taildown_params("wave", de = 1, range = 2)
  expect_error(cov_vector(params, list(hydro_pred_mat = 0, a_pred_mat = 0,
    b_pred_mat = 0, mask_pred_mat = 1)), "install.packages('pracma')", fixed = TRUE)
  expect_equal(cov_matrix(tailup_params("exponential", de = 1, range = 2),
                         list(hydro_mat = 0, mask_mat = 1, w_mat = 1)), 1)
})
