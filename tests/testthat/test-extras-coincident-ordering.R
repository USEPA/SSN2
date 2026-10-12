skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

conditional_ordering_data <- function() {
  sf::st_as_sf(data.frame(
    pid = c(3, 1, 1, 2), netID = c(1, 2, 1, 2),
    x = c(0, 0, 1, 2), y = c(0, 0, 1, 0)
  ), coords = c("x", "y"))
}

conditional_ordering_rows <- function(data) {
  data.frame(pid = data$pid, NetworkID = data$netID, row = seq_len(NROW(data)))
}

test_that("pid ordering keeps coincident rows and needs no optional package", {
  data <- conditional_ordering_data()
  rows <- conditional_ordering_rows(data)
  fun <- SSN2:::get_decorrelate_order
  environment(fun) <- new.env(parent = environment(fun))
  environment(fun)$requireNamespace <- function(...) stop("Unexpected dependency check")
  expect_identical(SSN2:::get_decorrelate_ordering(NULL), "pid")
  expect_identical(SSN2:::get_local_vecchia_settings(list())$ordering, "pid")
  expect_identical(fun(rows), c(3L, 2L, 4L, 1L))
  expect_identical(fun(rows, "pid"), fun(rows))
  expect_identical(fun(rows, "none"), seq_len(NROW(rows)))
  expect_error(fun(rows, "maxmin", data), "cannot be used with coincident locations", fixed = TRUE)
  expect_error(fun(rows, "invalid", data), "ordering must be")
  for (ordering in c("pid", "none", "random", "maxmin", "grts")) {
    expect_identical(fun(rows[FALSE, ], ordering), integer())
    expect_identical(fun(rows[1, , drop = FALSE], ordering), 1L)
  }
  set.seed(2)
  perm <- fun(rows, "random")
  expect_identical(sort(perm), seq_len(NROW(rows)))
  expect_identical(rows[perm, ][order(perm), ], rows)
})

test_that("ordering rejects illegal permutations", {
  rows <- conditional_ordering_rows(conditional_ordering_data())
  fun <- SSN2:::get_decorrelate_order
  environment(fun) <- new.env(parent = environment(fun))
  for (invalid in list(c(1L, 1L, 3L, 4L), 1:3, c(1, 2, 3, NA),
                       c(0, 1, 2, 3), c(1, 2, 3, Inf), c(1, 2, 3, 4.5),
                       letters[1:4])) {
    environment(fun)$order <- function(...) invalid
    expect_error(fun(rows, "pid"),
                 "Ordering did not properly rearrange the row indices.", fixed = TRUE)
  }
})

test_that("only explicitly requested spatial orderings require dependencies", {
  data <- conditional_ordering_data()[c(1, 3, 4), ]
  rows <- conditional_ordering_rows(data)
  fun <- SSN2:::get_decorrelate_order
  environment(fun) <- new.env(parent = environment(fun))
  environment(fun)$requireNamespace <- function(...) FALSE
  for (ordering in c("maxmin", "middleout", "outsidein", "coordinate")) {
    expect_error(fun(rows, ordering, data), "Install the GPvecchia package", fixed = TRUE)
  }
  expect_error(fun(rows, "grts", data), "Install the spsurvey package", fixed = TRUE)
})

test_that("explicit spatial orderings return row permutations", {
  set.seed(2)
  skip_if_not_installed("GPvecchia")
  data <- conditional_ordering_data()[c(1, 3, 4), ]
  rows <- conditional_ordering_rows(data)
  for (ordering in c("maxmin", "middleout", "outsidein", "coordinate")) {
    perm <- SSN2:::get_decorrelate_order(rows, ordering, data)
    expect_equal(sort(perm), seq_len(NROW(rows)))
    expect_identical(rows[perm, ][order(perm), ], rows)
  }
})

test_that("GRTS retains distinct records at coincident coordinates", {
  set.seed(2)
  skip_if_not_installed("spsurvey")
  data <- conditional_ordering_data()
  rows <- conditional_ordering_rows(data)
  perm <- SSN2:::get_decorrelate_order(rows, "grts", data)
  expect_equal(sort(perm), seq_len(NROW(rows)))
  expect_identical(rows[perm, ][order(perm), ], rows)
})

test_that("simulation and decorrelation callers share pid and duplicate safeguards", {
  set.seed(2)
  network <- mf04p
  network$obs <- network$obs[1:8, ]
  network$preds <- list(pred1km = network$preds$pred1km[1:4, ])
  network$obs$x <- seq(-1, 1, length.out = 8)
  network$obs$y <- c(1, 2, 3, 2, 4, 3, 2, 5)
  network$preds$pred1km$x <- c(-.5, .2, .8, 0)
  sf::st_geometry(network$obs)[2] <- sf::st_geometry(network$obs)[1]
  sf::st_geometry(network$preds$pred1km)[2] <- sf::st_geometry(network$preds$pred1km)[1]
  params <- list(tailup_params = tailup_params("none"),
                 taildown_params = taildown_params("none"),
                 euclid_params = euclid_params("exponential", de = .4, range = 1000),
                 nugget_params = nugget_params("nugget", nugget = .2))
  formula <- y ~ x
  transform <- function(ordering = NULL) {
    do.call(ssn_decorrelate_data, c(list(formula = formula, ssn.object = network,
      ordering = ordering, local = list(size = 3)), params))
  }
  expect_identical(transform(), transform("pid"))
  expect_identical(transform()$ordering, "pid")
  expect_error(transform("maxmin"), "coincident locations")

  simulate <- function(ordering = NULL) {
    do.call(ssn_rnorm, c(list(ssn.object = network, samples = 3,
      local = list(approximation = "vecchia", size = 3, ordering = ordering)), params))
  }
  set.seed(2)
  omitted <- simulate()
  set.seed(2)
  expect_identical(omitted, simulate("pid"))
  expect_true(all(is.finite(omitted)))
  expect_error(simulate("maxmin"), "coincident locations")

  for (glm in c(FALSE, TRUE)) {
    args <- list(formula = y ~ x, ssn.object = network, tailup_type = "none",
      euclid_type = "exponential",
      euclid_initial = euclid_initial("exponential", de = .4, range = 1000, known = "given"),
      nugget_initial = nugget_initial("nugget", nugget = .2, known = "given"))
    if (glm) args$family <- "poisson"
    fit <- do.call(if (glm) ssn_glm else ssn_lm, args)
    set.seed(2)
    omitted <- conditional(fit, "pred1km", samples = 3,
      local = list(approximation = "vecchia", size = 3))
    set.seed(2)
    explicit <- conditional(fit, "pred1km", samples = 3,
      local = list(approximation = "vecchia", size = 3, ordering = "pid"))
    expect_identical(omitted, explicit)
    expect_equal(dim(omitted), c(4, 3))
    expect_true(all(is.finite(omitted)))
    expect_error(conditional(fit, "pred1km", samples = 3,
      local = list(approximation = "vecchia", size = 3, ordering = "maxmin")),
      "coincident locations")
    if (!glm) {
      expect_error(predict(fit, "pred1km", block = TRUE,
        local = list(method = "covariance", method_new = "subset", size_new = 2,
                     ordering = "maxmin")),
        "coincident locations")
    }
  }
})
