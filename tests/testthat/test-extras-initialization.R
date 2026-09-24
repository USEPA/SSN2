skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

# cov initial grid
grid_test_inputs <- function(types = c("exponential", "exponential", "exponential", "nugget"),
                             anisotropy = FALSE, random = character(), dispersion = NULL) {
  data <- list(s2 = 10 / 1.2, tail_max = 100, euclid_max = 40,
               anisotropy = anisotropy, randcov_names = random)
  initial <- list(
    tailup_initial = tailup_initial_NA(tailup_initial(types[1])),
    taildown_initial = taildown_initial_NA(taildown_initial(types[2])),
    euclid_initial = euclid_initial_NA(euclid_initial(types[3]), data),
    nugget_initial = nugget_initial_NA(nugget_initial(types[4])),
    randcov_initial = if (length(random)) list(initial = setNames(rep(NA_real_, length(random)), random)) else NULL
  )
  if (!is.null(dispersion)) initial$dispersion_initial <- list(initial = c(dispersion = dispersion))
  list(initial = initial, data = data)
}

test_that("the medium grid covers the specified 48 allocation/range combinations", {
  fixture <- grid_test_inputs()
  grid <- build_cov_initial_grid(fixture$initial, fixture$data)
  vars <- c("tailup_de", "taildown_de", "euclid_de", "nugget")
  expect_equal(nrow(grid), 48)
  expect_equal(nrow(unique(grid)), 48)
  expect_equal(rowSums(grid[, vars]), rep(10, 48))
  expect_equal(sort(unique(grid$tailup_range)), c(25, 75) / 3)
  expect_equal(sort(unique(grid$taildown_range)), c(25, 75) / 3)
  expect_equal(sort(unique(grid$euclid_range)), c(10, 30) / 3)
  proportions <- as.matrix(grid[, vars]) / 10
  count <- function(w) sum(apply(abs(sweep(proportions, 2, w)), 1, max) < 1e-12)
  expect_equal(count(rep(0.25, 4)), 8)
  expect_equal(count(c(0.3, 0.3, 0.3, 0.1)), 8)
  expect_equal(count(c(0.45, 0.45, 0.05, 0.05)), 4)
  expect_equal(count(c(0.45, 0.05, 0.45, 0.05)), 4)
  expect_equal(count(c(0.05, 0.45, 0.45, 0.05)), 4)
  expect_equal(count(c(0.9, rep(0.1 / 3, 3))), 2)
  pair <- which(abs(proportions[, 1] - 0.45) < 1e-12 & abs(proportions[, 2] - 0.45) < 1e-12)
  expect_equal(unname(as.matrix(grid[pair, c("tailup_range", "taildown_range", "euclid_range")])),
               rbind(c(25, 25, 10), c(75, 75, 30), c(75, 75, 10), c(25, 25, 30)) / 3)
})

test_that("all active-component subsets have valid normalized starts", {
  combinations <- expand.grid(rep(list(c(FALSE, TRUE)), 4))
  names <- c("tailup_de", "taildown_de", "euclid_de", "nugget")
  for (i in seq_len(nrow(combinations))) {
    active <- as.logical(combinations[i, ])
    fixture <- grid_test_inputs(ifelse(active, c(rep("exponential", 3), "nugget"), "none"))
    grid <- build_cov_initial_grid(fixture$initial, fixture$data)
    expect_equal(rowSums(grid[, names]), rep(if (any(active)) 10 else 0, nrow(grid)))
    expect_true(all(as.matrix(grid[, names[!active], drop = FALSE]) == 0))
    expect_true(all(as.matrix(grid[, names[active], drop = FALSE]) > 0))
    expect_equal(nrow(unique(grid[, names])), if (any(active)) 2^sum(active) - 1 else 1)
    if (sum(active[1:3]) <= 2) {
      expect_equal(nrow(grid), max(1, 2^sum(active) - 1) * 2^sum(active[1:3]))
    }
  }
})

test_that("supplied starts override every expansion without changing known flags", {
  fixture <- grid_test_inputs(anisotropy = TRUE, random = "group", dispersion = 12)
  fixture$initial$tailup_initial <- tailup_initial_NA(tailup_initial("exponential", de = 20, known = "none"))
  fixture$initial$euclid_initial <- euclid_initial_NA(
    euclid_initial("exponential", range = 3, rotate = 0.2, scale = 0.7, known = "given"), fixture$data)
  fixture$initial$randcov_initial$initial <- c(group = 30)
  original <- fixture$initial
  grid <- build_cov_initial_grid(fixture$initial, fixture$data, is_glm = TRUE)
  expect_true(all(grid$tailup_de == 20 & grid$euclid_range == 3 & grid$group == 30))
  expect_true(all(grid$rotate == 0.2 & grid$scale == 0.7 & grid$dispersion == 12))
  expect_identical(fixture$initial, original)
  expect_false(fixture$initial$tailup_initial$is_known[["de"]])
  expect_true(fixture$initial$euclid_initial$is_known[["range"]])

  fixture$initial$tailup_initial <- tailup_initial_NA(tailup_initial("exponential", 20, 4, known = "none"))
  fixture$initial$taildown_initial <- taildown_initial_NA(taildown_initial("exponential", 2, 5, known = "given"))
  fixture$initial$euclid_initial$initial[["de"]] <- 3
  fixture$initial$nugget_initial <- nugget_initial_NA(nugget_initial("nugget", 7, known = "none"))
  expect_equal(nrow(build_cov_initial_grid(fixture$initial, fixture$data, is_glm = TRUE)), 1)
})

test_that("anisotropy and random effects replicate the full grid without a budget", {
  for (nrandom in 0:3) {
    random <- if (nrandom) paste0("group", seq_len(nrandom)) else character()
    fixture <- grid_test_inputs(anisotropy = TRUE, random = random)
    grid <- build_cov_initial_grid(fixture$initial, fixture$data)
    multiplier <- if (nrandom == 0) 1 else if (nrandom == 1) 3 else 1 + 2 * (nrandom + 1)
    expect_equal(nrow(grid), 96 * multiplier)
    variance <- c("tailup_de", "taildown_de", "euclid_de", "nugget", random)
    expect_equal(rowSums(grid[, variance]), rep(10, nrow(grid)))
    expect_equal(as.integer(table(grid$scale)), rep(48L * multiplier, 2))
    if (nrandom) {
      random_share <- round(rowSums(grid[, random, drop = FALSE]) / 10, 10)
      expect_equal(sort(unique(random_share)), c(0.1, 0.5, 0.9))
      for (share in c(0.1, 0.5, 0.9)) {
        block <- grid[random_share == share, ]
        normalized <- block[, 1:4] / (1 - share)
        expect_equal(nrow(unique(normalized)), 15)
      }
    }
  }
})

test_that("GLM dispersion and variance-start bounds are preserved", {
  fixture <- grid_test_inputs(dispersion = NA_real_)
  grid <- build_cov_initial_grid(fixture$initial, fixture$data, is_glm = TRUE)
  expect_equal(nrow(grid), 96)
  expect_equal(sort(unique(grid$dispersion)), c(1, 100))
  fixture$initial$dispersion_initial$initial[["dispersion"]] <- 1
  expect_equal(nrow(build_cov_initial_grid(fixture$initial, fixture$data, is_glm = TRUE)), 48)
  fixture$data$s2 <- 1e-5
  fixture$initial$taildown_initial <- taildown_initial_NA(taildown_initial("none"))
  grid <- build_cov_initial_grid(fixture$initial, fixture$data, is_glm = TRUE)
  expect_true(all(grid$taildown_de == 0))
  expect_true(all(as.matrix(grid[, c("tailup_de", "euclid_de", "nugget")]) >= 0.05))
})

test_that("shape expansion preserves spmodel's family-specific range corners", {
  for (type in c("matern", "cauchy", "pexponential")) {
    fixture <- grid_test_inputs(c("exponential", "exponential", type, "nugget"))
    grid <- build_cov_initial_grid(fixture$initial, fixture$data)
    expect_equal(nrow(grid), 96)
    for (shape in unique(grid$euclid_extra)) {
      rows <- grid[grid$euclid_extra == shape, ]
      ranges <- sort(unique(rows$euclid_range))
      expect_equal(ranges, c(10, 30) / if (type == "cauchy") 1 else 3)
    }
    fixture$initial$euclid_initial <- euclid_initial_NA(
      euclid_initial(type, extra = 1, range = 5, known = "none"), fixture$data)
    supplied <- build_cov_initial_grid(fixture$initial, fixture$data)
    expect_true(all(supplied$euclid_extra == 1 & supplied$euclid_range == 5))
  }
})

test_that("range scaling follows each Euclidean and stream kernel convention", {
  euclid <- c("exponential", "gaussian", "spherical", "circular", "cubic", "pentaspherical",
               "wave", "jbessel", "gravity", "rquad", "magnetic", "matern", "cauchy", "pexponential")
  factors <- c(3, sqrt(3), rep(1, 9), 3, 1, 3)
  for (i in seq_along(euclid)) {
    expect_equal(get_cov_initial_range(euclid_initial(euclid[i]), c(25, 75)), c(25, 75) / factors[i])
  }
  for (constructor in list(tailup_initial, taildown_initial)) {
    for (type in c("linear", "spherical", "epa", "mariah")) {
      expect_equal(get_cov_initial_range(constructor(type), c(25, 75)), c(25, 75))
    }
    expect_equal(get_cov_initial_range(constructor("exponential"), c(25, 75)), c(25, 75) / 3)
    ranges <- get_cov_initial_range(constructor("gaussian"), c(25, 75))
    u <- c(25, 75) / ranges
    expect_equal(2 * exp(-u^2) * (1 - pnorm(sqrt(2) * u)), rep(0.05, 2), tolerance = 1e-7)
  }
})

test_that("grid construction is deterministic and leaves the RNG untouched", {
  fixture <- grid_test_inputs(anisotropy = TRUE, random = c("a", "b"))
  set.seed(2)
  before <- .Random.seed
  grid <- build_cov_initial_grid(fixture$initial, fixture$data)
  expect_identical(.Random.seed, before)
  expect_identical(build_cov_initial_grid(fixture$initial, fixture$data), grid)
})

test_that("Gaussian selection agrees with independent GLS likelihood scores", {
  initial <- get_initial_object("exponential", "exponential", "exponential", "nugget",
                                NULL, NULL, NULL, NULL)
  data <- get_data_object(Summer_mn ~ ELEV_DEM, mf04p, "afvArea", FALSE,
                          initial, NULL, NULL, NULL, NULL)
  initial <- get_initial_NA_object(initial, data)
  grid <- build_cov_initial_grid(initial, data)
  X <- do.call(rbind, data$X_list)
  y <- do.call(rbind, data$y_list)
  score <- vapply(seq_len(nrow(grid)), function(i) {
    params <- get_params_object_grid(unlist(grid[i, ]), initial)
    V <- as.matrix(get_cov_matrix_list(params, data)[[1]])
    Vinv <- solve(V)
    beta <- solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv %*% y))
    resid <- y - X %*% beta
    as.numeric(nrow(X) * log(2 * pi) + determinant(V, logarithm = TRUE)$modulus +
                 crossprod(resid, Vinv %*% resid))
  }, numeric(1))
  selected <- cov_initial_search(initial, mf04p, data, "ml")$initial_object
  actual <- get_params_object_grid(unlist(grid[which.min(score), ]), initial)
  expect_equal(selected$tailup_initial$initial, unclass(actual$tailup))
  expect_equal(selected$taildown_initial$initial, unclass(actual$taildown))
  expect_equal(selected$euclid_initial$initial, unclass(actual$euclid))
  expect_equal(selected$nugget_initial$initial, unclass(actual$nugget))
})

# numerical policy
test_that("numerical fitting defaults follow spmodel", {
  expect_equal(get_optim_dotlist()$control$reltol, 1e-6)
  expect_equal(get_optim_dotlist(control = list(reltol = 1e-4))$control$reltol, 1e-4)
  expect_equal(get_optim_dotlist()$control$maxit, 2000)
  expect_equal(get_optim_dotlist(method = "Nelder-Mead")$control$maxit, 2000)
  expect_equal(get_optim_dotlist(control = list(maxit = 75))$control$maxit, 75)
  expect_null(get_optim_dotlist(method = "BFGS")$control$maxit)
  expect_equal(get_optim_dotlist(method = "BFGS", control = list(maxit = 75))$control$maxit, 75)
})

test_that("GLM data builders use the fixed diagonal tolerance", {
  s <- mf04p
  s$obs$policy_count <- pmax(0, round(s$obs$Summer_mn))
  initial <- get_initial_object_glm(
    "none", "none", "none", "nugget", NULL, NULL, NULL, NULL, "poisson", NULL
  )
  dense <- get_data_object_glm(
    policy_count ~ ELEV_DEM, s, "poisson", NULL, FALSE, initial,
    NULL, NULL, NULL, NULL
  )
  expect_equal(dense$diagtol, 1e-4)

  ssn_create_bigdist(s, overwrite = TRUE, no_cores = 1, verbose = FALSE)
  big <- get_data_object_bigdata_glm(
    policy_count ~ ELEV_DEM, s, "poisson", NULL, FALSE, initial,
    NULL, NULL, NULL, list(index = rep(1, nrow(s$obs)))
  )
  expect_equal(big$diagtol, 1e-4)
})

test_that("only estimated nuggets are reconciled to the covariance floor", {
  params <- list(
    tailup = tailup_params("exponential", de = 0.5, range = 100),
    taildown = taildown_params("none", de = 0, range = Inf),
    euclid = euclid_params("exponential", de = 1.5, range = 100, rotate = 0, scale = 1),
    nugget = nugget_params("nugget", nugget = 1e-12),
    randcov = NULL
  )
  estimated <- floor_estimated_nugget(params, c(nugget = FALSE), diagtol = 1e-4)
  known <- floor_estimated_nugget(params, c(nugget = TRUE), diagtol = 1e-4)

  expect_equal(estimated$nugget[["nugget"]], 2e-4)
  expect_equal(known$nugget[["nugget"]], 1e-12)
})

test_that("fitted estimated nuggets match the covariance matrix floor", {
  euclid <- euclid_initial("exponential", de = 2, range = 10000, known = "given")
  estimated <- ssn_lm(
    Summer_mn ~ ELEV_DEM, mf04p, euclid_type = "exponential", nugget_type = "nugget",
    euclid_initial = euclid, estmethod = "reml", control = list(reltol = 1e-6)
  )
  known <- ssn_lm(
    Summer_mn ~ ELEV_DEM, mf04p, euclid_type = "exponential", nugget_type = "nugget",
    euclid_initial = euclid,
    nugget_initial = nugget_initial("nugget", nugget = 1e-9, known = "given"),
    estmethod = "reml"
  )

  expect_equal(coef(estimated, type = "nugget")[["nugget"]], 2e-4, tolerance = 1e-12)
  expect_equal(covmatrix(estimated)[1, 1], 2.0002, tolerance = 1e-12)
  expect_equal(coef(known, type = "nugget")[["nugget"]], 1e-9, tolerance = 1e-15)
  expect_equal(covmatrix(known)[1, 1], 2.0002, tolerance = 1e-12)
})

