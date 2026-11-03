skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

simulation_split_network <- function() {
  network <- mf04p
  index <- rev(round(seq(1, NROW(network$obs), length.out = 16)))
  network$obs <- network$obs[index, ]
  network$obs$g <- factor(rep(c("a", "b"), 8))
  network$obs$h <- rep(c("u", "v", "w", "u"), 4)
  network$obs$x <- rep(c(-2, .5, 1, -2), 4)
  network$obs$part <- factor(rep(c("a|b", "a", "a|b", "a"), 4))
  network$obs$other <- rep(c("c", "b|c", "c", "b|c"), 4)
  network$obs$other[16] <- "singleton"
  network$path <- withr::local_tempdir(.local_envir = parent.frame())
  stopifnot(file.copy(file.path(mf04p$path, "binaryID.db"), network$path))
  ssn_create_distmat(network, overwrite = TRUE)
  network
}

simulation_split_params <- function(component = "mixed", effects = TRUE) {
  params <- list(tailup = tailup_params("none"),
    taildown = taildown_params("none"), euclid = euclid_params("none"),
    nugget = nugget_params("nugget", nugget = .2),
    randcov = if (effects) c("1 | g" = .4, "x | g" = .3, "1 | h" = .1) else NULL)
  if (component %in% c("up", "mixed", "wave"))
    params$tailup <- tailup_params("exponential", de = .4, range = 12000)
  if (component %in% c("down", "mixed", "wave"))
    params$taildown <- taildown_params("exponential", de = .3, range = 9000)
  if (component %in% c("euclid", "mixed", "anisotropic", "wave"))
    params$euclid <- euclid_params("exponential", de = .2, range = 16000)
  if (component == "anisotropic")
    params$euclid <- euclid_params("exponential", de = .7, range = 16000,
      rotate = .7, scale = .25)
  if (component %in% c("wave-up", "wave"))
    params$tailup <- tailup_params("wave", de = .4, range = 2500)
  if (component %in% c("wave-down", "wave"))
    params$taildown <- taildown_params("wave", de = .3, range = 3500)
  params
}

simulation_split_fit <- function(network, params, partition = NULL) {
  anisotropy <- !inherits(params$euclid, "euclid_none") &&
    (params$euclid[["rotate"]] != 0 || params$euclid[["scale"]] != 1)
  SSN2:::get_ssn_simulate_covariance_fit(network, params, "afvArea", anisotropy, partition)
}

simulation_split_native <- function(fit, backend = "dense") {
  params <- fit$coefficients$params_object
  types <- vapply(params[1:4], function(x) SSN2:::remove_covtype(class(x)), character(1))
  initial <- SSN2:::get_initial_object(types[1], types[2], types[3], types[4],
    NULL, NULL, NULL, NULL)
  distance <- SSN2:::get_dist_object(fit$ssn.object, initial, fit$additive,
    fit$anisotropy, backend = backend)
  random <- SSN2:::get_randcov_Zs(fit$ssn.object$obs, names(params$randcov))
  partition <- SSN2:::partition_matrix(fit$partition_factor, fit$ssn.object$obs)
  de <- sum(vapply(params[1:3], function(x) x[["de"]], numeric(1)))
  unname(as.matrix(SSN2:::get_cov_matrix(params, distance, random, partition,
    fit$anisotropy, de, diagtol = 0)))
}

simulation_split_args <- function(fit, samples, local, mean = 0) {
  params <- fit$coefficients$params_object
  random <- params$randcov
  if (length(random)) names(random) <- sub("^x ", "0 + x ", names(random))
  list(ssn.object = fit$ssn.object, tailup_params = params$tailup,
    taildown_params = params$taildown, euclid_params = params$euclid,
    nugget_params = params$nugget, randcov_params = random,
    partition_factor = fit$partition_factor, additive = fit$additive,
    mean = mean, samples = samples, local = local)
}

simulation_split_basis <- function(fit, size, ordering = "pid", mean = 0) {
  n <- NROW(fit$ssn.object$obs)
  design <- SSN2:::get_randcov_design(fit$ssn.object$obs,
    fit$coefficients$params_object$randcov, fit$partition_factor)
  q <- NCOL(design$observed)
  samples <- n + q
  # Disjoint columns make spatial and coefficient innovations orthogonal.
  innovations <- c(cbind(diag(n), matrix(0, n, q)),
    cbind(matrix(0, q, n), diag(q)))
  used <- 0L
  local_mocked_bindings(rnorm = function(n, mean = 0, sd = 1) {
    index <- used + seq_len(n)
    used <<- used + n
    stopifnot(used <= length(innovations))
    mean + sd * innovations[index]
  }, .package = "SSN2")
  set.seed(2)
  seed <- .Random.seed
  value <- do.call(ssn_rnorm, simulation_split_args(fit, samples,
    list(method = "covariance", size = size, ordering = ordering), mean))
  expect_identical(used, length(innovations))
  expect_identical(.Random.seed, seed)
  value
}

simulation_split_expected_neighbors <- function(covariance, permutation, size, groups) {
  lapply(seq_along(permutation), function(position) {
    pool <- seq_len(position - 1L)
    pool <- pool[groups[permutation[pool]] == groups[permutation[position]]]
    cross <- covariance[permutation[position], permutation[pool]]
    if (length(pool) > size) {
      keep <- order(-abs(cross), pool)[seq_len(size)]
      pool <- pool[keep]
      cross <- cross[keep]
    }
    list(positions = pool, covariance = unname(cross))
  })
}

test_that("full spatial neighborhoods and exact effects recover native total covariance", {
  network <- simulation_split_network()
  identity <- ssn_get_netgeom(network$obs, reformat = TRUE)
  expect_gt(length(unique(identity$NetworkID)), 1L)
  expect_true(is.unsorted(identity$pid))
  expect_true(any(diff(sort(identity$pid)) > 1))
  for (partition in list(NULL, ~part + other)) {
    for (component in c("up", "down", "euclid", "mixed", "anisotropic")) {
      fit <- simulation_split_fit(network, simulation_split_params(component), partition)
      original <- fit
      reference <- simulation_split_native(fit)
      draws <- simulation_split_basis(fit, NROW(network$obs) - 1L)
      expect_equal(tcrossprod(draws), reference, tolerance = 1e-12, info = component)
      expect_equal(sum(colMeans(draws)^2), mean(reference), tolerance = 1e-12)
      expect_identical(fit, original)
    }
  }
})

test_that("full wave neighborhoods recover native total covariance", {
  skip_if_not_installed("pracma")
  network <- simulation_split_network()
  for (partition in list(NULL, ~part + other)) {
    for (component in c("wave-up", "wave-down", "wave")) {
      fit <- simulation_split_fit(network, simulation_split_params(component), partition)
      reference <- simulation_split_native(fit)
      draws <- simulation_split_basis(fit, 15L)
      expect_equal(tcrossprod(draws), reference, tolerance = 1e-12, info = component)
      expect_equal(sum(colMeans(draws)^2), mean(reference), tolerance = 1e-12)
    }
  }
})

test_that("coincident Euclidean sites retain the native nugget floor", {
  network <- simulation_split_network()
  sf::st_geometry(network$obs)[2] <- sf::st_geometry(network$obs)[1]
  params <- simulation_split_params("euclid", effects = FALSE)
  for (nugget in list(nugget_params("none"), nugget_params("nugget", nugget = 0))) {
    params$nugget <- nugget
    fit <- simulation_split_fit(network, params)
    reference <- simulation_split_native(fit)
    draws <- simulation_split_basis(fit, 15L, "none", mean = seq_len(16) / 16)
    expect_equal(tcrossprod(draws - seq_len(16) / 16), reference, tolerance = 1e-12)
    expect_equal(reference[1, 1] - reference[1, 2], .2 * 1e-4, tolerance = 1e-12)
  }
})

test_that("spatial neighbors use native absolute covariance and partition predecessors", {
  network <- simulation_split_network()
  permutation <- c(16L, seq_len(15L))
  for (partition in list(NULL, ~part + other)) {
    for (component in c("up", "down", "mixed", "anisotropic")) {
      fit <- simulation_split_fit(network, simulation_split_params(component, FALSE), partition)
      covariance <- simulation_split_native(fit)
      groups <- if (is.null(partition)) rep(1L, 16) else
        as.integer(interaction(network$obs$part, network$obs$other, drop = TRUE))
      for (size in c(1L, 3L, 15L)) {
        reference <- simulation_split_expected_neighbors(covariance, permutation, size, groups)
        set.seed(2)
        seed <- .Random.seed
        actual <- SSN2:::get_ssn_simulate_neighbors(fit, permutation, size)
        expect_equal(actual, reference, tolerance = 1e-12, info = paste(component, size))
        expect_identical(.Random.seed, seed)
        expect_length(actual[[1]]$positions, 0L)
        if (!is.null(partition)) expect_length(actual[[2]]$positions, 0L)
      }
    }
  }
  initial <- SSN2:::get_initial_object("exponential", "exponential", "none", "nugget",
    NULL, NULL, NULL, NULL)
  distance <- SSN2:::get_dist_object(network, initial, "afvArea", FALSE)
  same <- as.matrix(distance$mask_mat) == 1 & row(distance$a_mat) != col(distance$a_mat)
  expect_true(any(distance$a_mat[same] == 0))
  expect_true(any(distance$a_mat[same] > 0))
})

test_that("wave neighbor ranking includes native negative lobes", {
  skip_if_not_installed("pracma")
  network <- simulation_split_network()
  permutation <- c(16L, seq_len(15L))
  for (partition in list(NULL, ~part + other)) {
    for (component in c("wave-up", "wave-down")) {
      fit <- simulation_split_fit(network, simulation_split_params(component, FALSE), partition)
      covariance <- simulation_split_native(fit)
      groups <- if (is.null(partition)) rep(1L, 16) else
        as.integer(interaction(network$obs$part, network$obs$other, drop = TRUE))
      for (size in c(1L, 3L)) {
        actual <- SSN2:::get_ssn_simulate_neighbors(fit, permutation, size)
        expect_equal(actual, simulation_split_expected_neighbors(covariance,
          permutation, size, groups), tolerance = 1e-12)
      }
      if (component == "wave-down" && is.null(partition))
        expect_true(any(vapply(actual, function(x) any(x$covariance < 0), logical(1))))
    }
  }
})

test_that("zero covariance ties retain earlier positions and whole pools retain order", {
  network <- simulation_split_network()
  params <- simulation_split_params("up", FALSE)
  params$tailup <- tailup_params("exponential", de = .4, range = 1e-6)
  fit <- simulation_split_fit(network, params, ~part + other)
  permutation <- rev(seq_len(16))
  groups <- as.integer(interaction(network$obs$part, network$obs$other, drop = TRUE))
  covariance <- simulation_split_native(fit)
  expect_equal(covariance - diag(diag(covariance)), matrix(0, 16, 16))
  for (size in c(1L, 3L, 15L)) {
    actual <- SSN2:::get_ssn_simulate_neighbors(fit, permutation, size)
    expect_equal(actual, simulation_split_expected_neighbors(covariance, permutation, size, groups))
    expect_true(all(vapply(actual, function(x) !is.unsorted(x$positions), logical(1))))
  }
})

test_that("anisotropy changes native neighbor ranking without a network filter", {
  network <- simulation_split_network()
  network$obs <- network$obs[c(1, 16, 2), ]
  sf::st_geometry(network$obs) <- sf::st_sfc(sf::st_point(c(1, 0)),
    sf::st_point(c(0, 2)), sf::st_point(c(0, 0)), crs = sf::st_crs(network$obs))
  params <- simulation_split_params("euclid", FALSE)
  params$euclid <- euclid_params("exponential", de = .7, range = 10)
  isotropic <- simulation_split_fit(network, params)
  params$euclid <- euclid_params("exponential", de = .7, range = 10, rotate = pi / 2, scale = .1)
  anisotropic <- simulation_split_fit(network, params)
  plain <- SSN2:::get_ssn_simulate_neighbors(isotropic, 1:3, 1L)
  transformed <- SSN2:::get_ssn_simulate_neighbors(anisotropic, 1:3, 1L)
  expect_identical(plain[[3]]$positions, 1L)
  expect_identical(transformed[[3]]$positions, 2L)
  ids <- ssn_get_netgeom(network$obs, reformat = TRUE)$NetworkID
  expect_false(ids[2] == ids[3])
  expect_gt(transformed[[3]]$covariance, 0)
  expect_equal(transformed, simulation_split_expected_neighbors(
    simulation_split_native(anisotropic), 1:3, 1L, rep(1L, 3)), tolerance = 1e-12)
})

test_that("zero cross covariance candidates can still affect spatial conditioning", {
  network <- simulation_split_network()
  network$obs <- network$obs[1:4, ]
  sf::st_geometry(network$obs) <- sf::st_sfc(lapply(c(.6, 1.2, 3, 0),
    function(x) sf::st_point(c(x, 0))), crs = sf::st_crs(network$obs))
  params <- simulation_split_params("euclid", FALSE)
  params$euclid <- euclid_params("spherical", de = 1, range = 1)
  fit <- simulation_split_fit(network, params)
  covariance <- simulation_split_native(fit)
  neighbors <- SSN2:::get_ssn_simulate_neighbors(fit, 1:4, 2L)
  expect_identical(neighbors[[4]]$positions, 1:2)
  expect_equal(neighbors[[4]]$covariance[2], 0)
  expect_gt(covariance[1, 2], 0)
  one <- tcrossprod(simulation_split_basis(fit, 1L, "none"))
  two <- tcrossprod(simulation_split_basis(fit, 2L, "none"))
  expect_gt(one[4, 2], 0)
  expect_equal(two[4, 2], 0, tolerance = 1e-12)
})

test_that("truncation preserves shared intercepts and signed slopes across networks", {
  network <- simulation_split_network()
  for (partition in list(NULL, ~part + other)) {
    params <- simulation_split_params("mixed")
    total <- simulation_split_fit(network, params, partition)
    params$randcov <- NULL
    spatial <- simulation_split_fit(network, params, partition)
    random <- simulation_split_native(total) - simulation_split_native(spatial)
    same <- outer(as.character(network$obs$g), as.character(network$obs$g), "==")
    expected <- same * (.4 + .3 * outer(network$obs$x, network$obs$x)) +
      .1 * outer(network$obs$h, network$obs$h, "==")
    if (!is.null(partition)) expected <- expected *
      outer(as.character(network$obs$part), as.character(network$obs$part), "==") *
      outer(network$obs$other, network$obs$other, "==")
    expect_equal(random, expected, tolerance = 1e-12)
    for (size in c(1L, 3L)) {
      draws <- simulation_split_basis(total, size)
      residual <- simulation_split_basis(spatial, size)
      difference <- tcrossprod(draws) - tcrossprod(residual)
      expect_equal(difference, expected, tolerance = 1e-12)
      expect_equal(sum(colSums(draws)^2) - sum(colSums(residual)^2),
        sum(expected), tolerance = 1e-12)
      expect_equal(draws[, seq_len(16)], residual, tolerance = 1e-12)
    }
    ids <- ssn_get_netgeom(network$obs, reformat = TRUE)$NetworkID
    expect_true(any(expected[outer(ids, ids, "!=")] > 0))
    expect_true(any(expected < 0))
  }
})

test_that("spatial setup sees no random effects and neighbors are cached once", {
  network <- simulation_split_network()
  fit <- simulation_split_fit(network, simulation_split_params())
  neighbors <- SSN2:::get_ssn_simulate_neighbors
  cross <- SSN2:::get_decorrelate_observed_cross_covariance
  pool <- SSN2:::get_decorrelate_observed_covariance
  marginal <- SSN2:::get_decorrelate_marginal_variance
  calls <- 0L
  spatial_only <- function(fit) {
    expect_null(fit$coefficients$params_object$randcov)
    expect_null(fit$random_xlev)
  }
  local_mocked_bindings(get_ssn_simulate_neighbors = function(covariance_fit, permutation, size) {
    spatial_only(covariance_fit)
    calls <<- calls + 1L
    neighbors(covariance_fit, permutation, size)
  }, get_decorrelate_observed_cross_covariance = function(covariance_fit, d1, d2) {
    spatial_only(covariance_fit)
    cross(covariance_fit, d1, d2)
  }, get_decorrelate_observed_covariance = function(covariance_fit, data) {
    spatial_only(covariance_fit)
    pool(covariance_fit, data)
  }, get_decorrelate_marginal_variance = function(covariance_fit, data) {
    spatial_only(covariance_fit)
    marginal(covariance_fit, data)
  }, .package = "SSN2")
  set.seed(2)
  draws <- do.call(ssn_rnorm, simulation_split_args(fit, 7L, list(size = 3L)))
  expect_equal(dim(draws), c(16L, 7L))
  expect_identical(calls, 1L)
})

test_that("dense and filematrix distances preserve simulation identity and row order", {
  skip_if_not_installed("filematrix")
  network <- simulation_split_network()
  fit <- simulation_split_fit(network, simulation_split_params("mixed"), ~part)
  reference <- simulation_split_native(fit)
  dense <- simulation_split_basis(fit, 15L)
  ssn_create_bigdist(network, overwrite = TRUE, no_cores = 1, verbose = FALSE)
  expect_equal(simulation_split_basis(fit, 15L), dense, tolerance = 1e-12)
  files <- list.files(file.path(network$path, "distance", "obs"),
    pattern = "[.]RData$", full.names = TRUE)
  expect_gt(length(files), 0L)
  expect_true(all(file.remove(files)))
  expect_identical(SSN2:::get_decorrelate_observed_backend(fit), "bigdata")
  expect_equal(simulation_split_native(fit, "bigdata"), reference, tolerance = 1e-12)
  expect_equal(tcrossprod(simulation_split_basis(fit, 15L)), reference, tolerance = 1e-12)
})

test_that("all dispatch is seeded exact simulation without ordering or local setup", {
  network <- simulation_split_network()
  fit <- simulation_split_fit(network, simulation_split_params("anisotropic"), ~other)
  unexpected <- function(...) stop("unexpected local setup")
  local_mocked_bindings(get_decorrelate_order = unexpected,
    get_ssn_simulate_covariance_fit = unexpected,
    get_ssn_simulate_neighbors = unexpected, .package = "SSN2")
  for (samples in c(1L, 3L)) {
    set.seed(2)
    exact <- do.call(ssn_rnorm, simulation_split_args(fit, samples, FALSE, 2))
    seed <- .Random.seed
    set.seed(2)
    all <- do.call(ssn_rnorm, simulation_split_args(fit, samples,
      list(method = "all", size = -1, ordering = "random"), 2))
    expect_identical(all, exact)
    expect_identical(.Random.seed, seed)
  }
})

simulation_split_iid <- function(n = 16L) {
  list(ssn.object = list(obs = data.frame(g = factor(seq_len(n)),
    h = rep(c("a", "b"), length.out = n), x = rep(c(-1, 2), length.out = n),
    part = rep(c("p", "q"), length.out = n), afvArea = rep(1, n))),
    tailup_params = tailup_params("none"), taildown_params = taildown_params("none"),
    euclid_params = euclid_params("none"), nugget_params = nugget_params("nugget", nugget = .4),
    randcov_params = c("1 | g" = .3, "0 + x | h" = .2), partition_factor = ~part,
    additive = "afvArea")
}

test_that("iid effects use sparse designs without covariance or distance setup", {
  unexpected <- function(...) stop("unexpected dense or spatial setup")
  validate <- SSN2:::get_randcov_Zs
  draw <- SSN2:::draw_ssn_simulate_random
  local_mocked_bindings(get_ssn_simulate_local = unexpected,
    get_ssn_simulate_covariance_fit = unexpected, get_decorrelate_order = unexpected,
    get_dist_object = unexpected, get_dist_object_bigdata_cross = unexpected,
    select_square_dist_backend = unexpected, get_cov_matrix = unexpected,
    partition_matrix = unexpected, chol = unexpected,
    get_randcov_Zs = function(data, randcov_names, ZZt = TRUE, ...) {
      expect_false(ZZt)
      validate(data, randcov_names, ZZt = ZZt, ...)
    }, draw_ssn_simulate_random = function(design, samples) {
      expect_s4_class(design$observed, "sparseMatrix")
      expect_lt(as.numeric(object.size(design$observed)), 1e6)
      draw(design, samples)
    }, .package = "SSN2")
  for (n in c(16L, 6001L)) for (samples in c(1L, 3L)) {
    args <- simulation_split_iid(n)
    design <- SSN2:::get_randcov_design(args$ssn.object$obs, args$randcov_params, ~part)
    q <- NCOL(design$observed)
    for (mean in list(2, seq_len(n) / n)) {
      set.seed(2)
      reference <- matrix(rnorm(n * samples, sd = sqrt(.4)), n, samples) + mean
      coefficients <- matrix(rnorm(q * samples), q, samples) * sqrt(design$variance)
      reference <- reference + as.matrix(design$observed %*% coefficients)
      if (samples == 1) reference <- as.vector(reference)
      seed <- .Random.seed
      for (local in list(FALSE, TRUE, list(size = -1, ordering = "unused"))) {
        set.seed(2)
        actual <- do.call(ssn_rnorm, c(args, list(samples = samples, mean = mean, local = local)))
        expect_equal(actual, reference, tolerance = 1e-12)
        expect_identical(.Random.seed, seed)
      }
    }
  }
  for (component in c("tailup", "taildown", "euclid")) {
    args <- simulation_split_iid()
    constructor <- get(paste0(component, "_params"), asNamespace("SSN2"))
    args[[paste0(component, "_params")]] <- constructor("exponential", de = 0, range = 1)
    set.seed(2)
    expect_length(do.call(ssn_rnorm, args), 16L)
  }
})

test_that("empty and zero-variance random designs preserve their draw contracts", {
  observed <- simulation_split_iid()$ssn.object$obs
  empty <- SSN2:::get_randcov_design(observed, NULL)
  set.seed(2)
  seed <- .Random.seed
  expect_null(SSN2:::draw_ssn_simulate_random(empty, 3L))
  expect_identical(.Random.seed, seed)
  design <- SSN2:::get_randcov_design(observed, c("1 | g" = 0, "x | h" = .2), ~part)
  q <- NCOL(design$observed)
  for (samples in c(1L, 3L)) {
    set.seed(2)
    expected <- unname(as.matrix(design$observed %*%
      (matrix(rnorm(q * samples), q, samples) * sqrt(design$variance))))
    seed <- .Random.seed
    set.seed(2)
    actual <- SSN2:::draw_ssn_simulate_random(design, samples)
    expect_identical(actual, expected)
    expect_identical(.Random.seed, seed)
    expect_equal(dim(actual), c(16L, samples))
  }
})

test_that("iid random effects retain group slope partition and network validation", {
  args <- simulation_split_iid()
  expect_error(do.call(ssn_rnorm, c(args, list(network = "pred1km"))), "network must be")
  for (column in c("g", "x", "part")) {
    invalid <- args
    invalid$ssn.object$obs[[column]][1] <- NA
    expect_error(suppressWarnings(do.call(ssn_rnorm, invalid)))
  }
  invalid <- args
  invalid$ssn.object$obs$g <- seq_len(16)
  expect_error(do.call(ssn_rnorm, invalid), "categorical|factor")
  invalid <- args
  invalid$ssn.object$obs$x <- factor(rep(letters[1:3], length.out = 16))
  expect_error(do.call(ssn_rnorm, invalid), "numeric|slope")
  invalid <- args
  invalid$ssn.object$obs$x[1] <- Inf
  expect_error(do.call(ssn_rnorm, invalid), "slope")
  invalid <- args
  invalid$partition_factor <- ~absent
  expect_error(do.call(ssn_rnorm, invalid), "absent")
  invalid <- args
  invalid$tailup_params <- tailup_params("exponential", de = 0, range = 1)
  invalid$additive <- NULL
  expect_error(do.call(ssn_rnorm, invalid), "additive")
  invalid$additive <- "absent"
  expect_error(do.call(ssn_rnorm, invalid), "additive|absent")
  invalid$additive <- "afvArea"
  invalid$ssn.object$obs$afvArea[1] <- Inf
  expect_error(do.call(ssn_rnorm, invalid), "additive|finite")
  for (component in c("tailup", "taildown", "euclid")) {
    constructor <- get(paste0(component, "_params"), asNamespace("SSN2"))
    expect_error(constructor("exponential", de = 0, range = 0), "range")
    invalid <- args
    invalid[[paste0(component, "_params")]] <- constructor("exponential", de = 0, range = 1)
    invalid[[paste0(component, "_params")]][["range"]] <- 0
    expect_error(do.call(ssn_rnorm, invalid), "range")
  }
})

test_that("family wrappers preserve the latent interface and response support", {
  skip_if_not_installed("statmod")
  network <- simulation_split_network()
  fit <- simulation_split_fit(network, simulation_split_params("mixed"), ~part)
  families <- c("gaussian", "poisson", "nbinomial", "binomial", "beta", "Gamma", "inverse.gaussian")
  for (samples in c(1L, 3L)) {
    args <- simulation_split_args(fit, samples, list(size = 1L, ordering = "none"), mean = -.5)
    set.seed(2)
    latent <- do.call(ssn_rnorm, args)
    for (family in families) {
      reference <- local({
        set.seed(2)
        do.call(ssn_rnorm, args)
        local_mocked_bindings(ssn_rnorm = function(ssn.object, mean, samples,
          randcov_params, partition_factor, ...) {
          expect_identical(ssn.object, network)
          expect_identical(mean, args$mean)
          expect_identical(samples, args$samples)
          expect_identical(randcov_params, args$randcov_params)
          expect_identical(partition_factor, args$partition_factor)
          expect_identical(list(...)$local, args$local)
          latent
        }, .package = "SSN2")
        do.call(ssn_simulate, c(list(family = family, size = 4, dispersion = 2), args))
      })
      set.seed(2)
      actual <- do.call(ssn_simulate, c(list(family = family, size = 4, dispersion = 2), args))
      expect_identical(actual, reference, info = family)
      expect_true(all(is.finite(actual)))
      if (family %in% c("poisson", "nbinomial", "binomial"))
        expect_true(all(actual >= 0 & actual == floor(actual)))
      if (family == "binomial") expect_true(all(actual <= 4))
      if (family == "beta") expect_true(all(actual > 0 & actual < 1))
      if (family %in% c("Gamma", "inverse.gaussian")) expect_true(all(actual > 0))
    }
  }
})
