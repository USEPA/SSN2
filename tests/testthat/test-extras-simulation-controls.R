skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("simulation resolvers reject malformed and unknown option lists", {
  object <- list(n = 1L, ssn.object = list(obs = data.frame(pid = 1)))
  resolvers <- list(
    function(x) SSN2:::get_ssn_simulate_local(x, 1L, NULL),
    function(x) SSN2:::get_local_list_conditional(x, object)
  )
  malformed <- list(list(1), list(size = 1, size = 2),
    setNames(list(1), ""), setNames(list(1), NA_character_))
  removed <- c("approximation", "method_base", "size_base", "reorder_base",
    "method_new", "size_new", "reorder_new", "kmeans_new", "parallel", "ncores")
  for (resolve in resolvers) {
    for (value in list(NA, logical(), c(TRUE, FALSE), 1, "vecchia")) {
      expect_error(resolve(value), "local must be FALSE, TRUE, or a list", fixed = TRUE)
    }
    for (value in malformed) {
      expect_error(resolve(value), "nonempty, unique, nonmissing names", fixed = TRUE)
    }
    for (name in removed) {
      for (value in list(NULL, FALSE, "vecchia")) {
        expect_error(resolve(setNames(list(value), name)), "Unknown local option", fixed = TRUE)
      }
    }
    expect_true(resolve(FALSE)$exact)
    expect_true(resolve(c(use = FALSE))$exact)
    expect_false(resolve(c(use = TRUE))$exact)
    expect_false(resolve(list())$exact)
    expect_error(resolve(list(method = "all", ordering = "invalid")), "ordering must be")
  }
})

test_that("exact dispatch avoids ordering and local covariance contexts", {
  set.seed(2)
  local_mocked_bindings(
    get_decorrelate_order = function(...) stop("unexpected ordering"),
    get_ssn_simulate_covariance_fit = function(...) stop("unexpected local covariance"),
    .package = "SSN2"
  )
  network <- mf04p
  network$obs <- network$obs[c(3, 1, 2), ]
  args <- list(ssn.object = network,
    tailup_params = tailup_params("none"), taildown_params = taildown_params("none"),
    euclid_params = euclid_params("exponential", de = .4, range = 12000),
    nugget_params = nugget_params("nugget", nugget = .2), samples = 3,
    mean = c(-.2, 0, .4), size = c(2, 3, 4))
  for (family in c("gaussian", "poisson", "binomial", "nbinomial", "beta", "Gamma",
                   "inverse.gaussian")) {
    set.seed(2)
    exact <- do.call(ssn_simulate, c(args, list(family = family, local = FALSE)))
    exact_seed <- .Random.seed
    set.seed(2)
    all <- do.call(ssn_simulate, c(args, list(family = family,
      local = list(method = "all", size = NA, ordering = "random"))))
    expect_identical(all, exact, info = family)
    expect_identical(.Random.seed, exact_seed, info = family)
  }
  object <- list(n = 3L, ssn.object = network)
  expect_true(SSN2:::get_local_list_conditional(
    list(method = "all", ordering = "random", size = NA, chunk_size = NA),
    object, network$obs)$exact)
  expect_error(SSN2:::get_local_list_conditional(
    list(method = "all", size_random = 2), object), "only used for Vecchia")
})
