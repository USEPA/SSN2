test_that("latent coefficient covariance does not depend on local variance adjustment", {
  set.seed(2)
  network <- mf04p
  n <- nrow(network$obs)
  network$obs$x <- as.numeric(scale(network$obs$ELEV_DEM))
  index <- rep(seq_len(4L), each = ceiling(n / 4))[seq_len(n)]
  for (family in c("poisson", "binomial")) {
    network$obs$y <- if (family == "poisson") {
      rpois(n, exp(1 + .2 * network$obs$x))
    } else {
      rbinom(n, 1, plogis(.2 * network$obs$x))
    }
    fits <- lapply(c("none", "theoretical", "pooled", "empirical"), function(adjustment) {
      ssn_glm(y ~ x, network, family = as.character(family),
        euclid_initial = euclid_initial("exponential", .5, 10000, known = "given"),
        nugget_initial = nugget_initial("nugget", .2, known = "given"),
        local = list(index = index, var_adjust = adjustment))
    })
    reference <- fits[[1]]
    latent <- vcov(reference) - vcov(reference, var_correct = FALSE)
    for (fit in fits[-1]) {
      expect_identical(coef(fit), coef(reference))
      expect_identical(fit$coefficients$params_object, reference$coefficients$params_object)
      expect_equal(vcov(fit) - vcov(fit, var_correct = FALSE), latent, tolerance = 1e-10)
    }
  }
})
