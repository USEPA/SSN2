#' @name ssn_simulate
#' @export
ssn_rpois <- function(ssn.object, network = "obs",
                      tailup_params, taildown_params, euclid_params, nugget_params,
                      mean = 0, samples = 1, additive,
                      randcov_params, partition_factor, ...) {
  if (any(!(network %in% "obs"))) {
    stop("network must be \"obs\".", call. = FALSE)
  }

  # re-dispatch this call to ssn_rnorm() to simulate the shared latent
  # Gaussian process, then transform it below into the target distribution
  call_val <- match.call()
  call_val[[1]] <- as.symbol("ssn_rnorm")
  ssn_rnorm_val <- eval(call_val, envir = parent.frame())
  mu <- exp(ssn_rnorm_val)


  n <- NROW(ssn.object$obs)
  if (is.matrix(mu)) {
    # multiple samples: split the matrix into one mean vector per column
    # (transpose first so split() walks columns instead of rows), simulating
    # each sample's Poisson draws independently
    mu_list <- split(t(mu), seq_len(NCOL(mu)))
    ssn_rpois_val <- vapply(mu_list, function(x) rpois(n, x), numeric(n))
  } else {
    ssn_rpois_val <- rpois(n, mu)
  }
  ssn_rpois_val
}
