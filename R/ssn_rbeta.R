#' @name ssn_simulate
#' @export
ssn_rbeta <- function(ssn.object, network = "obs",
                      tailup_params, taildown_params, euclid_params, nugget_params,
                      dispersion = 1, mean = 0, samples = 1, additive,
                      randcov_params, partition_factor, ...) {
  if (any(!(network %in% "obs"))) {
    stop("network must be \"obs\".", call. = FALSE)
  }

  call_val <- match.call()
  call_val[[1]] <- as.symbol("ssn_rnorm")
  call_list <- as.list(call_val)
  if ("dispersion" %in% names(call_list)) {
    call_list <- call_list[-which(names(call_list) == "dispersion")]
  }
  call_val <- as.call(call_list)
  ssn_rnorm_val <- eval(call_val, envir = parent.frame())
  mu <- expit(ssn_rnorm_val)

  n <- NROW(ssn.object$obs)
  if (is.matrix(mu)) {
    mu_list <- split(t(mu), seq_len(NCOL(mu)))
    ssn_rbeta_val <- vapply(mu_list, function(x) {
      a <- x * dispersion
      b <- (1 - x) * dispersion
      val <- rbeta(n, shape1 = a, shape2 = b)
      val <- pmax(1e-4, val)
      val <- pmin(1 - 1e-4, val)
    }, numeric(n))
  } else {
    a <- mu * dispersion
    b <- (1 - mu) * dispersion
    ssn_rbeta_val <- rbeta(n, shape1 = a, shape2 = b)
    ssn_rbeta_val <- pmax(1e-4, ssn_rbeta_val)
    ssn_rbeta_val <- pmin(1 - 1e-4, ssn_rbeta_val)
  }
  ssn_rbeta_val
}
