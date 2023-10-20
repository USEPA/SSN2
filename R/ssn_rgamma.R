#' @name ssn_simulate
#' @export
ssn_rgamma <- function(ssn.object, network = "obs",
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
  mu <- exp(ssn_rnorm_val)


  n <- NROW(ssn.object$obs)
  if (is.matrix(mu)) {
    mu_list <- split(t(mu), seq_len(NCOL(mu)))
    ssn_rgamma_val <- vapply(mu_list, function(x) rgamma(n, shape = dispersion, scale = x / dispersion), numeric(n))
  } else {
    ssn_rgamma_val <- rgamma(n, shape = dispersion, scale = mu / dispersion)
  }
  ssn_rgamma_val
}
