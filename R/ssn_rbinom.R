#' @name ssn_simulate
#' @export
ssn_rbinom <- function(ssn.object, network = "obs",
                       tailup_params, taildown_params, euclid_params, nugget_params,
                       mean = 0, size = 1, samples = 1, additive,
                       randcov_params, partition_factor, ...) {
  if (any(!(network %in% "obs"))) {
    stop("network must be \"obs\".", call. = FALSE)
  }

  call_val <- match.call()
  call_val[[1]] <- as.symbol("ssn_rnorm")
  call_list <- as.list(call_val)
  if ("size" %in% names(call_list)) {
    call_list <- call_list[-which(names(call_list) == "size")]
  }
  call_val <- as.call(call_list)
  ssn_rnorm_val <- eval(call_val, envir = parent.frame())
  mu <- expit(ssn_rnorm_val)


  n <- NROW(ssn.object$obs)
  if (is.matrix(mu)) {
    mu_list <- split(t(mu), seq_len(NCOL(mu)))
    ssn_rbinom_val <- vapply(mu_list, function(x) rbinom(n, size, x), numeric(n))
  } else {
    ssn_rbinom_val <- rbinom(n, size, mu)
  }
  ssn_rbinom_val
}
