#' @name ssn_simulate
#' @export
ssn_rnorm <- function(ssn.object, network = "obs",
                      tailup_params, taildown_params, euclid_params, nugget_params,
                      mean = 0, samples = 1, additive,
                      randcov_params, partition_factor, local, ...) {
  if (any(!(network %in% "obs"))) {
    stop("network must be \"obs\".", call. = FALSE)
  }

  if (missing(local)) local <- NULL

  # fix additive depending on format
  if (missing(additive)) additive <- NULL
  if (is.symbol(substitute(additive))) { # or is.language
    additive <- deparse1(substitute(additive))
  }

  # save covariance types (for distance object later)
  tailup_type <- remove_covtype(class(tailup_params))
  taildown_type <- remove_covtype(class(taildown_params))
  euclid_type <- remove_covtype(class(euclid_params))
  nugget_type <- remove_covtype(class(nugget_params))

  # check on additive
  if (is.null(additive)) {
    if (tailup_type != "none") {
      stop("additive must be specified.", call. = FALSE)
    }
  }

  # compute the random effects design matrices
  if (missing(randcov_params)) {
    randcov_params <- NULL
    randcov_Zs <- NULL
  } else {
    names(randcov_params) <- get_randcov_names(reformulate(paste("(", names(randcov_params), ")", sep = "")))
    randcov_Zs <- get_randcov_Zs(data = ssn.object$obs, names(randcov_params))
  }

  # create the covariance parameter object
  params_object <- list(
    tailup = tailup_params, taildown = taildown_params,
    euclid = euclid_params, nugget = nugget_params, randcov = randcov_params
  )

  # perform relevant anisotropy checks if euclidean anisotropy required
  if (euclid_type != "none") {
    # find anisotropy
    if (params_object$euclid[["rotate"]] == 0 && params_object$euclid[["scale"]] == 1) {
      anisotropy <- FALSE
    } else {
      anisotropy <- TRUE
    }
  } else {
    anisotropy <- FALSE
  }

  # compute the partition matrix
  if (missing(partition_factor)) {
    partition_factor <- NULL
  }

  n <- NROW(ssn.object$obs)
  if (is.null(local)) {
    if (n > 5000) {
      local <- TRUE
      message(
        "Because the sample size exceeds 5,000, we are using a low-rank big-data approximation. ",
        "Set local = FALSE to use the exact simulation, or local = list(approximation = \"vecchia\") ",
        "for neighbor-truncated sequential conditioning instead."
      )
    } else {
      local <- FALSE
    }
  }
  conditioning <- get_ssn_simulate_local(local, n, ssn.object$obs)

  if (!identical(conditioning$method, "exact")) {
    covariance_fit <- get_ssn_simulate_covariance_fit(ssn.object, params_object, additive, anisotropy, partition_factor)
    ssn_rnorm_val <- if (identical(conditioning$approximation, "vecchia")) {
      get_ssn_simulate_vecchia(covariance_fit, samples, conditioning) + mean
    } else {
      get_ssn_simulate_lowrank(covariance_fit, samples, conditioning) + mean
    }
    if (samples == 1) {
      ssn_rnorm_val <- as.vector(ssn_rnorm_val)
    }
    return(ssn_rnorm_val)
  }

  partition_matrix_val <- partition_matrix(partition_factor, ssn.object$obs)

  # create the distance object required
  initial_object <- get_initial_object(
    tailup_type = tailup_type,
    taildown_type = taildown_type,
    euclid_type = euclid_type,
    nugget_type = nugget_type,
    tailup_initial = NULL,
    taildown_initial = NULL,
    euclid_initial = NULL,
    nugget_initial = NULL
  )
  # matches covmatrix()'s obs.obs backend resolution: prefer the standard
  # dense (.RData) distance matrices when they exist, falling back to the
  # big-data (.bmat) filematrix reader only when they do not
  tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
  taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")
  backend <- select_square_dist_backend(ssn.object, "obs", tailup_none, taildown_none)
  dist_object <- get_dist_object(ssn.object, initial_object, additive, anisotropy, backend = backend)

  # The covariance builder applies the same anisotropy transformation as fitting.
  de_scale <- sum(params_object$tailup[["de"]], params_object$taildown[["de"]], params_object$euclid[["de"]])
  cov_matrix_val <- get_cov_matrix(params_object, dist_object, randcov_Zs, partition_matrix_val,
    anisotropy, de_scale,
    diagtol = 0
  )

  # find the lower Cholesky needed for normal sim
  cov_matrix_lowchol <- t(chol(cov_matrix_val))

  # simulate n random normal vectors
  ssn_rnorm_val <- vapply(seq_len(samples), function(x) mean + as.numeric(cov_matrix_lowchol %*% rnorm(n)), numeric(n))

  # store as a vector if only one sample required
  if (samples == 1) {
    ssn_rnorm_val <- as.vector(ssn_rnorm_val)
  }

  # return
  ssn_rnorm_val
}
