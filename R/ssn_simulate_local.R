get_local_vecchia_settings <- function(local, default_ordering = "pid") {
  if (is.null(local$method)) {
    local$method <- "covariance"
  } else if (length(local$method) != 1L || is.na(local$method) || !local$method %in% c("all", "covariance")) {
    stop("local$method must be \"all\" or \"covariance\".", call. = FALSE)
  }
  if (is.null(local$size)) local$size <- 30L
  if (identical(local$method, "covariance") &&
      (length(local$size) != 1L || !is.numeric(local$size) ||
       !is.finite(local$size) || local$size < 1 || local$size != floor(local$size))) {
    stop("local$size must be a positive integer.", call. = FALSE)
  }
  if (is.null(local$ordering)) local$ordering <- default_ordering
  local$ordering <- get_decorrelate_ordering(local$ordering)
  local
}

get_ssn_simulate_local <- function(local, n, data) {
  if (is.null(local)) {
    local <- n > 5000
    if (local) {
      message("Because the desired simulation size exceeds 5,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    }
  }
  local <- validate_simulation_local(local, c("method", "size", "ordering"))
  if (identical(local, FALSE)) return(list(exact = TRUE, method = "all"))
  if (identical(local, TRUE)) local <- list()
  local <- get_local_vecchia_settings(local)
  if (local$method == "all") return(list(exact = TRUE, method = "all"))
  local$exact <- FALSE
  local
}

get_ssn_simulate_covariance_fit <- function(ssn.object, params_object, additive, anisotropy, partition_factor) {
  randcov_names <- names(params_object$randcov)
  random_xlev <- if (is.null(randcov_names)) NULL else get_randcov_xlev(randcov_names, ssn.object$obs)
  partition_xlev <- NULL
  if (!is.null(partition_factor)) {
    formula <- reformulate(labels(terms(partition_factor)), intercept = FALSE,
      env = asNamespace("SSN2"))
    frame <- model.frame(formula, ssn.object$obs, na.action = na.pass)
    partition_xlev <- .getXlevels(terms(frame), frame)
  }
  list(
    ssn.object = ssn.object,
    coefficients = list(params_object = params_object),
    additive = additive,
    anisotropy = anisotropy,
    partition_factor = partition_factor,
    diagtol = 0,
    random_xlev = random_xlev,
    partition_xlev = partition_xlev
  )
}

get_ssn_simulate_vecchia <- function(covariance_fit, samples, conditioning) {
  observed <- covariance_fit$ssn.object$obs
  n <- NROW(observed)
  design <- get_randcov_design(observed, covariance_fit$coefficients$params_object$randcov,
    covariance_fit$partition_factor)
  covariance_fit$coefficients$params_object$randcov <- NULL
  covariance_fit$random_xlev <- NULL
  rows <- get_decorrelate_rows(observed)
  permutation <- get_decorrelate_order(rows, conditioning$ordering, observed)
  neighbors <- get_ssn_simulate_neighbors(covariance_fit, permutation, conditioning$size)
  # Each neighborhood lies within one partition, so its covariance needs no mask.
  covariance_fit$partition_factor <- NULL
  covariance_fit$partition_xlev <- NULL

  Z <- matrix(rnorm(n * samples), n, samples)
  Y_ordered <- matrix(NA_real_, n, samples)

  for (position in seq_len(n)) {
    current <- permutation[[position]]
    variance <- get_decorrelate_marginal_variance(covariance_fit, observed[current, , drop = FALSE])

    neighbor_positions <- neighbors[[position]]$positions
    if (!length(neighbor_positions)) {
      Y_ordered[position, ] <- sqrt(variance) * Z[position, ]
      next
    }

    neighbor_index <- permutation[neighbor_positions]

    covariance_neighbors <- get_decorrelate_observed_covariance(
      covariance_fit, observed[neighbor_index, , drop = FALSE]
    )
    chol_pool <- tryCatch(chol(covariance_neighbors), error = function(error) NULL)
    if (is.null(chol_pool)) {
      stop(
        "The neighbor covariance matrix for simulated record ", position,
        " (of the requested ordering) is not numerically positive definite; check covariance parameters and duplicate locations.",
        call. = FALSE
      )
    }

    c_vec <- neighbors[[position]]$covariance
    w <- backsolve(chol_pool, forwardsolve(t(chol_pool), c_vec))
    cond_var <- max(variance - sum(w * c_vec), 0)
    pool_val <- Y_ordered[neighbor_positions, , drop = FALSE]
    cond_mean <- as.numeric(crossprod(w, pool_val))

    Y_ordered[position, ] <- cond_mean + sqrt(cond_var) * Z[position, ]
  }

  Y <- matrix(NA_real_, n, samples)
  Y[permutation, ] <- Y_ordered
  random <- draw_ssn_simulate_random(design, samples)
  if (!is.null(random)) Y <- Y + random
  Y
}

get_ssn_simulate_partition <- function(observed, partition_factor) {
  if (is.null(partition_factor)) return(rep("", NROW(observed)))
  formula <- reformulate(labels(terms(partition_factor)), intercept = FALSE,
    env = asNamespace("SSN2"))
  frame <- model.frame(formula, observed, na.action = na.pass)
  if (anyNA(frame)) {
    stop("Missing values in partition_factor are not supported.", call. = FALSE)
  }
  randcov_design_key(frame)
}

get_ssn_simulate_neighbors <- function(covariance_fit, permutation, size) {
  observed <- covariance_fit$ssn.object$obs
  partition <- get_ssn_simulate_partition(observed, covariance_fit$partition_factor)[permutation]
  covariance_fit$partition_factor <- NULL
  covariance_fit$partition_xlev <- NULL
  neighbors <- rep(list(list(positions = integer(), covariance = numeric())), length(permutation))
  groups <- split(seq_along(permutation), partition)
  for (group in groups) {
    if (length(group) < 2L) next
    for (j in seq.int(2L, length(group))) {
      position <- group[j]
      pool <- group[seq_len(j - 1L)]
      covariance <- get_decorrelate_observed_cross_covariance(covariance_fit,
        observed[permutation[position], , drop = FALSE], observed[permutation[pool], , drop = FALSE])
      if (length(pool) > size) {
        keep <- order(-abs(covariance), pool)[seq_len(size)]
        pool <- pool[keep]
        covariance <- covariance[keep]
      }
      neighbors[[position]] <- list(positions = pool, covariance = covariance)
    }
  }
  neighbors
}

draw_ssn_simulate_random <- function(design, samples) {
  q <- NCOL(design$observed)
  if (!q) return(NULL)
  effects <- matrix(rnorm(q * samples), q, samples) * sqrt(design$variance)
  unname(as.matrix(design$observed %*% effects))
}

is_ssn_simulate_independent <- function(params_object) {
  all(vapply(params_object[c("tailup", "taildown", "euclid")],
    function(x) isTRUE(x[["de"]] == 0), logical(1)))
}

validate_ssn_simulate_independent <- function(params_object, observed, additive) {
  for (component in c("tailup", "taildown", "euclid")) {
    params <- params_object[[component]]
    type <- remove_covtype(class(params))
    if (type == "none") next
    if (component == "euclid") {
      check_euclid_type(type)
      extra <- if (euclid_has_extra(type)) params[["extra"]] else NULL
      check_euclid_extra_parameters(type, params[["de"]], params[["range"]],
        extra, params[["rotate"]], params[["scale"]])
    } else {
      if (component == "tailup") check_tailup_type(type) else check_taildown_type(type)
      check_tailup_taildown_parameters(params[["de"]], params[["range"]])
    }
  }
  if (!inherits(params_object$tailup, "tailup_none")) {
    values <- suppressWarnings(as.numeric(observed[[additive]]))
    if (length(values) != NROW(observed) || any(!is.finite(values))) {
      stop("additive must name a column of finite values.", call. = FALSE)
    }
  }
  invisible(NULL)
}
