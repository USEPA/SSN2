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
