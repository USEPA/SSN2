#' @method AICc ssn_lm
#' @export
AICc.ssn_lm <- function(object, ..., k = 2) {
  # set k as 2
  k <- 2
  # store object and ...
  object_list <- list(object, ...)
  # see if ... has any elements
  if (length(object_list) == 1) {
    # number of estimated parameters
    if (object$estmethod == "ml") {
      n_est_param <- object$npar + object$p
    } else {
      n_est_param <- object$npar
    }
    # error if not ml or reml
    if (!object$estmethod %in% c("ml", "reml")) {
      stop("AICc is only defined if estmethod is \"ml\" or \"reml\".", call. = FALSE)
    }
    # compute AICc
    AICc_val <- -2 * as.numeric(logLik(object)) + 2 * object$n * (n_est_param) / (object$n - n_est_param - 1)
  } else {
    # warning if ml and reml in same call
    est_methods <- vapply(object_list, function(x) x$estmethod, character(1))
    if ("ml" %in% est_methods && "reml" %in% est_methods) {
      warning("AICc and AICcc should not compare models fit with
             \"ml\" to models fit with \"reml\"", call. = FALSE)
    }
    # warning if reml and fixed effects change
    est_methods_reml <- which(est_methods == "reml")
    if (length(est_methods_reml) > 1) {
      if (any(vapply(
        est_methods_reml,
        function(x) !identical(formula(object_list[[x]]), formula(object_list[[1]])), logical(1)
      ))) {
        warning("AIC and AICc should not be used to compare models fit with \"reml\" whose fixed effect formulas differ.", call. = FALSE)
      }
    }
    # find model names provided
    object_list_names <- as.character(c(substitute(object), (as.list(substitute(list(...)))[-1])))
    # error if any names duplicated
    if (any(duplicated(object_list_names))) {
      stop("Each model object must have a unique name", call. = FALSE)
    }
    # iterate through each model
    object_AICc <- lapply(object_list, function(x) {
      # warning if estmethod not ml or reml
      if (!object$estmethod %in% c("ml", "reml")) {
        stop("AICc is only defined is estmethod is \"ml\" or \"reml\".", call. = FALSE)
      }
      if (x$estmethod == "ml") {
        n_est_param <- x$npar + x$p
      } else {
        n_est_param <- x$npar
      }
      # store degrees of freedom (parames estimated) and AICc
      data.frame(df = n_est_param, AICc = -2 * logLik(x) + 2 * x$n * (n_est_param) / (x$n - n_est_param - 1))
    })
    # put all AICc data frames together
    AICc_val <- do.call("rbind", object_AICc)
    # set rownames as model names
    row.names(AICc_val) <- object_list_names
  }
  # return AICc value
  AICc_val
}

#' @method AICc ssn_glm
#' @export
AICc.ssn_glm <- AICc.ssn_lm
