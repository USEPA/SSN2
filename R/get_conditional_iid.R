is_conditional_iid <- function(object) {
  is_iid_covariance(coef(object, "tailup"), coef(object, "taildown"),
    coef(object, "euclid"), object$random)
}

get_conditional_iid_variance <- function(object, simulate_covparams = FALSE) {
  if (isTRUE(simulate_covparams) || !is_conditional_iid(object)) return(NULL)
  nugget <- coef(object, "nugget")
  variance <- max(if (inherits(nugget, "nugget_none")) 0 else nugget[["nugget"]],
                  object$diagtol)
  if (!is.finite(variance) || variance < 0) stop("Invalid nugget variance.", call. = FALSE)
  if (variance == 0) return(NULL)
  variance
}

draw_conditional_iid_error <- function(n, samples, variance, chunk_size = 1000L) {
  value <- matrix(NA_real_, n, samples)
  for (start in seq.int(1L, samples, by = chunk_size)) {
    cols <- seq.int(start, min(samples, start + chunk_size - 1L))
    value[, cols] <- sqrt(variance) * matrix(rnorm(n * length(cols)), n, length(cols))
  }
  value
}

get_conditional_iid_glm <- function(object, Xnew, variance, samples, latent = FALSE, chunk_size = 1000L, center = get_latent_center(object)) {
  b <- center$betahat
  delta <- t(chol(vcov(object))) %*% matrix(rnorm(length(b) * samples), length(b), samples)
  beta <- sweep(delta, 1L, b, "+")
  rownames(beta) <- names(b)
  observed <- NULL
  if (latent) {
    X <- model.matrix(object)
    eta <- fitted(object, type = "link")
    D <- diag(get_D(object$family, eta, object$y, object$size,
      as.vector(coef(object, type = "dispersion"))))
    precision <- 1 / variance - D
    if (any(!is.finite(precision) | precision <= 0)) {
      stop("The conditional GLM latent covariance is not positive definite; check the fitted model and its convergence.", call. = FALSE)
    }
    q <- 1 / precision
    A <- (q / variance) * X
    offset <- model.offset(model.frame(object))
    w <- if (is.null(offset)) center$eta else center$eta - offset
    observed <- matrix(NA_real_, NROW(X), samples)
    for (start in seq.int(1L, samples, by = chunk_size)) {
      cols <- seq.int(start, min(samples, start + chunk_size - 1L))
      observed[, cols] <- sweep(A %*% delta[, cols, drop = FALSE] +
        sqrt(q) * matrix(rnorm(NROW(X) * length(cols)), NROW(X), length(cols)), 1L, w, "+")
    }
  }
  value <- draw_conditional_iid_error(NROW(Xnew), samples, variance, chunk_size)
  for (start in seq.int(1L, samples, by = chunk_size)) {
    cols <- seq.int(start, min(samples, start + chunk_size - 1L))
    value[, cols] <- Xnew %*% beta[, cols, drop = FALSE] + value[, cols, drop = FALSE]
  }
  rownames(value) <- rownames(Xnew)
  list(beta = beta, newdata = value, latent = observed)
}

draw_conditional_iid <- function(object, newdata_name, variance, samples, output,
                                 type = "link", newdata_size = NULL,
                                 center = if (inherits(object, "ssn_glm")) get_latent_center(object) else NULL) {
  glm <- inherits(object, "ssn_glm")
  context <- if (glm) get_conditional_context_glm(object, newdata_name, local = TRUE, center = center) else
    get_conditional_context(object, newdata_name, local = TRUE)
  if (!is.null(object$partition_factor) && NROW(context$newdata)) {
    model.frame(object$partition_factor, context$newdata)
  }
  val <- list()
  if ("object" %in% output) {
    val$object <- matrix(if (glm) context$w else context$y_object,
                         nrow = NROW(context$Xmat), ncol = samples)
  }
  if (!any(c("newdata", "beta", "latent") %in% output)) return(val)
  if (glm) {
    draws <- get_conditional_iid_glm(object, context$x0, variance, samples,
      latent = "latent" %in% output, center = center)
    beta <- draws$beta
    value <- draws$newdata
    if ("latent" %in% output) {
      val$latent <- draws$latent
      if (!is.null(context$offset)) val$latent <- sweep(val$latent, 1L, context$offset, "+")
      rownames(val$latent) <- rownames(context$Xmat)
    }
  } else {
    b <- context$betahat
    beta <- sweep(t(chol(context$cov_betahat)) %*% matrix(rnorm(length(b) * samples), length(b), samples),
                  1L, b, "+")
    rownames(beta) <- names(b)
    if ("newdata" %in% output) {
      value <- draw_conditional_iid_error(NROW(context$x0), samples, variance)
      for (start in seq.int(1L, samples, by = 1000L)) {
        cols <- seq.int(start, min(samples, start + 999L))
        value[, cols] <- context$x0 %*% beta[, cols, drop = FALSE] + value[, cols, drop = FALSE]
      }
      rownames(value) <- rownames(context$x0)
    }
  }
  if ("beta" %in% output) val$beta <- beta
  if ("newdata" %in% output) {
    if (!is.null(context$newdata_offset)) value <- sweep(value, 1L, context$newdata_offset, "+")
    if (glm && type != "link") {
      if (identical(context$family, "binomial") && is.null(newdata_size)) newdata_size <- rep(1, NROW(value))
      mu <- invlink(value, context$family, size = 1)
      value <- if (type == "response") {
        if (identical(context$family, "binomial")) mu * newdata_size else mu
      } else draw_glm_response(context$family, mu, context$dispersion, newdata_size)
    }
    if (context$add_newdata_rows) rownames(value) <- object$missing_index
    val$newdata <- value
  }
  val
}
