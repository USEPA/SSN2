skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")
skip_if_not_installed("pracma")

test_that("stream wave kernels agree with independent finite integrals", {
  sinc <- function(x) ifelse(x == 0, 1, sin(x) / x)
  si <- function(x) {
    if (x == 0) return(0)
    integrate(sinc, 0, x, subdivisions = 10000, rel.tol = 1e-11, abs.tol = 1e-12)$value
  }
  F <- function(x) {
    if (x == 0) return(0)
    integrate(function(u) 2 * sin(u / 2)^2 / u, 0, x,
              subdivisions = 10000, rel.tol = 1e-11, abs.tol = 1e-12)$value
  }
  h <- c(0, 1e-12, 1e-8, 1e-4, .1, .4999, .5, .5001, 1, 3.47, 5,
         12.3, 31.999, 32, 32.001, 40, 100)
  reference <- vapply(h, function(x) {
    if (x == 0) return(1)
    (cos(x) * F(2 * x) / x + sinc(x) * (pi - si(2 * x))) / pi
  }, numeric(1))
  expect_equal(stream_wave_connected(h, 1), reference, tolerance = 1e-10)
  pairs <- expand.grid(a = c(0, 1e-12, 1e-6, .009, .01, .5, 1, 3, 32, 100),
                       d = c(0, 1e-12, 1e-8, .000999, .001001, .1, 1, 10))
  pairs$b <- pairs$a + pairs$d
  reference <- mapply(function(a, b) {
    d <- b - a
    average <- if (d == 0) 2 * a * sinc(a)^2 else {
      integrate(function(t) 2 * (a + d * t) * sinc(a + d * t)^2, 0, 1,
                subdivisions = 10000, rel.tol = 1e-11, abs.tol = 1e-12)$value
    }
    (cos(d) * average + sinc(d) * (pi - si(2 * a) - si(2 * b))) / pi
  }, pairs$a, pairs$b)
  expect_equal(stream_wave_unconnected(pairs$a, pairs$b, 1), reference, tolerance = 1e-10)
  for (scale in c(1e-8, 1e8)) {
    expect_equal(stream_wave_connected(h * scale, scale), stream_wave_connected(h, 1), tolerance = 1e-12)
    expect_equal(stream_wave_unconnected(pairs$a * scale, pairs$b * scale, scale),
                 stream_wave_unconnected(pairs$a, pairs$b, 1), tolerance = 1e-12)
  }
  z <- c(63.999, 64, 64.001, 100, 1000, 1e6, 1e10, 1e12)
  tail <- vapply(z, function(x) {
    f <- integrate(function(u) exp(-u) / (1 + (u / x)^2), 0, 50, rel.tol = 1e-12)$value / x
    g <- integrate(function(u) u * exp(-u) / (1 + (u / x)^2), 0, 50, rel.tol = 1e-12)$value / x^2
    cos(x) * f + sin(x) * g
  }, numeric(1))
  expect_equal(stream_wave_si_tail(z), tail, tolerance = 1e-13)
  expect_true(all(is.finite(stream_wave_connected(c(1e-300, 1e300), 1e-10))))
  expect_equal(stream_wave_unconnected(c(0, 1), c(1, 1), 1e-320), c(0, 0))
  large <- c(1e100, 1e160, 1e200, 1e300)
  equal <- stream_wave_unconnected(large, large, 1)
  expect_true(all(equal > 0))
  expect_equal(equal * large, rep(1 / pi, length(large)), tolerance = 1e-12)
  for (distance in large) {
    params <- taildown_params("wave", de = distance, range = 1)
    observed <- list(a_mat = distance, b_mat = distance, hydro_mat = 2 * distance, mask_mat = 1)
    expect_equal(as.numeric(cov_matrix(params, observed)), 1 / pi, tolerance = 1e-12)
  }
})

test_that("stream wave overlap integrals have behaved truncation errors", {
  sinc <- function(x) ifelse(x == 0, 1, sin(x) / x)
  pairs <- data.frame(a = c(0, 0, 0, .2, 1, 1, 2), b = c(0, .5, 3.47, 4, 1, 1.0000001, 12))
  upper <- 20000
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$a[i]
    b <- pairs$b[i]
    d <- b - a
    boundaries <- c(a, seq(a + 100 * pi, upper, by = 100 * pi), upper)
    value <- sum(vapply(seq_len(length(boundaries) - 1), function(j) {
      integrate(function(u) sinc(u) * sinc(u + d), boundaries[j], boundaries[j + 1],
                subdivisions = 1000, rel.tol = 1e-9, abs.tol = 1e-12)$value
    }, numeric(1)))
    tail <- if (d == 0) 1 / (2 * upper) else cos(d) / (2 * d) * log1p(d / upper)
    reference <- 2 / pi * (value + tail)
    bound <- 1 / (pi * upper * (upper + d))
    expect_lt(abs(stream_wave_unconnected(a, b, 1) - reference), bound + 1e-10)
  }
})

test_that("stream wave matrices match a two-junction moving-average Gram-matrix", {
  parent <- c(A = "C", B = "C", C = "out", D = "out")
  afv <- c(A = 1, B = 3, C = 4, D = 2, out = 6)
  length_C <- 2
  sites <- data.frame(segment = rep(names(afv), each = 2),
                      position = rep(c(.4, 1.2), length(afv)))
  sites <- rbind(sites, sites[1, ])
  n <- nrow(sites)
  downstream <- function(from, x, to, y) {
    if (from == to) {
      d <- if (from == "out") y - x else x - y
      return(ifelse(d >= 0, d, Inf))
    }
    if (from == "out") return(rep(Inf, length(x + y)))
    d <- x
    segment <- parent[[from]]
    while (segment != to && segment != "out") {
      d <- d + length_C
      segment <- parent[[segment]]
    }
    if (segment != to) return(rep(Inf, length(x + y)))
    d + if (to == "out") y else length_C - y
  }
  a <- b <- h <- w <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    si <- sites$segment[i]
    sj <- sites$segment[j]
    x <- sites$position[i]
    y <- sites$position[j]
    d <- min(downstream(si, x, sj, y), downstream(sj, y, si, x))
    if (is.finite(d)) {
      h[i, j] <- b[i, j] <- d
      w[i, j] <- sqrt(min(afv[c(si, sj)]) / max(afv[c(si, sj)]))
    } else {
      junction <- if (all(c(si, sj) %in% c("A", "B"))) "C" else "out"
      position <- if (junction == "C") length_C else 0
      arms <- c(downstream(si, x, junction, position), downstream(sj, y, junction, position))
      a[i, j] <- min(arms)
      b[i, j] <- max(arms)
      h[i, j] <- sum(arms)
    }
  }
  ids <- c(19, 2, 41, 7, 33, 11, 59, 4, 23, 6, 61)
  dimnames(h) <- dimnames(a) <- dimnames(b) <- dimnames(w) <- list(as.character(ids), as.character(ids))
  observed <- list(hydro_mat = h, a_mat = a, b_mat = b, mask_mat = matrix(1, n, n), w_mat = w)
  covariance <- lapply(c("tailup", "taildown"), function(type) {
    as.matrix(cov_matrix(get(paste0(type, "_params"))("wave", de = 1, range = 1), observed))
  })
  nodes <- 12L
  j <- seq_len(nodes - 1L)
  jacobi <- matrix(0, nodes, nodes)
  jacobi[cbind(j, j + 1L)] <- jacobi[cbind(j + 1L, j)] <- j / sqrt(4 * j^2 - 1)
  eig <- eigen(jacobi, symmetric = TRUE)
  upper <- 1024 * pi
  gram <- list(matrix(0, n, n), matrix(0, n, n))
  sinc <- function(x) { value <- numeric(length(x)); ok <- is.finite(x); value[ok] <- ifelse(x[ok] == 0, 1, sin(x[ok]) / x[ok]); value }
  for (segment in names(afv)) {
    limit <- if (segment == "C") length_C else upper
    breaks <- sort(unique(c(0, sites$position, seq(0, limit, length.out = ceiling(limit / pi) + 1), limit)))
    mid <- (head(breaks, -1) + tail(breaks, -1)) / 2
    half <- diff(breaks) / 2
    x <- as.vector(outer(half, eig$values) + mid)
    weight <- as.vector(outer(half, 2 * eig$vectors[1, ]^2))
    ku <- kd <- matrix(0, length(x), n)
    for (i in seq_len(n)) {
      ku[, i] <- sinc(downstream(segment, x, sites$segment[i], sites$position[i])) *
        sqrt(afv[[segment]] / afv[[sites$segment[i]]])
      kd[, i] <- sinc(downstream(sites$segment[i], sites$position[i], segment, x))
    }
    gram[[1]] <- gram[[1]] + crossprod(ku * sqrt(weight)) * 2 / pi
    gram[[2]] <- gram[[2]] + crossprod(kd * sqrt(weight)) * 2 / pi
  }
  bound <- 2 / (pi * (upper - length_C - max(sites$position)))
  for (k in 1:2) {
    type <- c("tailup", "taildown")[k]
    sigma <- covariance[[k]]
    expect_lt(max(abs(sigma - gram[[k]])), bound)
    expect_equal(diag(sigma), rep(1, n), ignore_attr = TRUE)
    expect_equal(sigma, t(sigma))
    expect_true(all(abs(sigma) <= 1 + 1e-12))
    expect_gt(min(eigen(sigma, symmetric = TRUE)$values), -n * norm(sigma, "2") * 1e-12)
    params <- get(paste0(type, "_params"))("wave", de = 1, range = 1)
    for (rows in list(c(8, 2, 5), 3L, integer())) {
      cols <- c(10, 6, 1, 4)
      cross <- lapply(observed, function(x) Matrix::Matrix(x[rows, cols, drop = FALSE]))
      names(cross) <- sub("_mat$", "_pred_mat", names(cross))
      expect_equal(as.matrix(cov_vector(params, cross)), sigma[rows, cols, drop = FALSE])
    }
  }
  expect_true(all(covariance[[1]][1:2, 3:4] == 0))
  expect_equal(covariance[[1]][1, 5], covariance[[2]][1, 5] * sqrt(1 / 4))
  expect_equal(covariance[[1]][1, 9], covariance[[2]][1, 9] * sqrt(1 / 6))
  expect_lt(covariance[[2]][1, 3], 1)
})

test_that("stream wave sparse and dense grid anchors use the kernel scale", {
  local_mocked_bindings(get_effective_range = function(...) stop("unexpected effective range"))
  for (dense in c(FALSE, TRUE)) {
    wave <- ssn_decorrelate_grid(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "wave",
      taildown_type = "wave", euclid_type = "wave", additive = "afvArea", dense_grid = dense)
    reference <- ssn_decorrelate_grid(Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "spherical",
      taildown_type = "spherical", euclid_type = "wave", additive = "afvArea", dense_grid = dense)
    for (component in c("tailup", "taildown", "euclid")) {
      expect_equal(wave[[paste0(component, "_range")]], reference[[paste0(component, "_range")]])
    }
    expect_equal(nrow(wave), nrow(reference))
  }
})
