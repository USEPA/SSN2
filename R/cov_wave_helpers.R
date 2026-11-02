stream_wave_dependency <- function() {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Package 'pracma' is required for tailup and taildown wave covariance. Install it with install.packages('pracma').",
         call. = FALSE)
  }
}

stream_wave_sinc <- function(x) {
  out <- rep(1, length(x))
  nonzero <- x != 0
  out[nonzero] <- sin(x[nonzero]) / x[nonzero]
  out
}

stream_wave_F <- function(z) {
  original <- z
  z <- unique(z)
  out <- numeric(length(z))
  small <- z <= 1
  if (any(small)) {
    zz <- z[small]^2
    term <- zz / 4
    value <- term
    for (k in 1:8) {
      term <- -term * zz * (2 * k) / ((2 * k + 2)^2 * (2 * k + 1))
      value <- value + term
    }
    out[small] <- value
  }
  if (any(!small)) {
    zz <- z[!small]
    out[!small] <- -digamma(1) + log(zz) - pracma::Ci(zz)
  }
  out[match(original, z)]
}

stream_wave_si_tail <- function(z) {
  original <- z
  z <- unique(z)
  out <- numeric(length(z))
  large <- z >= 64
  if (any(!large)) out[!large] <- pi / 2 - pracma::Si(z[!large])
  if (any(large)) {
    zz <- z[large]
    inverse <- 1 / zz
    f_term <- inverse
    g_term <- inverse^2
    f <- f_term
    g <- g_term
    for (k in 1:10) {
      f_term <- -f_term * (2 * k) * (2 * k - 1) * inverse^2
      g_term <- -g_term * (2 * k + 1) * (2 * k) * inverse^2
      f <- f + f_term
      g <- g + g_term
    }
    out[large] <- cos(zz) * f + sin(zz) * g
  }
  out[match(original, z)]
}

stream_wave_ci <- function(z) {
  values <- unique(z)
  if (!length(values)) return(numeric())
  pracma::Ci(values)[match(z, values)]
}

stream_wave_connected <- function(h, range) {
  stream_wave_dependency()
  stream_wave_check(h, range)
  original <- as.numeric(h) / range
  x <- unique(original)
  out <- numeric(length(x))
  out[x == 0] <- 1
  active <- x > 0 & is.finite(2 * x)
  if (any(active)) {
    xx <- x[active]
    out[active] <- cos(xx) * stream_wave_F(2 * xx) / (pi * xx) +
      stream_wave_sinc(xx) * (0.5 + stream_wave_si_tail(2 * xx) / pi)
  }
  out[match(original, x)]
}

stream_wave_unconnected <- function(a, b, range) {
  stream_wave_dependency()
  stream_wave_check(a, range)
  stream_wave_check(b, range)
  if (length(a) != length(b) || any(a > b)) {
    stop("Stream wave junction distances must have equal lengths and satisfy a <= b.", call. = FALSE)
  }
  aa <- as.numeric(a) / range
  bb <- as.numeric(b) / range
  d <- (as.numeric(b) - as.numeric(a)) / range
  out <- numeric(length(aa))
  active <- is.finite(2 * bb) & is.finite(d)
  aa <- aa[active]
  bb <- bb[active]
  dd <- d[active]
  divided <- numeric(length(dd))
  near <- dd <= 0.001
  if (any(near)) {
    m <- aa[near] + dd[near] / 2
    g <- 2 * sin(m) * stream_wave_sinc(m)
    g_second <- numeric(length(m))
    small <- m < 0.01
    mm <- m[small]
    g_second[small] <- -4 * mm + 16 / 9 * mm^3 - 4 / 15 * mm^5
    mm <- m[!small]
    g_second[!small] <- 4 * cos(2 * mm) / mm - 4 * sin(2 * mm) / mm^2 +
      4 * sin(mm)^2 / mm^3
    divided[near] <- g + dd[near]^2 * g_second / 24
  }
  regular <- !near
  if (any(regular)) {
    ar <- aa[regular]
    br <- bb[regular]
    dr <- dd[regular]
    difference <- numeric(length(ar))
    small <- ar < 0.5
    difference[small] <- stream_wave_F(2 * br[small]) - stream_wave_F(2 * ar[small])
    if (any(!small)) {
      difference[!small] <- log1p(dr[!small] / ar[!small]) +
        stream_wave_ci(2 * ar[!small]) - stream_wave_ci(2 * br[!small])
    }
    divided[regular] <- difference / dr
  }
  out[active] <- (cos(dd) * divided + stream_wave_sinc(dd) *
    (stream_wave_si_tail(2 * aa) + stream_wave_si_tail(2 * bb))) / pi
  out
}

stream_wave_check <- function(distance, range) {
  if (!is.numeric(range) || length(range) != 1L || !is.finite(range) || range <= 0) {
    stop("The stream wave range must be finite and positive.", call. = FALSE)
  }
  if (any(!is.finite(distance) | distance < 0)) {
    stop("Stream wave distances must be finite and nonnegative.", call. = FALSE)
  }
}

stream_wave_covariance <- function(params, h, m, w = NULL, a = NULL, b = NULL) {
  stream_wave_dependency()
  value <- h
  value[] <- 0
  if (!is.null(w)) {
    active <- which(m != 0 & w != 0)
    value[active] <- params[["de"]] * w[active] *
      stream_wave_connected(h[active], params[["range"]])
  } else {
    active <- which(m != 0)
    stream_wave_check(a[active], params[["range"]])
    stream_wave_check(b[active], params[["range"]])
    connected <- active[a[active] == 0]
    unconnected <- active[a[active] != 0]
    value[connected] <- params[["de"]] *
      stream_wave_connected(h[connected], params[["range"]])
    value[unconnected] <- params[["de"]] *
      stream_wave_unconnected(a[unconnected], b[unconnected], params[["range"]])
  }
  value
}
