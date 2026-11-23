#' @rdname print.SSN2
#' @method print ssn_lmRF
#' @export
print.ssn_lmRF <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Random forest mean model:\n")
  print(x$ranger, ...)
  cat("\nSSN residual model:\n")
  print(x$ssn_lm, digits = digits, ...)
  invisible(x)
}

#' @rdname print.SSN2
#' @method print summary.ssn_lmRF
#' @export
print.summary.ssn_lmRF <- function(x, digits = max(3L, getOption("digits") - 3L),
                                   signif.stars = getOption("show.signif.stars"), ...) {
  cat("Random forest mean model:\n")
  print(x$ranger, ...)
  cat("\nSSN residual model:\n")
  print(x$ssn_lm, digits = digits, signif.stars = signif.stars, ...)
  invisible(x)
}
