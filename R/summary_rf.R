#' @rdname summary.SSN2
#' @method summary ssn_lmRF
#' @export
summary.ssn_lmRF <- function(object, ...) {
  structure(
    list(ranger = object$ranger, ssn_lm = summary(object$ssn_lm)),
    class = "summary.ssn_lmRF"
  )
}
