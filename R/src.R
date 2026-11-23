# Native sources are archived in src_deprecated/ and excluded from package builds.

#' Find the longest common binary-ID prefix length against a reference
#'
#' Uses base-R binary search: for each ID, finds the length of
#' its longest common prefix with \code{reference}, negated when the shorter
#' of the two IDs is itself a complete prefix of the other (signaling flow
#' connectivity in \code{\link{get.rid.fc}()}).
#'
#' @param ids A character vector of binary stream IDs to compare.
#' @param reference A single reference binary stream ID.
#'
#' @return An integer vector (length \code{length(ids)}) of common-prefix
#'   lengths, negated when the shorter ID is a complete prefix of the longer.
#'
#' @noRd
get_binary_id_match <- function(ids, reference) {
  if (!is.character(ids) || !is.character(reference)) stop("invalid arguments")
  if (length(reference) != 1L) stop("reference must have length one")
  if (anyNA(ids) || is.na(reference)) stop("missing binary ID")
  n <- length(ids)
  shorter <- pmin.int(nchar(ids), nchar(reference))
  lower <- integer(n)
  upper <- shorter + 1L
  reference <- rep.int(reference, n)
  active <- which(upper - lower > 1L)
  while (length(active)) {
    middle <- (lower[active] + upper[active]) %/% 2L
    matches <- substr(ids[active], 1L, middle) ==
      substr(reference[active], 1L, middle)
    lower[active[matches]] <- middle[matches]
    upper[active[!matches]] <- middle[!matches]
    active <- active[upper[active] - lower[active] > 1L]
  }
  complete <- lower == shorter
  lower[complete] <- -shorter[complete]
  lower
}
