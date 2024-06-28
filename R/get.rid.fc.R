#' Get reach ids (in C)
#'
#' @param binIDs Binary ID values in the stream network
#' @param referenceBinID A reference set of Binary ID values in the stream network
#'
#' @return The reach ids are unique ids for each stream segment in the stream network.
#' @noRd
get.rid.fc <- function(binIDs, referenceBinID) {
  ind.match <- .Call("test_fc", binIDs, referenceBinID)
  data.frame(
    fc = ind.match < 0,
    binaryID = substr(binIDs, 1, abs(ind.match)),
    stringsAsFactors = FALSE
  )
}
