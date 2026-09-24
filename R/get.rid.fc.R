#' Get flow connectivity and common downstream binary IDs
#'
#' @param binIDs Binary ID values in the stream network
#' @param referenceBinID A single reference binary ID in the stream network
#'
#' @return A data frame with flow connectivity and common downstream binary IDs.
#' @noRd
get.rid.fc <- function(binIDs, referenceBinID) {
  # this is where src was used
  # ind.match <- .Call("test_fc", binIDs, referenceBinID)
  ind.match <- get_binary_id_match(binIDs, referenceBinID)
  data.frame(
    fc = ind.match < 0,
    binaryID = substr(binIDs, 1, abs(ind.match)),
    stringsAsFactors = FALSE
  )
}
