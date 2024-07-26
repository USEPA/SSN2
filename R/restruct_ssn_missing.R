#' Restructure a SSN object with missing response data
#'
#' @param ssn.object SSN object
#' @param observed_index Index of observed (non-NA) response values
#' @param missing_index Index of missing (NA) response values
#'
#' @return An SSN object with missing data stored as a prediction data set and
#'   observed data adjusted accordingly
#'
#' @noRd
restruct_ssn_missing <- function(ssn.object, observed_index, missing_index) {
  if (length(missing_index) > 0) {
    # only return .missing object if necessary
    ssn.object$preds$.missing <- ssn.object$obs[missing_index, , drop = FALSE]
  }
  ssn.object$obs <- ssn.object$obs[observed_index, , drop = FALSE]
  ssn.object
}
