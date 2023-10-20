restruct_ssn_missing <- function(ssn.object, observed_index, missing_index) {
  if (length(missing_index) > 0) {
    # only return .missing object if necessary
    ssn.object$preds$.missing <- ssn.object$obs[missing_index, , drop = FALSE]
  }
  ssn.object$obs <- ssn.object$obs[observed_index, , drop = FALSE]
  ssn.object
}
