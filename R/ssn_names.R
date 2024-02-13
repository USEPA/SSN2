#' Return names of data in an SSN object
#'
#' @description Extract and print names from the \code{edges}, \code{sites}
#' and \code{preds} elements of an SSN object.
#'
#' @param ssn.object An \code{SSN} object.
#'
#' @return Print variable names to console
#'
#' @name ssn_names
#' @export
ssn_names <- function(ssn.object) {
  d <- ssn.object$obs
  no <- names(d)
  np <- length(ssn.object$preds)
  namesList <- vector("list", np + 1)
  namesList[[1]] <- no
  names4List <- "obs"
  if (np > 0) {
    for (i in seq_len(np)) {
      names4List <- c(names4List, names(ssn.object$preds)[[i]])
      d <- ssn.object$preds[[i]]
      namesList[[i + 1]] <- names(d)
    }
  }
  names(namesList) <- names4List
  return(namesList)
}
