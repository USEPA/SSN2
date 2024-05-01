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

  if (length(ssn.object$obs) == 0) {
    no <- 0
    nameso <- NULL
  } else {
    # is obs present
    no <- 1
    nameso <- names(ssn.object$obs)
  }
  np <- length(ssn.object$preds)
  nt <- no + np
  if (nt == 0) {
    return(cat("There are no observed or prediction data and hence, no variables names."))
  }

  namesList <- vector("list", nt)
  if (no == 1) {
    namesList[[1]] <- nameso
    names4List <- "obs"
  } else {
    names4List <- NULL
  }
  if (np > 0) {
    prednames <- names(ssn.object$preds)
    for (i in seq_len(np)) {
      namespi <- names(ssn.object$preds[[i]])
      names4List <- c(names4List, prednames[i])
      namesList[[i + no]] <- namespi
    }
  }
  names(namesList) <- names4List
  print(namesList)
}
