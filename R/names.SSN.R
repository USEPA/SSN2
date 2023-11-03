#' names SSN object
#'
#' @description Extract and print names from the SSN object
#'
#' @param x An SSN object.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return Print variable names to console
#'
#' @name names.SSN
#' @method names SSN
#' @export


names.SSN <- function(x, ...) {
  d <- x$obs
  no <- names(d)
  np <- length(x$preds)
  namesList <- vector("list", np + 1)
  namesList[[1]] <- no
  names4List <- "obs"
  if (np > 0) {
    for (i in seq_len(np)) {
      names4List <- c(names4List, names(x$preds)[[i]])
      d <- x$preds[[i]]
      namesList[[i + 1]] <- names(d)
    }
  }
  names(namesList) <- names4List
  return(namesList)
}
