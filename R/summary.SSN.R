#' @title Summarize an SSN object
#'
#' @description Summarize data found in an SSN object.
#'
#' @param object An SSN object.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details \code{summary.SSN()} creates a summary of a SSN object
#'   intended to be printed using \code{print()}. This summary
#'   contains information about the number of observed and prediction
#'   locations, as well as the column names found in their respective
#'   sf data.frames.
#'
#' @return A list with several fitted model quantities used to create
#'   informative summaries when printing.
#'
#' @name summary.SSN
#' @method summary SSN
#' @export
summary.SSN <- function(object, ...) {
  print(object)
}
