#' Update path in an SSN object
#'
#' @description Update the local path in an existing SSN object based
#'   on an user defined file.
#' @param x An SSN, ssn_lm or ssn_glm object.
#' @param path Filepath to the .ssn folder associated with the SSN
#'   object.
#' @param verbose A logical that indicates if the new path should be printed
#'   to the console.
#'
#' @details At times, it may be necessary to move a .ssn directory,
#'   which is linked to an SSN object in an R workspace. If the .ssn
#'   directory is moved, the path must be updated before using the
#'   \code{ssn_glmssn} function and other functions that read/write
#'   to the .ssn directory. The \command{ssn_update_path} is a helper
#'   function that serves this purpose.
#'
#' @return An SSN object with a new path list element.
#'
#' @export
#' @examples
#' ## Use mf04p SSN object provided in SSN2
#' data(mf04p)
#'
#' ## For examples only, make sure mf04p has the correct path
#' ## If you use ssn_import(), the path will be correct
#' newpath <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_update_path(mf04p, newpath)
ssn_update_path <- function(x, path, verbose = FALSE) {
  file <- path

  if (inherits(x, "SSN")) {
    x$path <- file
  }

  if (inherits(x, c("ssn_lm", "ssn_glm"))) {
    x$ssn.object$path <- file
  }

  if (verbose) {
    message(paste("SSN path updated to", file, sep = " "))
  }

  return(x)
}
