#' @title Copy LSN to temporary directory
#'
#' @description Copies the LSN directory MiddleFork04.ssn to R's temporary
#'   directory so the examples in SSN2 do not write to the local
#'   library or any other places.
#' @return A copy of MiddleFork04.ssn residing in R's temporary directory
#'
#' @details Copies the LSN directory MiddleFork04.ssn to R's temporary directory
#' @export
#' @examples
#' copy_lsn_to_temp()
#' # getwd()
#' # setwd(tempdir())
#' # getwd()
#' # if unix-alike, list temporary directory contents using: system('ls')
#' # if windows, list temporary directory contents using: shell('dir')
copy_lsn_to_temp <- function() {
  ## Create temporary .ssn directory to work with
  if (dir.exists(paste0(tempdir(), "/MiddleFork04.ssn"))) {
    return(invisible())
  }
  file.copy(
    from = system.file("lsndata/MiddleFork04.ssn", package = "SSN2"),
    to = tempdir(), recursive = TRUE,
    overwrite = FALSE, copy.mode = FALSE
  )
  if (.Platform$OS.type == "unix") {
    system(paste0("chmod -R 777 ", tempdir(), "/MiddleFork04.ssn"))
  }
  invisible()
}
