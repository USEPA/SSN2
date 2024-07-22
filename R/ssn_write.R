#' write an SSN object
#'
#' @description This function writes an \code{SSN} object to a local
#'   .ssn directory
#'
#' @param ssn An \code{SSN} object.
#' @param path filepath to the local .ssn directory to write to.
#' @param overwrite If \code{TRUE}, overwrite existing files in \code{file}
#'   (if it exists). Defaults to \code{FALSE}.
#' @param copy_dist If \code{TRUE}, copy distance matrices to \code{file}
#'   (if they exist). Defaults to \code{FALSE}.
#' @param import If \code{TRUE}, import and return the \code{SSN} object
#'   after writing to file. Defaults to \code{FALSE}.
#'
#' @return{ssn_write} creates an .ssn directory that contains the
#'   spatial, topological, and attribute information stored in the
#'   original \code{SSN} object. Spatial datasets found in the
#'   \code{SSN} object (e.g. edges, obs, and prediction sites) are
#'   saved in geopackage format. When \code{import = TRUE}, the
#'   \code{SSN} object is imported and returned.
#'
#' @export
#' @examples
#' ## For examples only, copy MiddleFork04.ssn directory to R's
#' # temporary directory
#' copy_lsn_to_temp()
#' ## Import SSN object with prediction sites
#' mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
#'   predpts = "pred1km",
#'   overwrite = TRUE
#' )
#'
#' ## Write SSN to new .ssn directory
#' ssn_write(mf04p,
#'   path = paste0(tempdir(), "/tempSSN.ssn"),
#'   overwrite = TRUE
#' )
#'
#' ## Write SSN to .ssn directory and return SSN object
#' tempSSN <- ssn_write(mf04p, path = paste0(
#'   tempdir(),
#'   "/tempSSN.ssn"
#' ), overwrite = TRUE, import = TRUE)
ssn_write <- function(ssn, path, overwrite = FALSE,
                      copy_dist = FALSE, import = FALSE) {
  ## Add .ssn extension if necessary
  if (substr(path, nchar(path) - 3, nchar(path)) != ".ssn") {
    print(paste0(
      "path must include .ssn extension. ",
      path, " updated to ", paste0(path, ".ssn")
    ))
    path <- paste0(path, ".ssn")
  }

  suppressWarnings({
    if (!file.exists(path)) {
      dir.create(path)
    } else {
      if (overwrite == FALSE) stop(paste0(path, " exists and overwrite = FALSE"))
      if (overwrite == TRUE) {
        unlink(path, recursive = TRUE)
        dir.create(path)
      }
    }

    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(path)

    #######################################################################
    ## Get vector of filenames not associated with shapefiles or geopackages
    ssn.tmp <- ssn
    pred.len <- length(ssn.tmp$preds)

    ssn.tmp$path <- getwd()

    ssn.files <- list.files(ssn$path)

    ind.gpkg <- substr(ssn.files, nchar(ssn.files) - 4, nchar(ssn.files)) == ".gpkg"
    ind.shp <- substr(ssn.files, nchar(ssn.files) - 3, nchar(ssn.files)) == ".shp"

    sub.list <- substr(ssn.files[ind.gpkg], 1, nchar(ssn.files[ind.gpkg]) - 5)
    sub.list <- c(sub.list, substr(ssn.files[ind.shp], 1, nchar(ssn.files[ind.shp]) - 4))

    ind.list <- vector(mode = "logical", length = length(ssn.files))

    for (j in seq_len(length(ssn.files))) {
      tmp <- unlist(strsplit(ssn.files[j], "[.]"))
      if ((tmp[1] %in% sub.list) == TRUE) ind.list[j] <- TRUE
    }

    ssn.files <- ssn.files[!ind.list]

    ## Copy files to new .ssn directory
    for (i in seq_len(length(ssn.files))) {
      fn.old <- file.path(ssn$path, ssn.files[i])
      if (basename(fn.old) != "distance") {
        fn.new <- file.path(ssn.tmp$path, ssn.files[i])
        file.copy(fn.old, fn.new, overwrite = TRUE)
      }
    }
    rm(fn.old, fn.new)

    ## Copy distance matrices
    if (copy_dist == TRUE & ("distance" %in% ssn.files)) {
      file.copy(paste0(ssn$path, "/distance"), getwd(), recursive = TRUE)
    }

    ## Copy observed sites if they exist
    if (class(ssn.tmp$obs)[1] == c("sf")) {
      st_write(ssn$obs, paste0(ssn.tmp$path, "/sites.gpkg"), quiet = TRUE)
    }

    ## Copy edges
    st_write(ssn$edges, paste0(ssn.tmp$path, "/edges.gpkg"), quiet = TRUE)

    ## Copy prediction sites
    if (pred.len > 0) {
      ## Check prediction datasets for duplicate names
      if (sum(duplicated(attributes(ssn.tmp$preds)$names)) > 0) {
        stop("Duplicated prediction dataset names are not allowed.")
      }
      pred.name.vec <- attributes(ssn.tmp$preds)$names
      for (i in seq_len(pred.len)) {
        pred.name <- pred.name.vec[i]
        st_write(ssn$preds[[pred.name]], paste0(
          ssn.tmp$path, "/",
          pred.name, ".gpkg"
        ),
        quiet = TRUE
        )
        rm(pred.name)
      }
    }

    ## Import SSN without prediction sites
    if (import == TRUE & pred.len == 0) {
      ssn.tmp <- ssn_import(ssn.tmp$path, overwrite = FALSE)
      return(ssn.tmp2)
    }
    ## Import SSN with all prediction sites
    if (import == TRUE & pred.len > 0) {
      ssn.tmp2 <- ssn_import(ssn.tmp$path, overwrite = FALSE)
      for (j in seq_len(pred.len)) {
        ssn.tmp2 <- ssn_import_predpts(ssn.tmp2, pred.name.vec[j])
      }
      return(ssn.tmp2)
    }
  })
}
