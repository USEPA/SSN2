#' Import prediction points into an SSN, ssn_lm, or ssn_glm object
#'
#' @description A shapefile of prediction points found in the .ssn
#'   directory are imported into an existing object of class
#'   \code{SSN}, \code{ssn_lm}, or \code{ssn_glm}.
#' @param x An object of class\code{SSN}, \code{ssn_lm}, or
#'   \code{ssn_glm}.
#' @param predpts Name of the prediction point shapefile to import in
#'   character format, without the .shp extension.
#' @param format_additive Logical indicating whether the columns containing
#'   the addtive function values should be formated for
#'   \code{SSN2}. Default = \code{FALSE}.
#' @param names_additive Character vector of column names in observed and
#'   prediction site datasets containing additive function
#'   values. Must be defined if \code{format_additive = TRUE}. Default =
#'   \code{NULL}.
#'
#' @details \command{ssn_import_predpts} imports a shapefile of
#'   prediction points residing in the .ssn directory into an existing
#'   \code{SSN}, \code{ssn_lm}, or \code{ssn_glm} object. The
#'   prediction dataset must reside in the ssn.object$path
#'   directory. The path for an \code{SSN} object can be updated using
#'   \command{ssn_update_path()} prior to importing prediction
#'   datasets. Note that, the prediction dataset must contain the
#'   spatial, topological and attribute information needed to make
#'   predictions using an ssn_lm or ssn_glm object.  This information
#'   can be generated using a number of proprietary and open source
#'   software tools: \itemize{ \item{The Spatial Tools for the
#'   Analysis of River Systems (STARS) tools for ArcGIS Desktop
#'   versions 9.3x-10.8x (Peterson and Ver Hoef 2014). This custom
#'   ArcGIS toolset is designed to work with existing streams data in
#'   vector format.}  \item{The openSTARS package (Kattwinkel et
#'   al. 2020) extends the functionality of the STARS toolset, which
#'   makes use of R and GRASS GIS. It is open source and designed to
#'   derive streams in raster format from a digital elevation model
#'   (DEM).}  \item{The SSNbler package (currently in development as
#'   of September 2023) is an open source version of the STARS
#'   toolset, which makes use of the functionality found in the sf
#'   package to process streams data in vector format.}  }
#' @return an object of class \code{SSN}, \code{ssn_lm}, or
#'   \code{ssn_glm} which contains the new prediction dataset. The
#'   name of the prediction dataset in the preds list corresponds to
#'   the basenames of the prediction site shapefiles (without the .shp
#'   extension) specified in \code{predpts}. See
#'   \code{\link[SSN2]{ssn_import}} for a detailed description of
#'   the prediction dataset format within the \code{SSN} class object.
#'
#' @references
#' Kattwinkel, M., Szocs, E., Peterson, E., and Schafer,
#'   R.B. (2020) Preparing GIS data for analysis of stream monitoring
#'   data: The R package openSTARS. \emph{PLOS One} \bold{15(9)},
#'   e0239237.
#' Peterson, E., and Ver Hoef, J.M. (2014) STARS: An
#'   ArcGIS toolset used to calculate the spatial information needed
#'   to fit spatial statistical stream network models to stream
#'   network data. \emph{Journal of Statistical Software}
#'   \bold{56(2)}, 1--17.
#' @export
#' @examples
#' ## Create local temporary copy of MiddleFork04.ssn found in
#' # SSN2/lsndata folder. Only necessary for this example.
#' copy_lsn_to_temp()
#'
#' ## Import SSN object with no prediction sites
#' mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
#'   overwrite = TRUE
#' )
#'
#' ## Import pred1km prediction dataset into SSN object
#' mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"))
#' mf04p <- ssn_import_predpts(mf04p, predpts = "pred1km")
#' names(mf04p$preds)
#'
#' ## Import pred1km prediction dataset into a ssn_glm object
#' ssn_gmod <- ssn_glm(Summer_mn ~ netID, mf04p,
#'   family = "Gamma",
#'   tailup_type = "exponential", additive = "afvArea"
#' )
#' ssn_gmod <- ssn_import_predpts(ssn_gmod, predpts = "CapeHorn")
#' names(ssn_gmod$ssn.object$preds)
ssn_import_predpts <- function(x, predpts,
                               format_additive = FALSE, names_additive = NULL) {
  obj.type <- class(x)

  old_wd <- getwd()
  on.exit(setwd(old_wd))

  ## For fitted model objects- check if predpts already exists
  if (obj.type %in% c("ssn_lm", "ssn_glm")) {
    setwd(x$ssn.object$path)
    count <- 0

    if (length(x$ssn.object$preds) > 0) {
      for (m in seq_len(length(x$ssn.object$preds))) {
        if (names(x$ssn.object$preds[m]) == predpts) {
          pred.num <- m
          count <- count + 1
        }
      }
    }

    if (count > 0) {
      stop("Fitted model object already contains predpoints named ", predpts)
    }
  }

  ## For SSN objects - check if predpts already exits
  if (obj.type == "SSN") {
    setwd(x$path)
    count <- 0

    if (length(x$preds) > 0) {
      for (m in seq_len(length(x$preds))) {
        if (names(x$preds[m]) == predpts) {
          pred.num <- m
          count <- count + 1
        }
      }
    }

    if (count > 0) {
      stop("SSN already contains a prediction dataset named ", predpts)
    }
  }

  if (file.exists(paste(predpts, ".shp", sep = "")) == 0) {
    stop(paste(predpts, ".shp data is missing from ", old_wd,
      " folder",
      sep = ""
    ))
  }
  ## Read in prediction point shapefile
  predpoints <- st_read(paste0(predpts, ".shp"), quiet = TRUE)

  ## Check geometry type
  if (sum(st_geometry_type(predpoints, by_geometry = TRUE) == "POINT") != nrow(predpoints)) {
    stop(paste0(predpts, " does not have POINT geometry"))
  }

  if (format_additive == TRUE) {
    if (is.null(names_additive)) {
      stop("names_additive is required when format_additive = TRUE")
    }
    for (q in seq_len(length(names_additive))) {
      if (!names_additive[q] %in% colnames(predpoints)) {
        warning(paste0(names_additive[q], " is not found in ", predpts), call. = FALSE)
      } else {
        ## Extract afv column
        tmp.col <- st_drop_geometry(predpoints[, names_additive[q]])
        ## convert to character if column is numeric
        if (!inherits(tmp.col[, 1], "numeric")) {
          message(paste0(names_additive[q], " is not numeric. AFV values were not modified"))
        } else {
          tmp.col2 <- formatC(tmp.col[, 1],
            digits = 254, format = "f",
            drop0trailing = TRUE
          )
          predpoints[, names_additive[q]] <- tmp.col2
        }
      }
    }
  }

  ## Add netgeometry column
  predpoints[, "netgeometry"] <- paste0("SNETWORK (", paste(
    predpoints$netID,
    predpoints$rid,
    predpoints$upDist,
    predpoints$ratio,
    predpoints$pid,
    predpoints$locID
  ),
  ")",
  sep = ""
  )

  ## Put prediction points in SSN object
  if (obj.type == "SSN") {
    pred.num <- length(x$preds)
    x$preds[[pred.num + 1]] <- predpoints
    names(x$preds)[[pred.num + 1]] <- predpts
  }
  ## Put prediction points in fitted model object
  if (obj.type %in% c("ssn_glm", "ssn_lm")) {
    pred.num <- length(x$ssn.object$preds)
    x$ssn.object$preds[[pred.num + 1]] <- predpoints
    names(x$ssn.object$preds)[[pred.num + 1]] <- predpts
  }
  return(x)
}
