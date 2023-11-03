#' Import \code{SSN} object
#'
#' @description This function reads spatial data from a .ssn folder
#'   and creates an \code{SSN} object.
#' @param path Filepath to the .ssn directory. See details.
#' @param include_obs default = \code{TRUE}. Logical indicating
#'   whether observed sites should be included in the SSN object.
#' @param predpts Vector of shapefile basenames for prediction sites
#'   found within the .ssn folder.
#' @param format_additive Logical indicating whether the columns containing
#'   the addtive function values should be formated for
#'   \code{SSN2}. Default = \code{FALSE}.
#' @param names_additive Character vector of column names in observed and
#'   prediction site datasets containing additive function values.
#'   Must be defined if \code{format_additive = TRUE}. Default =
#'   \code{NULL}.
#' @param overwrite default = \code{FALSE}. If \code{TRUE}, overwrite
#'   existing binaryID.db files.
#' @details The \command{importSSN} function imports spatial data from a .ssn
#'   folder to create an \code{SSN} object. The information contained in the
#'   .ssn folder can be generated using a number of proprietary and
#'   open source software tools:
#'   \itemize{
#'     \item{The Spatial Tools for the Analysis of River Systems
#'     (STARS) tools for ArcGIS Desktop versions 9.3x-10.8x (Peterson
#'     and Ver Hoef 2014). This custom ArcGIS toolset is designed to
#'     work with existing streams data in vector format.}
#'     \item{The openSTARS package (Kattwinkel et al. 2020) extends
#'     the functionality of the STARS toolset, which makes use of R
#'     and GRASS GIS. It is open source and designed to derive streams
#'     in raster format from a digital elevation model (DEM).}
#'     \item{The SSNbler package (currently in development as of
#'     September 2023) is an open source version of the STARS toolset,
#'     which makes use of the functionality found in the sf package to
#'     process streams data in vector format.}
#'   }
#'   When spatial data are processed using one of these software
#'   tools, a .ssn directory is output which contains all of the
#'   spatial, topological and attribute data needed to fit a spatial
#'   statistical stream network model to streams data.  This includes:
#'   \itemize{
#'     \item{An edges shapefile of lines that represent the stream
#'     network.}
#'     \item{A sites shapefile of points where observed data were
#'     collected on the stream network.}
#'     \item{Prediction sites shapefile(s) of locations where
#'     predictions will be made.}
#'     \item{netID.dat files for each distinct network, which store
#'     the topological relationships of the line segments in edges.}
#'   }
#'   A more detailed description of the .ssn directory and its
#'   contents is provided in Peterson and Ver Hoef (2014).
#'
#'   The \command{ssn_import} imports the edges, observed sites, and
#'   prediction sites as \code{sf data.frame} objects. A new column named 'netgeometry'
#'   is created to store important data that represents
#'   topological relationships in a spatial stream network
#'   model. These data are stored in character format, which is less
#'   likely to be inadvertantly changed by users. See
#'   \code{\link[SSN2]{ssn_get_netgeometry}} for a more detailed description of
#'   the format and contents of 'netgeometry'.
#'
#'   The information contained in the netID text files is imported
#'   into an SQLite database, binaryID.db, which is stored in the .ssn
#'   directory. This information is used internally by
#'   \code{\link[SSN2]{ssn_create_distmat}},
#'   \code{\link[SSN2]{ssn_lm}} and
#'   \code{\link[SSN2]{ssn_glm}} to calculate the data necessary
#'   to fit a spatial statistical model to stream network data. If
#'   \code{overwrite = TRUE} (\code{overwrite = FALSE} is the default) and a binaryID.db
#'   file already exists within the .ssn directory, it will be
#'   overwriten when the \code{SSN} object is created.
#'
#'   At a minimum, an \code{SSN} object must always contain streams, which
#'   are referred to as edges. The \code{SSN} object would also typically
#'   contain a set of observed sites, where measurements have been
#'   collected and only one observed dataset is permitted. When
#'   \code{include_obs=FALSE}, an \code{SSN} object is created without
#'   observations. This option provides flexibility for users who
#'   would like to simulate data on a set of artifical sites on an
#'   existing stream network. Note that observation sites must be
#'   included in the \code{SSN} object in order to fit models using
#'   \command{ssn_lm} or \command{ssn_glm}. The \code{SSN} object may contain
#'   multiple sets of prediction points (or none), which are stored as
#'   separate shapefiles in the .ssn directory. The
#'   \code{\link[SSN2]{ssn_import_predpts}} function allows users to import additional
#'   sets of prediction sites to a an existing \code{SSN} object.
#'
#' @return \code{ssn_import} returns an object of class SSN, which is a list
#'   with four elements containing:
#'   \itemize{
#'     \item{\code{edges}: An \code{sf data.frame} containing the stream network,
#'     with an additional 'netgeometry' column.}
#'     \item{\code{obs}: An sf data.frame containing observed site locations,
#'     with an additional 'netgeometry' column. NA if \code{include_obs =
#'     FALSE}.}
#'     \item{\code{preds}: A list of sf data.frames containing prediction
#'     site locations. The names of the preds list correspond to the
#'     basenames of the prediction site shapefiles (without the .shp
#'     extension) specified in \code{predpts}. Empty list if \code{predpts} is not provided.}
#'     \item{path: The local file to the .ssn directory associated with the \code{SSN}
#'     object.}
#'   }
#' @seealso [ssn_get_netgeometry]
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
#' mf04 <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
#'   overwrite = TRUE
#' )
#'
#' ## Import SSN object with 3 sets of prediction sites
#' mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
#'   predpts = c(
#'     "pred1km.shp",
#'     "CapeHorn.shp",
#'     "Knapp.shp"
#'   ),
#'   overwrite = TRUE
#' )
#'
ssn_import <- function(path, include_obs = TRUE, predpts,
                       format_additive = FALSE, names_additive = NULL,
                       overwrite = FALSE) {
  ## Extract dirname
  file <- path
  ssn_folder <- file

  # Get wd
  old_wd <- getwd()
  on.exit(setwd(old_wd)) # if the function crashes or finishes, it will restore the initial working directory

  ## Adjust if relative pathname is supplied and fix edges
  if (substr(ssn_folder, start = 1, stop = 2) == "./") {
    rel.path <- substr(ssn_folder, start = 2, stop = nchar(ssn_folder))
    ssn_folder <- paste0(old_wd, rel.path)
  }

  if (!dir.exists(ssn_folder)) stop("Cannot find the .ssn folder.")

  # Set wd
  setwd(ssn_folder)

  ## Check inputs, generally
  missingsites <- !include_obs
  missingpreds <- missing(predpts)

  # Check observation sites exist and import
  if (!missingsites) {
    obs_sites <- "sites.shp"

    # Check the file exists and then ingest, process
    if (file.exists(obs_sites)) {
      sfsites <- st_read(obs_sites, quiet = TRUE)

      ## Check geometry type
      if (sum(st_geometry_type(sfsites, by_geometry = TRUE) == "POINT") != nrow(sfsites)) {
        stop("Observed sites do not have POINT geometry")
      }
      ## Add network geometry column
      sfsites[, "netgeometry"] <- paste0("SNETWORK (", paste(
        sfsites$netID, sfsites$rid, sfsites$upDist,
        sfsites$ratio, sfsites$pid, sfsites$locID
      ),
      ")",
      sep = ""
      )

      if (format_additive == TRUE) {
        if (is.null(names_additive)) {
          stop("names_additive is required when format_additive = TRUE")
        }
        for (q in seq_len(length(names_additive))) {
          if (!names_additive[q] %in% colnames(sfsites)) {
            warning(paste0(names_additive[q], " is not found in obs"), call. = FALSE)
          } else {
            ## Extract afv column
            tmp.col <- st_drop_geometry(sfsites[, names_additive[q]])
            ## convert to character if column is numeric
            if (!inherits(tmp.col[, 1], "numeric")) {
              message(paste0(names_additive[q], " is not numeric. AFV values were not modified"))
            } else {
              tmp.col2 <- formatC(tmp.col[, 1],
                digits = 254, format = "f",
                drop0trailing = TRUE
              )
              sfsites[, names_additive[q]] <- tmp.col2
            }
          }
        }
      }
    } else {
      stop("This shapefile does not exist: ", obs_sites)
    }
  }

  if (!missingpreds) {
    # Append .shp if necessary
    nofe <- !grepl(".shp$", predpts)
    if (any(nofe)) predpts[nofe] <- paste0(predpts[nofe], ".shp")

    ## Check for relative pathnames
    ## rel.name <- substr(predpts, start = 1, stop = 2)== "./"
    ## if(any(rel.name)) predpts[rel.name]<- basename(predpts[rel.name])

    # Check the files exist and then import prediction points
    if (all(file.exists(predpts))) {
      sfpreds <- vector(mode = "list", length = length(predpts))

      for (m in seq_len(length(predpts))) {
        tmp.preds <- st_read(paste0(file, "/", predpts[m]), quiet = TRUE)

        ## Check geometry type
        if (sum(st_geometry_type(tmp.preds, by_geometry = TRUE) == "POINT") != nrow(tmp.preds)) {
          stop(paste0(predpts[m], " do not have POINT geometry"))
        }

        if (format_additive == TRUE) {
          if (is.null(names_additive)) {
            stop("names_additive is required when format_additive = TRUE")
          }
          for (q in seq_len(length(names_additive))) {
            if (!names_additive[q] %in% colnames(tmp.preds)) {
              warning(paste0(names_additive[q], " is not found in ", predpts[m]), call. = FALSE)
            } else {
              ## Extract afv column
              tmp.col <- st_drop_geometry(tmp.preds[, names_additive[q]])
              ## convert to character if column is numeric
              if (!inherits(tmp.col[, 1], "numeric")) {
                message(paste0(names_additive[q], " is not numeric. AFV values were not modified"))
              } else {
                tmp.col2 <- formatC(tmp.col[, 1],
                  digits = 254, format = "f",
                  drop0trailing = TRUE
                )
                tmp.preds[, names_additive[q]] <- tmp.col2
              }
            }
          }
        }

        ## Add network geometry column
        tmp.preds[, "netgeometry"] <- paste0("SNETWORK (", paste(
          tmp.preds$netID,
          tmp.preds$rid,
          tmp.preds$upDist,
          tmp.preds$ratio,
          tmp.preds$pid,
          tmp.preds$locID
        ),
        ")",
        sep = ""
        )
        sfpreds[[m]] <- tmp.preds
        names(sfpreds)[m] <- substr(basename(predpts[m]),
          start = 1,
          stop = nchar(basename(predpts[m])) - 4
        )
        rm(tmp.preds)
      }
    } else {
      stop("At least one of these shapefiles does not exist: ", predpts)
    }
  }
  ## Import edges
  if (file.exists("edges.shp")) {
    # Get the edges
    sfedges <- st_read("edges.shp", quiet = TRUE)

    ## Check geometry type
    if (sum(st_geometry_type(sfedges, by_geometry = TRUE) == "LINESTRING") != nrow(sfedges)) {
      stop("Edges do not have LINESTRING geometry")
    }

    if (format_additive == TRUE) {
      if (is.null(names_additive)) {
        stop("names_additive is required when format_additive = TRUE")
      }
      for (q in seq_len(length(names_additive))) {
        if (!names_additive[q] %in% colnames(sfedges)) {
          warning(paste0(names_additive[q], " is not found in obs"), call. = FALSE)
        } else {
          ## Extract afv column
          tmp.col <- st_drop_geometry(sfedges[, names_additive[q]])
          ## convert to character if column is numeric
          if (!inherits(tmp.col[, 1], "numeric")) {
            message(paste0(names_additive[q], " is not numeric. AFV values were not modified"))
          } else {
            tmp.col2 <- formatC(tmp.col[, 1],
              digits = 254, format = "f",
              drop0trailing = TRUE
            )
            sfedges[, names_additive[q]] <- tmp.col2
          }
        }
      }
    }


    ## Add network geometry column to edges
    sfedges[, "netgeometry"] <- paste0("ENETWORK (", paste(
      sfedges$netID,
      sfedges$rid,
      sfedges$upDist
    ),
    ")",
    sep = ""
    )
    ssnlist <- list(edges = sfedges)
  } else {
    stop(paste0("Edges is missing from ", file))
  }

  ## }

  ## )

  ## Add observation sites and prediction sites if present
  if (missingsites) {
    ssnlist$obs <- NA
  } else {
    ssnlist$obs <- sfsites
  }

  if (!missingpreds) {
    ssnlist$preds <- sfpreds
  } else {
    ssnlist$preds <- list()
  }

  # Add the path
  ## ssnlist$path <- dirname(edges)
  ssnlist$path <- ssn_folder

  # Set custom class
  class(ssnlist) <- "SSN"

  ## Create Binary ID database
  createBinaryID(ssnlist, overwrite = overwrite)

  ## Warning when observation sites are not included if(include_obs ==
  if (include_obs == FALSE) {
    warning("SSN does not include observed sites, which are needed to fit models. If this was a mistake, set include_obs == TRUE")
  }


  # Return
  return(ssnlist)
}
