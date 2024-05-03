#' Import \code{SSN} object
#'
#' @description This function reads spatial data from a .ssn folder
#'   and creates an \code{SSN} object.
#' @param path Filepath to the .ssn directory. See details.
#' @param include_obs default = \code{TRUE}. Logical indicating
#'   whether observed sites should be included in the SSN object.
#' @param predpts Vector of prediction site dataset names
#'   found within the .ssn folder. See details.
#' @param overwrite default = \code{FALSE}. If \code{TRUE}, overwrite
#'   existing binaryID.db files.
#'
#' @details The \command{ssn_import} function imports spatial data (shapefile or geopackage format)
#'   from a .ssn folder generated using the
#'   \code{SSNbler} package function \code{SSNbler::lsn_to_ssn}. The .ssn folder contains all of the spatial, topological and
#'   attribute data needed to fit a spatial statistical stream network
#'   model to streams data.  This includes:
#' \itemize{
#'     \item{An edges dataset with LINESTRING geometry representing the stream network.}
#'     \item{A sites dataset with POINT geometry where observed data were collected on
#'   the stream network.}
#'     \item{Prediction sites dataset(s) representing
#'   locations where predictions will be made.}
#'     \item{netID.dat file(s) for each distinct network, which stores the topological
#'   relationships of the line features in edges.}}
#'
#'   A more detailed description of the .ssn directory and its
#'   contents is provided in Peterson and Ver Hoef (2014).
#'
#'   The \command{ssn_import} imports the edges, observed sites (optional), and
#'   prediction sites (optional) as \code{sf data.frame} objects. A new column named 'netgeom'
#'   is created to store important data representing
#'   topological relationships in a spatial stream network
#'   model. These data are stored in character format, which is less
#'   likely to be inadvertantly changed by users. See
#'   \code{\link[SSN2]{ssn_get_netgeom}} for a more detailed description of
#'   the format and contents of 'netgeom'.
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
#'   At a minimum, an \code{SSN} object must always contain streams,
#'   which are referred to as edges. The \code{SSN} object would also
#'   typically contain a set of observed sites, where measurements
#'   have been collected. Only one observed dataset is permitted in an
#'   \code{SSN} object. When \code{include_obs=FALSE}, an \code{SSN}
#'   object is created without observations. This option provides
#'   flexibility for users who would like to simulate data on a set of
#'   artifical sites on an existing stream network. Note that
#'   observation sites must be included in the \code{SSN} object in
#'   order to fit models using \command{ssn_lm} or
#'   \command{ssn_glm}. The \code{SSN} object may contain multiple
#'   sets of prediction points (or none), which are stored as separate
#'   datasets in the .ssn directory. If \code{predpts} is a named
#'   vector, the names of the \code{preds} list in the \code{SSN}
#'   object correspond to the vector names. Otherwise, they are set to
#'   the basename of the prediction site dataset file(s) specified in
#'   \code{predpts}. The \code{\link[SSN2]{ssn_import_predpts}}
#'   function allows users to import additional sets of prediction
#'   sites to a an existing \code{SSN} object.
#'
#' @return \code{ssn_import} returns an object of class SSN, which is a list
#'   with four elements containing:
#'   \itemize{
#'     \item{\code{edges}: An \code{sf data.frame} containing the stream network,
#'     with an additional 'netgeom' column.}
#'     \item{\code{obs}: An sf data.frame containing observed site locations,
#'     with an additional 'netgeom' column. NA if \code{include_obs =
#'     FALSE}.}
#'     \item{\code{preds}: A list of sf data.frames containing prediction
#'     site locations. An empty list is returned if \code{predpts} is not provided.}
#'     \item{path: The local filepath for the .ssn directory associated with the \code{SSN}
#'     object.}
#'   }
#'
#' @references
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
#'     "pred1km",
#'     "CapeHorn"
#'   ),
#'   overwrite = TRUE
#' )
#'

ssn_import <- function(path, include_obs = TRUE, predpts = NULL,
                        overwrite = FALSE) {

  if(!dir.exists(path)) stop("Cannot find the .ssn folder.")

  # Get wd
  old_wd <- getwd()
  on.exit(setwd(old_wd)) # if the function crashes or finishes, it will restore the initial working directory

  local_dir(path)

  ## Adjust if relative pathname is supplied in path
  if(substr(path, start = 1, stop = 2) == "./") {
    rel.path <- substr(path, start = 2, stop = nchar(path))
    path <- paste0(old_wd, rel.path)
  }

#################################################################
  ## Check format of predpts
################################################################

  ## If names are provided, use them
  if(is.vector(predpts) & !is.null(names(predpts))) {
    p.names<- names(predpts)
  }

  ## If no names provided assign based on name of file without extension
  if(is.vector(predpts) & is.null(names(predpts))) {
    p.names <- NULL

    for(d in 1:length(predpts)) {

      shp.ext<- substr(predpts[d], nchar(predpts[d])-3, nchar(predpts[d])) == ".shp"
      gpkg.ext <- substr(predpts[d], nchar(predpts[d])-4, nchar(predpts[d])) == ".gpkg"

      p.names[d] <- ifelse(
        shp.ext == TRUE, substr(predpts[d], 1, nchar(predpts[d])-4),

                 ifelse(gpkg.ext == TRUE,
                        substr(predpts[d], 1, nchar(predpts[d])-5),

                 ifelse(shp.ext == FALSE & gpkg.ext == FALSE, predpts[d])))
    }
  }


  ## Remove path to predpts files, if included
  if(!is.null(predpts)) {
    predpts<- basename(predpts)
  }


  ## ----------------------------------------------------
  ## Import edges
  ## ----------------------------------------------------
  sfedges <- get_sf_obj("edges")

  if(exists("sfedges")) {
    ## Check geometry type
    if(sum(st_geometry_type(sfedges, by_geometry = TRUE) ==
            "LINESTRING") != nrow(sfedges)) {
      stop("Edges do not have LINESTRING geometry")
    }

    ## Ensure geometry column is named geometry
    if(!"geometry" %in% colnames(sfedges)) {
      sf::st_geometry(sfedges) <- "geometry"
    }

    ## Add network geometry column to edges
    sfedges[, "netgeom"] <-
      paste0("ENETWORK (", paste(
                             sfedges$netID,
                             sfedges$rid,
                             sfedges$upDist),
             ")", sep = "")

  } else {
    stop(paste0("Edges is missing from ", path))
  }

  ## ----------------------------------------------------
  ## Import obs
  ## ----------------------------------------------------

  ## Check observation sites exist and import
  if(include_obs == TRUE) {

    ## Check the file exists and then ingest, process
    sfsites <- get_sf_obj("sites")

    ## Check geometry type
    if(sum(st_geometry_type(sfsites, by_geometry = TRUE) ==
           "POINT") != nrow(sfsites)) {
      stop("Observed sites do not have POINT geometry")
    }

    ## Ensure geometry column is named geometry
    if(!"geometry" %in% colnames(sfsites)) {
      sf::st_geometry(sfsites) <- "geometry"
    }

    ## ## Add network geometry column
    ##sfsites<- create_netgeom2(sfsites, type = "point")
    sfsites[, "netgeom"] <- paste0("SNETWORK (",
                                   paste(
                                     sfsites$netID, sfsites$rid, sfsites$upDist,
                                     sfsites$ratio, sfsites$pid, sfsites$locID),
                                   ")", sep = "")
  } else {
    sfsites <- NA
    # sfsites <- list()
  }

  ## ----------------------------------------------------
  ## Import preds
  ## ----------------------------------------------------
  if(!is.null(predpts)) {
    sfpreds <- vector(mode = "list", length = length(predpts))

    for (m in seq_len(length(predpts))) {
      ##tmp.preds <- st_read(paste0(file, "/", predpts[m]), quiet = TRUE)
      tmp.preds <- get_sf_obj(paste0(path, "/", predpts[m]))

       ## Check geometry type
       if (sum(st_geometry_type(tmp.preds, by_geometry = TRUE) == "POINT") != nrow(tmp.preds)) {
         stop(paste0(predpts[m], " do not have POINT geometry"))
       }

      ## Add network geometry column
      tmp.preds[, "netgeom"] <- paste0("SNETWORK (", paste(
                                                       tmp.preds$netID,
                                                       tmp.preds$rid,
                                                       tmp.preds$upDist,
                                                       tmp.preds$ratio,
                                                       tmp.preds$pid,
                                                       tmp.preds$locID
                                                     ),")",sep = "")


      sfpreds[[m]] <- tmp.preds
      ## names(sfpreds)[m] <- substr(basename(predpts[m]),
      ##                             start = 1,
      ##                             stop = nchar(basename(predpts[m])) - 4)

      names(sfpreds)[m]<-p.names[m]
      rm(tmp.preds)
      }
    ## } else {
    ##   stop("At least one of these shapefiles does not exist: ", predpts)
    ## }
  } else {
    sfpreds <- list()
  }


  #################################
  ##if(preds.exist) {

  ##   ## Check the files exist and then import prediction points
  ##   if(all(file.exists(unlist(predpts)))) {

  ##     ## Create empty list to hold sf objects and add names if missing
  ##     sfpreds <- vector(mode = "list", length = length(predpts))

  ##     if(is.null(names(predpts))) {
  ##       names(predpts) <- sub("\\.\\w+$", "", predpts)
  ##     } else {
  ##       names(sfpreds) <- names(predpts)
  ##     }

  ##     for(m in seq_len(length(predpts))) {
  ##       ## Check filename
  ##       if(grepl(".gpkg", predpts[m])) {
  ##         pred.format <- ".gpkg"
  ##       }
  ##       if(grepl(".shp", predpts[m])) {
  ##         pred.format = ".shp"
  ##       }
  ##       if(!exists("pred.format")) {
  ##         stop("All predpts must be in shapefile or geopackage format and include the file extension")
  ##       }

  ##       tmp.preds <- st_read(paste0(path, "/", predpts[[m]]), quiet = TRUE)

  ##       ## Check geometry type
  ##       if(sum(st_geometry_type(tmp.preds, by_geometry = TRUE) == "POINT") != nrow(tmp.preds)) {
  ##         stop(paste0(predpts[[m]], " do not have POINT geometry"))
  ##       }

  ##       ## Ensure geometry column is named geometry
  ##       if(!"geometry" %in% colnames(tmp.preds)) {
  ##         sf::st_geometry(tmp.preds) <- "geometry"
  ##       }

  ##       ## Add network geometry column
  ##       tmp.preds<- create_netgeom2(tmp.preds, type = "point")

  ##       sfpreds[[m]] <- tmp.preds

  ##       rm(tmp.preds)
  ##     }
  ##   } else {
  ##     stop(paste("At least one of these files does not exist: ", paste(predpts, collapse = ", ")))
  ##   }
  ## } else {
  ##   sfpreds <- list()
  ## }

  ## ---------------------------------------------------------
  ## Create SSN object and return
  ## ---------------------------------------------------------
  ssnlist <- list(edges = sfedges, obs = sfsites, preds = sfpreds, path = path)
  class(ssnlist) <- "SSN"

  ## Create Binary ID database
  createBinaryID(ssnlist, overwrite = overwrite)

  ## Warning when observation sites are not included in SSN
  if(is.logical(ssnlist$obs)) {
    warning("SSN does not include observed sites, which are needed to fit models. If this was a mistake, run ssn_import() with obs correctly defined")
  }

  ## Return
  return(ssnlist)
}

