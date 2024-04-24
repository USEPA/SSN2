#' Import prediction points into an SSN, ssn_lm, or ssn_glm object
#'
#' @description A shapefile of prediction points found in the .ssn
#'   directory are imported into an existing object of class
#'   \code{SSN}, \code{ssn_lm}, or \code{ssn_glm}.
#' @param x An object of class\code{SSN}, \code{ssn_lm}, or
#'   \code{ssn_glm}.
#' @param predpts Name of the prediction point dataset to import in
#'   character format. See details.
#'
#' @details \command{ssn_import_predpts} imports one set of prediction
#'   points residing in the .ssn directory into an existing
#'   \code{SSN}, \code{ssn_lm}, or \code{ssn_glm} object. The
#'   prediction dataset must be in shapefile or geopackage format
#'   (.shp or .gpkg, respectively) and reside in the ssn.object$path
#'   directory. The path for an \code{SSN} object can be updated using
#'   \command{ssn_update_path()} prior to importing prediction
#'   datasets. The argument \code{predpts} accepts the name of the
#'   prediction point dataset, with or without the file extension. If
#'   it is passed as a named vector (of length 1), then the name
#'   provided is used as the prediction dataset name in the \code{SSN}
#'   object prediction sites list
#'   (e.g. \code{names(ssn.obj$preds)}). Otherwise, the file basename
#'   is used in the names attribute. See
#'   \code{\link[SSN2]{ssn_import}} for a detailed description of the
#'   prediction dataset format within the \code{SSN} class object.
#'
#'   The prediction dataset specified in \code{predpts} must contain the
#'   spatial, topological and attribute information needed to make
#'   predictions using an ssn_lm or ssn_glm object. This information
#'   is generated using the \code{SSNbler} package, which makes use of
#'   the functionality found in the \code{sf} and \code{igraph}
#'   packages to process streams data in vector format.
#' @return an object of class \code{SSN}, \code{ssn_lm}, or
#'   \code{ssn_glm} which contains the new prediction dataset.
#'
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
#' ## Import pred1km prediction dataset into SSN object and assign the
#' ## name preds1
#' mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"))
#' mf04p <- ssn_import_predpts(mf04p, predpts = c(preds1 = "pred1km"))
#' names(mf04p$preds)
#'
#' ## Import CapeHorn prediction dataset into a ssn_glm object, using
#' ## the default file basename as the name
#' ssn_gmod <- ssn_glm(Summer_mn ~ netID, mf04p,
#'   family = "Gamma",
#'   tailup_type = "exponential", additive = "afvArea"
#' )
#' ssn_gmod <- ssn_import_predpts(ssn_gmod, predpts = "CapeHorn")
#' names(ssn_gmod$ssn.object$preds)
#' 
ssn_import_predpts <- function(x, predpts) {
  
  obj.type <- class(x)

  old_wd <- getwd()
  on.exit(setwd(old_wd))

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
      
      shp.ext<- substr(predpts[d], nchar(predpts[d])-3,
                       nchar(predpts[d])) == ".shp"
      gpkg.ext <- substr(predpts[d], nchar(predpts[d])-4,
                         nchar(predpts[d])) == ".gpkg"
      
      p.names[d] <- ifelse(
        shp.ext == TRUE, substr(predpts[d], 1, nchar(predpts[d])-4),

                 ifelse(gpkg.ext == TRUE,
                        substr(predpts[d], 1, nchar(predpts[d])-5),
                        
                 ifelse(shp.ext == FALSE & gpkg.ext == FALSE,
                        predpts[d])))
    }
  }

  ## Remove path to predpts files, if included
  predpts<- basename(predpts)

  ################################################


  ## For fitted model objects- check if predpts already exists
  if (obj.type %in% c("ssn_lm", "ssn_glm")) {
    setwd(x$ssn.object$path)
    count <- 0

    if (length(x$ssn.object$preds) > 0) {
      for (m in seq_len(length(x$ssn.object$preds))) {
        if (names(x$ssn.object$preds[m]) == p.names) {
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
        if (names(x$preds[m]) == p.names) {
          pred.num <- m
          count <- count + 1
        }
      }
    }

    if (count > 0) {
      stop("SSN already contains a prediction dataset named ", predpts)
    }
  }

  ## Import predpts as sf object
  predpoints <- get_sf_obj(predpts)
 
  ## Check geometry type
  if (sum(st_geometry_type(predpoints, by_geometry = TRUE) == "POINT") != nrow(predpoints)) {
    stop(paste0(predpts, " does not have POINT geometry"))
  }

  ## Add netgeom column
  predpoints[, "netgeom"] <- paste0("SNETWORK (", paste(
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
    names(x$preds)[[pred.num + 1]] <- p.names
  }
  ## Put prediction points in fitted model object
  if (obj.type %in% c("ssn_glm", "ssn_lm")) {
    pred.num <- length(x$ssn.object$preds)
    x$ssn.object$preds[[pred.num + 1]] <- predpoints
    names(x$ssn.object$preds)[[pred.num + 1]] <- p.names
  }
  return(x)
}
