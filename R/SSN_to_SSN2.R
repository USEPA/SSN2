#' Convert object from \code{SpatialStreamNetwork} class to \code{SSN} class
#'
#' @description Convert an S4 \code{SpatialStreamNetwork} object
#'   created in the SSN package to an S3 \code{SSN} object used in the
#'   SSN2 package.
#' @param object A SpatialStreamNetwork object
#' @details \command{SSN_to_SSN2()} has been made available to help users
#'   migrate from the \code{SSN} package to the updated \code{SSN2}
#'   package. It is used to convert existing S4 \code{SpatialStreamNetwork}
#'   objects *stored in saved workspaces* to the S3 \code{SSN} class object
#'   used in the \code{SSN2} package. Note that \code{\link[SSN2]{ssn_import}} is
#'   used to create an S3 \code{SSN} object from data stored locally in a .ssn
#'   directory.
#'
#' @return An S3 \code{SSN} class object.
#'
#' @name SSN_to_SSN2
#' @export
#'
SSN_to_SSN2 <- function(object) {
  .Deprecated(new = "ssn_import")

  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Install the sp package before using SSN_to_SSN2", call. = FALSE)
  } else {
    if (!inherits(object, "SpatialStreamNetwork")) {
      stop("object is not of class SpatialStreamNetwork")
    }

    ## ---------------------------------------------------
    ## Convert edges to sf and add netgeom column
    ## ---------------------------------------------------

    sl <- sp::SpatialLines(object@lines, proj4string = object@proj4string)
    edges.data <- object@data
    edges <- st_as_sf(sp::SpatialLinesDataFrame(sl, edges.data, match.ID = FALSE))

    nl.coords <- object@network.line.coords

    edges <- create_netgeom(edges, type = "LINESTRING", overwrite = TRUE) 

    ## edges$netgeom <- paste("ENETWORK",
    ##   paste("(",
    ##     paste(
    ##       nl.coords$NetworkID,
    ##       nl.coords$SegmentID,
    ##       nl.coords$DistanceUpstream,
    ##       sep = " "
    ##     ),
    ##     ")",
    ##     sep = ""
    ##   ),
    ##   sep = " "
    ## )

    ## ## Convert additive function values to text
    ## if (!is.null(edge_additive)) {
    ##   if (sum(edge_additive %in% colnames(edges)) != length(edge_additive)) {
    ##     warning(paste0(
    ##       "AFV column '",
    ##       edge_additive[!edge_additive %in% colnames(edges)],
    ##       "' not found in edges"
    ##     ), call. = FALSE)
    ##   } else {
    ##     for (j in seq_len(length(edge_additive))) {
    ##       tmp <- edges[, edge_additive[j]]
    ##       tmp <- st_set_geometry(tmp, NULL)
    ##       ## tmp<- 'st_geometry<-'(tmp, NULL)
    ##       ## edges[,edge_additive[j]]<- as.character(tmp[,1])
    ##       edges[, edge_additive[j]] <- formatC(tmp[, 1],
    ##         digits = 254, format = "f",
    ##         drop0trailing = TRUE
    ##       )
    ##     }
    ##   }
    ## }


    ## ------------------------------------------------
    ## Convert observed sites to sf and add netgeom
    ## ------------------------------------------------
    ## sites<- st_as_sf(sp::SpatialPointsDataFrame(object@obspoints@SSNPoints[[1]]@point.coords,
    ##       object@obspoints@SSNPoints[[1]]@point.data,
    ##       proj4string = object@proj4string))

    sites.data <- cbind(
      object@obspoints@SSNPoints[[1]]@point.data,
      object@obspoints@SSNPoints[[1]]@point.coords
    )
    sites <- st_as_sf(sites.data, coords = c("coords.x1", "coords.x2"))
    sites <- st_set_crs(sites, object@proj4string)

    np.coords <- object@obspoints@SSNPoints[[1]]@network.point.coords
    np.coords <- cbind(np.coords, sites[, c("ratio", "locID")])
    np.coords$pid <- rownames(object@obspoints@SSNPoints[[1]]@network.point.coords)

    sites <- create_netgeom(sites, type = "POINT", overwrite = TRUE)
    ## sites$netgeom <- paste(
    ##   "SNETWORK",
    ##   paste(
    ##     "(",
    ##     paste(
    ##       np.coords$NetworkID,
    ##       np.coords$SegmentID,
    ##       np.coords$DistanceUpstream,
    ##       np.coords$ratio,
    ##       np.coords$pid,
    ##       np.coords$locID,
    ##       sep = " "
    ##     ),
    ##     ")",
    ##     sep = ""
    ##   ),
    ##   sep = " "
    ## )

    rm(np.coords)

    ## ## Convert additive function values to text
    ## if (!is.null(site_additive)) {
    ##   if (sum(site_additive %in% colnames(sites)) != length(site_additive)) {
    ##     warning(paste0(
    ##       "AFV column '",
    ##       site_additive[!site_additive %in% colnames(sites)],
    ##       "' not found in observed sites"
    ##     ), call. = FALSE)
    ##   } else {
    ##     for (j in seq_len(length(site_additive))) {
    ##       ## tmp<- 'st_geometry<-'(sites[,site_additive[j]], NULL)
    ##       tmp <- st_set_geometry(sites[, site_additive[j]], NULL)
    ##       ## sites[,site_additive[j]]<- as.character(tmp[,])
    ##       sites[, site_additive[j]] <- formatC(tmp[, ],
    ##         digits = 254, format = "f",
    ##         drop0trailing = TRUE
    ##       )
    ##     }
    ##   }
    ## }


    ## ------------------------------------------------------------- If
    ## prediction sites are present, convert to list of sf data.frames
    ## and add netgeom column
    ## -------------------------------------------------------------

    if (length(object@predpoints@ID) > 0) {
      pred.list <- list()
      for (i in seq_len(length(object@predpoints@ID))) {
        pred.name <- object@predpoints@ID[i]
        ## tmp.sf <- st_as_sf(sp::SpatialPointsDataFrame(object@predpoints@SSNPoints[[i]]@point.coords,
        ##   object@predpoints@SSNPoints[[i]]@point.data,
        ##   proj4string = object@proj4string))

        tmp.data <- cbind(
          object@predpoints@SSNPoints[[i]]@point.data,
          object@predpoints@SSNPoints[[i]]@point.coords
        )
        tmp.sf <- st_as_sf(tmp.data, coords = c("coords.x1", "coords.x2"))
        tmp.sf <- st_set_crs(tmp.sf, object@proj4string)

        np.coords <- cbind(
          object@predpoints@SSNPoints[[i]]@network.point.coords,
          tmp.sf[, c("ratio", "locID")]
        )
        np.coords$pid <- rownames(object@predpoints@SSNPoints[[i]]@network.point.coords)

        tmp.sf<- create_netgeom(tmp.sf, type = "POINT", overwrite = TRUE)
        ## tmp.sf$netgeom <- paste(
        ##   "SNETWORK",
        ##   paste(
        ##     "(",
        ##     paste(
        ##       np.coords$NetworkID,
        ##       np.coords$SegmentID,
        ##       np.coords$DistanceUpstream,
        ##       np.coords$ratio,
        ##       np.coords$pid,
        ##       np.coords$locID,
        ##       sep = " "
        ##     ),
        ##     ")",
        ##     sep = ""
        ##   ),
        ##   sep = " "
        ## )


        ## ## Convert additive function values to text
        ## if (!is.null(site_additive)) {
        ##   if (sum(site_additive %in% colnames(tmp.sf)) != length(site_additive)) {
        ##     warning(paste0(
        ##       "AFV column '",
        ##       site_additive[!site_additive %in% colnames(tmp.sf)],
        ##       "' not found in prediction dataset ", pred.name
        ##     ), call. = FALSE)
        ##   } else {
        ##     for (j in seq_len(length(site_additive))) {
        ##       ## tmp<- 'st_geometry<-'(tmp.sf[,site_additive[j]], NULL)
        ##       tmp <- st_set_geometry(tmp.sf[, site_additive[j]], NULL)
        ##       ## tmp.sf[,site_additive[j]]<- as.character(tmp[,])
        ##       tmp.sf[, site_additive[j]] <- formatC(tmp[, ],
        ##         digits = 254, format = "f",
        ##         drop0trailing = TRUE
        ##       )
        ##     }
        ##   }
        ## }

        pred.list[i] <- list(tmp.sf)
        names(pred.list)[i] <- pred.name

        rm(tmp.sf, pred.name, np.coords)
      }
    }
    ## -----------------------------------------------
    ## Construct SSN object
    ## -----------------------------------------------

    ssnlist <- list(
      edges = edges,
      obs = sites
    )

    if (exists("pred.list")) {
      ssnlist$preds <- pred.list
    } else {
      ssnlist$preds <- list()
    }

    ssnlist$path <- object@path

    class(ssnlist) <- "SSN"

    return(ssnlist)
  }
}

# # example
# ## Load saved workspace that contains S4 SpatialStreamNetwork
# # object created using SSN package
# load("old_SSN.rda")
#
# ## Load SSN2 library
# library(SSN2)
#
# ## Convert S4 SpatialStreamNetworkObject to S3 SSN object
# new_SSN <- SSN_to_SSN2(old_SSN)
# )
