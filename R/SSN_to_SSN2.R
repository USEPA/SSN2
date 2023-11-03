#' Convert object from \code{SpatialStreamNetwork} class to \code{SSN} class
#'
#' @description Convert an S4 \code{SpatialStreamNetwork} object
#'   created in the SSN package to an S3 \code{SSN} object used in the
#'   SSN2 package.
#' @param object A SpatialStreamNetwork object
#' @param edge_additive A character vector of additive function value
#'   column names found in edges. Default is NULL. See Details for
#'   more information.
#' @param site_additive A character vector of additive function value
#'   column names found in the observed sites and prediction
#'   sites. See Details for more information. Default is NULL.
#' @details \command{SSN_to_SSN2()} has been made available to help users
#'   migrate from the \code{SSN} package to the updated \code{SSN2}
#'   package. It is used to convert existing S4 \code{SpatialStreamNetwork}
#'   objects *stored in saved workspaces* to the S3 \code{SSN} class object
#'   used in the \code{SSN2} package. Note that \code{\link[SSN2]{ssn_import}} is
#'   used to create an S3 \code{SSN} object from data stored locally in a .ssn
#'   directory.
#'
#'   Additive function values are used to generate spatial weights for
#'   the tail-up covariance function used in \code{ssn_glm}. The range
#'   of additive function values are restricted to \eqn{0 \le AFV \le
#'   1}. In the \code{SSN2} package, columns containing additive
#'   function values are stored as text, rather than numeric
#'   format. This prevents values less than 1 with more than 10 digits
#'   from being truncated when writing/reading shapefiles (and their
#'   .dbf tables). The columns containing additive function values are
#'   specified using the \code{edge_additive} and \code{site_additive} arguments
#'   and converted to character format in the \code{SSN} class object
#'   returned. The arguments \code{edge_additive} and \code{site_additive}
#'   accept a single column name in character format, or a vector
#'   containing multiple column names. Note that, column names for
#'   additive function values on the edges, sites, and prediction
#'   sites may differ. If a column specified in \code{edge_additive} or
#'   \code{site_additive} is not present, the function will return a
#'   warning, rather than an error. Columns containing additive
#'   function values can also be converted to text manually using the
#'   \code{\link[base]{formatC}} function, which provides the
#'   flexibility needed to store the values with their full precision.
#'
#' @return An S3 \code{SSN} class object, with additive function value
#'   columns converted to text format.
#'
#' @name SSN_to_SSN2
#' @export
SSN_to_SSN2 <- function(object, edge_additive = NULL, site_additive = NULL) {
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Install the sp package before using SSN_to_SSN2", call. = FALSE)
  } else {
    if (!inherits(object, "SpatialStreamNetwork")) {
      stop("object is not of class SpatialStreamNetwork")
    }

    ## ---------------------------------------------------
    ## Convert edges to sf and add netgeometry column
    ## ---------------------------------------------------

    sl <- sp::SpatialLines(object@lines, proj4string = object@proj4string)
    edges.data <- object@data
    edges <- st_as_sf(sp::SpatialLinesDataFrame(sl, edges.data, match.ID = FALSE))

    nl.coords <- object@network.line.coords

    edges$netgeometry <- paste("ENETWORK",
      paste("(",
        paste(
          nl.coords$NetworkID,
          nl.coords$SegmentID,
          nl.coords$DistanceUpstream,
          sep = " "
        ),
        ")",
        sep = ""
      ),
      sep = " "
    )

    ## Convert additive function values to text
    if (!is.null(edge_additive)) {
      if (sum(edge_additive %in% colnames(edges)) != length(edge_additive)) {
        warning(paste0(
          "AFV column '",
          edge_additive[!edge_additive %in% colnames(edges)],
          "' not found in edges"
        ), call. = FALSE)
      } else {
        for (j in seq_len(length(edge_additive))) {
          tmp <- edges[, edge_additive[j]]
          tmp <- st_set_geometry(tmp, NULL)
          ## tmp<- 'st_geometry<-'(tmp, NULL)
          ## edges[,edge_additive[j]]<- as.character(tmp[,1])
          edges[, edge_additive[j]] <- formatC(tmp[, 1],
            digits = 254, format = "f",
            drop0trailing = TRUE
          )
        }
      }
    }


    ## ------------------------------------------------
    ## Convert observed sites to sf and add netgeometry
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
    sites$netgeometry <- paste(
      "SNETWORK",
      paste(
        "(",
        paste(
          np.coords$NetworkID,
          np.coords$SegmentID,
          np.coords$DistanceUpstream,
          np.coords$ratio,
          np.coords$pid,
          np.coords$locID,
          sep = " "
        ),
        ")",
        sep = ""
      ),
      sep = " "
    )

    rm(np.coords)

    ## Convert additive function values to text
    if (!is.null(site_additive)) {
      if (sum(site_additive %in% colnames(sites)) != length(site_additive)) {
        warning(paste0(
          "AFV column '",
          site_additive[!site_additive %in% colnames(sites)],
          "' not found in observed sites"
        ), call. = FALSE)
      } else {
        for (j in seq_len(length(site_additive))) {
          ## tmp<- 'st_geometry<-'(sites[,site_additive[j]], NULL)
          tmp <- st_set_geometry(sites[, site_additive[j]], NULL)
          ## sites[,site_additive[j]]<- as.character(tmp[,])
          sites[, site_additive[j]] <- formatC(tmp[, ],
            digits = 254, format = "f",
            drop0trailing = TRUE
          )
        }
      }
    }


    ## ------------------------------------------------------------- If
    ## prediction sites are present, convert to list of sf data.frames
    ## and add netgeometry column
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
        tmp.sf$netgeometry <- paste(
          "SNETWORK",
          paste(
            "(",
            paste(
              np.coords$NetworkID,
              np.coords$SegmentID,
              np.coords$DistanceUpstream,
              np.coords$ratio,
              np.coords$pid,
              np.coords$locID,
              sep = " "
            ),
            ")",
            sep = ""
          ),
          sep = " "
        )


        ## Convert additive function values to text
        if (!is.null(site_additive)) {
          if (sum(site_additive %in% colnames(tmp.sf)) != length(site_additive)) {
            warning(paste0(
              "AFV column '",
              site_additive[!site_additive %in% colnames(tmp.sf)],
              "' not found in prediction dataset ", pred.name
            ), call. = FALSE)
          } else {
            for (j in seq_len(length(site_additive))) {
              ## tmp<- 'st_geometry<-'(tmp.sf[,site_additive[j]], NULL)
              tmp <- st_set_geometry(tmp.sf[, site_additive[j]], NULL)
              ## tmp.sf[,site_additive[j]]<- as.character(tmp[,])
              tmp.sf[, site_additive[j]] <- formatC(tmp[, ],
                digits = 254, format = "f",
                drop0trailing = TRUE
              )
            }
          }
        }

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
      ssnlist$preds <- NA
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
# new_SSN <- SSN_to_SSN2(old_SSN,
#   edge_additive = "afvArea",
#   site_additive = "afvArea"
# )
