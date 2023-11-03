#' Get stream distance matrices from an \code{SSN} object
#'
#' @description Extracts the stream network distance matrices for the
#'   observation or prediction data from an \code{SSN} object.
#' @param x An \code{SSN} object
#' @param name Internal name of the dataset in the object
#'   \code{x}. For observed values, this will always be "obs", the
#'   default.  To get a stream network distance matrix for a
#'   prediction data set, the name of the dataset must be given, in
#'   quotes.
#' @details
#'
#' The internal \code{name} for observed data in objects of class
#' \code{SSN} is "obs" and it is the default. If another \code{name}
#' is specified, it must represent a prediction data set in the
#' \code{SSN} object. For \code{SSN}
#' objects, these names are obtained using the call \code{names(x$preds)}.
#'
#' Note that these are not traditional symmetric distance
#' matrices. First, distances in an \code{SSN} object represent stream
#' distance, or hydrologic distance, which is the distance between two
#' locations when movement is restricted to the branching stream
#' network. Another important difference is the distance matrices for
#' \code{SSN} objects contain the \emph{downstream only} stream
#' distance between two locations, making them asymmetric. This
#' asymmetry provides a way to store two types of spatial
#' relationships based on stream distance:
#'
#' \itemize{
#' \item{Flow-connected: Water flows from an upstream site to a
#' downstream site.}
#' \item{Flow-unconnected: Two sites reside on the
#' same stream network, but do not share flow.}
#' }
#'
#' For example, if two sites are flow-connected the downstream
#' distance from the upstream site to the downstream site is > 0,
#' while the downstream distance between the downstream site and the
#' upstream site = 0. For flow-unconnected sites, the downstream
#' distance represents the distance from each site to the closest
#' downstream junction and will be > 0 in both directions. Direction
#' is preserved, with columns representing the FROM site and rows
#' representing the TO site. Row and column names correspond to the
#' unique point identifier "pid" for each site. From this matrix, it
#' is also possible to get total stream distance (downstream +
#' upstream) between any two sites on the same network (see examples
#' for additional details).
#'
#' Stream distances are only calculated within a network and so the
#' asymmetric matrices are also stored by network.  For observation
#' data, a single square matrix of distances is returned for each
#' network, with the names based on the netID value (e.g. "dist.net1",
#' "dist.net2", etc.). However, two distance matrices ("a" and "b")
#' are required to store the downstream only distance between observed
#' and prediction sites. The label "a" represents the downstream
#' stream distance \emph{from} prediction sites \emph{to} observation
#' sites, and the label "b" represents the distance \emph{from}
#' observation sites \emph{to} predictions sites.  Thus, the list of
#' prediction matrices are labeled "dist.net1.a" for the downstream
#' only distance from prediction sites in the columns, to observation
#' sites in the rows, for the first network. A prediction matrix
#' labeled "dist.net1.b" contains downstream distances \emph{from}
#' observation sites in the columns \emph{to} prediction sites in the
#' rows, for the first network. The downstream only distance matrices
#' for observations and predictions will be rectangular, unless the
#' number of observation and prediction locations are equal.  If the
#' argument \code{amongPreds = TRUE} was used in the function
#' \code{ssn_create_distmat}, then the distance between prediction sites
#' themselves is also returned, using the same labelling convention as
#' for among observation sites. That is, the matrices for each network
#' will be labeled "dist.net1", "dist.net2", etc., for the first and
#' second network, etc.
#'
#' @return A \code{\link[base]{list}} of asymmetric downstream only
#'   stream distance matrices, by network.
#' @seealso [ssn_create_distmat()]
#' @references
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving
#' average approach to spatial statistical models of stream
#' networks. The Journal of the American Statistical Association,
#' \bold{105(489)}, 22--24
#' @export
#' @examples
#' ## For this example only, copy MiddleFork04.ssn directory to R's
#' ## temporary directory
#' copy_lsn_to_temp()
#' ## Create an SSN object with prediction sites
#' mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
#'   predpts = "pred1km", overwrite = TRUE
#' )
#'
#' ## Create distance matrices for obs x obs, obs x preds, and preds x
#' ## preds
#' \dontrun{
#' ssn_create_distmat(mf04p,
#'   predpts = "pred1km", among_predpts = TRUE,
#'   overwrite = TRUE
#' )
#' }
#'
#' ## Check names of prediction datasets
#' names(mf04p$preds)
#'
#' ## Get list of stream distance matrices for observations
#' dist_obs <- ssn_get_stream_distmat(mf04p)
#' ## Display structure of list and names of the matrices
#' str(dist_obs)
#' names(dist_obs)
#' ## Look at first 5 rows and columns in asymmetric
#' ## downstream only distance matrix for netID == 1
#' dist_obs$dist.net1[1:5, 1:5]
#'
#' ## Create symmetric total stream distance matrix between
#' ## observations
#' strdist_2 <- dist_obs$dist.net2 + t(dist_obs$dist.net2)
#' strdist_2[5:10, 5:10]
#'
#' ## Get maximum downstream only distance between
#' ## observations on netID == 2
#' a.mat <- pmax(dist_obs$dist.net2, t(dist_obs$dist.net2))
#' a.mat[5:10, 5:10]
#'
#' ## Get minimum downstream only distance between observations. If
#' ## minimum distance == 0, sites are flow-connected
#' b.mat <- pmin(dist_obs$dist.net2, t(dist_obs$dist.net2))
#' b.mat[5:10, 5:10]
#'
#' ## Get distance matrices for pred1km
#' dist_pred1km <- ssn_get_stream_distmat(mf04p, name = "pred1km")
#' str(dist_pred1km)
#' names(dist_pred1km)
#' ## Look at first 5 rows and columns of downstream only distances
#' ## FROM prediction sites TO observed sites on netID == 1
#' dist_pred1km$dist.net1.a[1:5, 1:5]
#'
#' ## Look at downstream only stream distances among prediction
#' ## sites in pred1km on netID == 1. This is useful for block
#' ## prediction
#' dist_pred1km$dist.net1[1:5, 1:5]
ssn_get_stream_distmat <- function(x, name = "obs") {
  if (name == "Obs") name <- "obs"
  if (class(x)[[1]] != "SSN") {
    return("Object not of class SSN")
  }
  path <- paste0(x$path, "/distance/", name)
  flist <- list.files(path)
  distMats <- vector("list", length(flist))
  for (i in seq_len(length(flist))) {
    path1 <- paste0(path, "/", flist[i])
    file_handle <- file(path1, open = "rb")
    distmat <- unserialize(file_handle)
    close(file_handle)
    ordrow <- order(as.numeric(rownames(distmat)))
    ordcol <- order(as.numeric(colnames(distmat)))
    distmat <- distmat[ordrow, ordcol, drop = FALSE]
    distMats[[i]] <- distmat
    nameSplit <- unlist(strsplit(flist[i], "[.]"))
    tname <- nameSplit[1]
    for (j in 2:(length(nameSplit) - 1)) {
      tname <- paste(tname, nameSplit[j], sep = ".")
    }
    flist[i] <- tname
  }
  names(distMats) <- flist
  return(distMats)
}
