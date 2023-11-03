#' Calculate Hydrologic Distances for an \code{SSN} object
#'
#' @description Creates a collection of (non-symmetric) matrices
#'   containing pairwise downstream hydrologic distances between sites
#'   in an \code{SSN} object
#'
#' @param ssn.object An \code{SSN} object
#' @param predpts name of prediction points in an \code{SSN} object.
#'   When a vector with length greater than one, each name is iterated upon.
#'   Default is NULL.
#' @param overwrite Logical. If \code{TRUE}, overwrite existing distance
#'   matrices. Defaults to \code{FALSE}.
#' @param among_predpts Logical. If \code{TRUE}, compute the pairwise distances
#'   between the prediction sites. Defaults to \code{FALSE}.
#' @param only_predpts Logical. If \code{TRUE}, only compute distances for
#'   prediction sites. Defaults to \code{FALSE}.
#'
#' @details A distance matrix that contains the hydrologic distance
#'  between any two sites in \code{SSN} object is needed to fit a spatial
#'  statistical model using the tail-up and tail-down autocovariance
#'  functions described in Ver Hoef and Peterson (2010). These models
#'  are implemented in R via \command{ssn_lm} and \command{ssn_glm} in
#'  the\code{SSN2} package. The hydrologic distance information needed to
#'  model the covariance between flow-connected (i.e. water flows
#'  from one location to the other) and flow-unconnected (i.e. water
#'  does not flow from one location to the other, but they reside on
#'  the same network) locations differs. The total hydrologic
#'  distance is a directionless measure; it represents the hydrologic
#'  distance between two sites, ignoring flow direction. The
#'  hydrologic distance from each site to a common downstream stream
#'  junction is used when creating models for flow-unconnected pairs,
#'  which we term downstream hydrologic distance. In contrast, the
#'  total hydrologic distance is used for modeling flow-connected
#'  pairs, which we term total hydrologic distance.
#'
#'  A downstream hydrologic distance matrix provides enough
#'  information to meet the data requirements for both the tail-up and
#'  tail-down models. When two locations are flow-connected, the
#'  downstream hydrologic distance from the upstream location to the
#'  downstream location is greater than zero, but it is zero in the
#'  other direction. When two locations are flow-unconnected the
#'  downstream hydrologic distance will be greater than zero in both
#'  directions. A site's downstream hydrologic distance to itself is
#'  equal to zero. The format of the downstream hydrologic distance
#'  matrix is efficient because distance information needed to fit
#'  both the tail-up and tail-down models is only stored once. As an
#'  example, a matrix containing the total hydrologic distance between
#'  sites is easily calculated by adding the downstream distance
#'  matrix to its transpose.
#'
#'  The downstream hydrologic distances are calculated based on the
#'  binaryIDs and stored as matrices. The matrices are stored in a
#'  directory named \sQuote{distance}, which is created by the
#'  \command{ssn_create_distmat} function within the .ssn directory. The distance
#'  directory will always contain at least one directory named
#'  \sQuote{obs}, which contains a number of .RData files, one for each
#'  network that has observed sites residing on it. The naming
#'  convention for the files is based on the netID number
#'  (e.g. dist.net1.RData). Each matrix in the \sQuote{obs} folder
#'  contains the information to form a square matrix, which contains
#'  the downstream hydrologic distance between each pair of observed
#'  sites on the network. Direction is preserved, with columns
#'  representing the FROM site and rows representing the TO site. Row
#'  and column names correspond to the pid attribute for each site.
#'
#'  If the argument \code{predpts} is specified in the call to the
#'  function, the downstream hydrologic distances between the observed
#'  and prediction sites will also be computed. A new directory is
#'  created within the distance directory, with the name corresponding
#'  to the names attribute for the preds
#'  (e.g. \code{attributes(ssn.object$preds)$names}). A sequence of
#'  .RData files is created within this directory, similar to the
#'  structure for the observed sites, except that two objects are
#'  stored for each network that contains \emph{both} observed and
#'  prediction sites. The letters \code{a} and \code{b} are used in
#'  the naming convention to distinguish between the two objects
#'  (e.g. dist.net1.a and dist.net1.b). The matrices that these
#'  objects represent are not necessarily square. In matrices of type
#'  \code{a}, rows correspond to observed locations and columns to
#'  prediction locations. In contrast, rows correspond to prediction
#'  locations and columns to observed locations in matrices of type
#'  \code{b}. Direction is also preserved, with columns representing
#'  the FROM site and rows representing the TO site in both object
#'  types. Again, row and column names correspond to the pid attribute
#'  for each site.
#'
#'  If \code{among_predpts = TRUE}, the downstream
#'  hydrologic distances will also be computed between prediction
#'  sites, for each network. Again these are stored within the distance
#'  directory with the name corresponding to the prediction points
#'  dataset. The naming convention for these prediction to prediction
#'  site distance matrices is the same as the distance matrices stored
#'  in the \sQuote{obs} directory (e.g. dist.net1.RData). These extra
#'  distance matrices are needed to perform block Kriging using
#'  \code{\link[SSN2]{predict.ssn_lm}}.
#'
#'  If \code{only_predpts = TRUE}, the downstream
#'  hydrologic distances will not be calculated between observed sites
#'  themselves. Pairwise distances will only be calculated for observed
#'  and prediction locations and. Pairwise distances between prediction
#'  locations will also be calculated if \code{among_predpts = TRUE}.
#'
#'
#' @return The \command{ssn_create_distmat} function creates a collection
#'   of hierarchical directories in the \code{ssn$path} directory,
#'   which store the pairwise distances between sites associated with
#'   the \code{SSN} object. See details section for additional information.
#'
#' @export
#'
#' @examples
#' ## Copy the MiddleForke04.ssn data to a local temporary directory.
#' ## Only needed for this example.
#' copy_lsn_to_temp()
#' ## Import SSN data
#' mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
#'   predpts = c("pred1km.shp", "Knapp"),
#'   overwrite = TRUE
#' )
#'
#' ## Create distance matrices for observations and one set of prediction sites
#' ## Include hydrologic distance matrices among prediction sites.
#' ssn_create_distmat(mf04p,
#'   predpts = "pred1km", overwrite = TRUE,
#'   among_predpts = TRUE
#' )
#'
#' ## Create distance matrices for an additional set of prediction points.
#' ## Distance matrices for observations and pred1km prediction sites are
#' ## not recalculated.
#' ssn_create_distmat(mf04p,
#'   predpts = "Knapp", overwrite = TRUE,
#'   among_predpts = TRUE, only_predpts = TRUE
#' )
ssn_create_distmat <- function(ssn.object, predpts = NULL, overwrite = FALSE,
                               among_predpts = FALSE, only_predpts = FALSE) {

  # iterate if predpts is not null
  if (length(predpts) > 1) {
    if (only_predpts) {
      x <- lapply(predpts, function(x) ssn_create_distmat(ssn.object,
                                                          predpts = x,
                                                          overwrite,
                                                          among_predpts,
                                                          only_predpts))
    } else {
      x1 <- ssn_create_distmat(ssn.object, overwrite = overwrite)
      x2 <- lapply(predpts, function(x) ssn_create_distmat(ssn.object,
                                                           predpts = x,
                                                           overwrite,
                                                           among_predpts,
                                                           only_predpts = TRUE))
    }
    return(invisible(NULL))
  }


  # recreate function argument name
  ssn <- ssn.object

  ## ------------------------------------------------------------------------
  ## Check arguments, format data
  ## ------------------------------------------------------------------------

  ## Check that arguments are valid
  if (among_predpts && (missing(predpts) || is.null(predpts))) {
    stop("A named collection of prediction points must be specified via the predpts option when among_predpts is TRUE")
  }
  ## Check to see whether distance folder exists...
  if (!file.exists(file.path(ssn$path, "distance"))) {
    dir.create(file.path(ssn$path, "distance"))
  }

  ## Check whether obs folder exists and if not, create it
  if (only_predpts == FALSE) {
    if (!file.exists(file.path(ssn$path, "distance", "obs"))) {
      dir.create(file.path(ssn$path, "distance", "obs"))
    }
  } else {
    ## Send warning if observed distance matrices are missing
    if (!file.exists(file.path(ssn$path, "distance", "obs"))) {
      warning("only_predpts == TRUE and distance matrices for observed sites are missing. Set only_predpts = FALSE to calculate distance matrices for observed and prediction sites.", call. = FALSE)
    }
  }

  ## If predpts exists
  if (!is.null(predpts)) {
    ## Check whether predpts exists in the SSN
    if (!predpts %in% attributes(ssn$preds)$names) {
      stop(predpts, " does not exist in SSN")
    }
    ## Check whether there are two sets of preds with the name predpts
    if (sum(predpts %in% attributes(ssn$preds)$names) > 1) {
      stop("SSN contains more than one copy of ", predpts)
    }
    ## Create predpts folder in distance matrix if it does not exist
    if (!file.exists(file.path(ssn$path, "distance", predpts))) {
      dir.create(file.path(ssn$path, "distance", predpts))
    }
    ## Extract netgeometry data from predpts and format
    tmp.df <- ssn_get_data(ssn, predpts)
    n.geom <- ssn_get_netgeometry(ssn$preds[[predpts]])
    colnames(n.geom)[4:6] <- paste0("ng.", colnames(n.geom[4:6]))
    n.geom$NetworkID <- as.factor(n.geom$NetworkID)
    n.geom$DistanceUpstream <- as.numeric(n.geom$DistanceUpstream)
    tmp.df <- cbind(tmp.df, n.geom)
    ssn <- ssn_put_data(tmp.df, ssn, predpts)
    rm(tmp.df, n.geom)
  }

  ## Extract netgeometry and format obs data
  tmp.df <- ssn_get_data(ssn)
  n.geom <- ssn_get_netgeometry(ssn$obs)
  colnames(n.geom)[4:6] <- paste0("ng.", colnames(n.geom[4:6]))
  n.geom$NetworkID <- as.factor(n.geom$NetworkID)
  n.geom$DistanceUpstream <- as.numeric(n.geom$DistanceUpstream)
  tmp.df <- cbind(tmp.df, n.geom)
  ssn <- ssn_put_data(tmp.df, ssn)
  rm(tmp.df)

  ## Get netID with observed or predicted sites
  if (!is.null(predpts)) {
    site.nets<- unique(c(levels(ssn$obs$NetworkID),
                         levels(ssn$preds[[predpts]]$NetworkID)))
  } else {
    site.nets<- unique(c(levels(ssn$obs$NetworkID)))
  }

  net.count <- length(site.nets)
  warned.overwrite <- FALSE

  ## Extract netgeometry and format edges data
  ssn$edges <- cbind(ssn$edges, ssn_get_netgeometry(ssn$edges))
  ssn$edges$NetworkID <- as.factor(ssn$edges$NetworkID)
  ssn$edges$DistanceUpstream <- as.numeric(ssn$edges$DistanceUpstream)

  ## ------------------------------------------------------------------
  ## Initialise binaryID.db
  ## ------------------------------------------------------------------

  if (file.exists(file.path(ssn$path, "binaryID.db")) == FALSE) {
    stop("binaryID.db is missing from SSN object. Use ssn_import() to create it.")
  }

  ## Connect to SQLite database
  connect.name <- file.path(ssn$path, "binaryID.db")

  connect <- dbConnect(SQLite(), connect.name)

  ## close sqlite connection upon function exit
  on.exit({
    dbDisconnect(connect)
  })


  ## ----------------------------------------------------------------
  ## FOR EACH NETWORK
  ## ----------------------------------------------------------------
  for (i in seq_len(net.count)) {
    ## Set network number and name
    ##net.num <- levels(ssn$edges$NetworkID)[i]
    net.num <- site.nets[i]
    net.name <- paste("net", net.num, sep = "")

    ## Get indicator for sites on this network
    ind.obs <- ssn$obs$NetworkID == net.num

    ## figure out how many observed and prediction sites there are in the network
    site.no <- nrow(ssn$obs[ind.obs, ])

    if (!is.null(predpts)) {
      ind.preds <- ssn$preds[[predpts]]$NetworkID == net.num
      pred.site.no <- nrow(ssn$preds[[predpts]][ind.preds, ])
    } else {
      pred.site.no <- 0
    }

    ## -------------------------------------------------------------
    ## CASE 1: IF OBS and PREDS EXIST ON NETWORK
    ## -------------------------------------------------------------
    if (site.no > 0 & pred.site.no > 0) {
      ## get sorted pids to use as dim names
      obs.pids <- sort(as.numeric(ssn$obs$ng.pid[ind.obs]))

      ## Get pred pids and sort
      pred.pids <- sort(as.numeric(ssn$preds[[predpts]]$ng.pid[ind.preds]))

      ## create o x p distance matrix full of NA
      current_distance_matrix_a <- matrix(NA,
        nrow = site.no, ncol = pred.site.no,
        dimnames = list(obs.pids, pred.pids)
      )

      ## create p x o distance matrix full of NA
      current_distance_matrix_b <- matrix(NA,
        nrow = pred.site.no, ncol = site.no,
        dimnames = list(pred.pids, obs.pids)
      )

      ## Set distance matrix names for o x p
      workspace.name.a <- paste("dist.net", net.num, ".a.RData", sep = "")
      workspace.name.b <- paste("dist.net", net.num, ".b.RData", sep = "")

      ## Extract binaryID table for the network
      bin.table <- dbReadTable(connect, net.name)

      ## Set o x o distance matrix name

      ##
      if (!only_predpts) {
        workspace.name <- paste("dist.net", net.num, ".RData", sep = "")
      }


      ## If overwrite == FALSE, check whether distance matrices exist
      if (!overwrite) {
        ## Check for o x o
        if (!only_predpts) {
          exists <- file.exists(file.path(ssn$path, "distance", "obs", workspace.name))
        } else {
          exists <- FALSE
        }

        ## Check for o x p
        exists <- c(
          exists, file.exists(file.path(
            ssn$path,
            "distance", predpts, workspace.name.a
          )),
          file.exists(file.path(
            ssn$path, "distance",
            predpts, workspace.name.b
          ))
        )

        ## If they are all already there, warn but continue
        if (all(exists)) {
          if (!warned.overwrite) {
            warned.overwrite <- TRUE
            message("Distance matrices already exist with overwrite set to FALSE. Not overwriting existing matrices\n")
          }
          next

          ## if some exist but some don't that's an error
        } else if (any(exists) && any(!exists)) {
          stop("overwrite was set to FALSE and some (but not all) distance matrices already exist", call. = FALSE)
        }
      }

      ## Create empty o x o distance matrix
      if (!only_predpts) {
        current_distance_matrix <- matrix(NA,
          nrow = site.no,
          ncol = site.no,
          dimnames = list(obs.pids, obs.pids)
        )
        diag(current_distance_matrix) <- 0

        rownames(current_distance_matrix) <- obs.pids
        colnames(current_distance_matrix) <- obs.pids
      }

      ## Get locID for observed sites on all networks
      ## OLD CODE RETURNS attributes(network.point.coords)$locID
      ## In practice this returns all locIDs in numeric format
      ## This is returning a subset in character format
      ## locID.obi<- ssn$obs$ng.locID[ind.obs]
      ## locID.obi<- ssn$obs$ng.locID

      ## Create data.frame for obs with columns pid, rid, locID
      ob.i <- ssn_get_netgeometry(ssn$obs[ind.obs, ], c("pid", "SegmentID", "locID"),
                                  reformat = TRUE)
      ##ob.i <- as.data.frame(sapply(ob.i, as.numeric))
      colnames(ob.i) <- c("pid", "rid", "locID")
      ob.i$locID <- as.factor(ob.i$locID)

      ## Get binaryID for rid observed sites reside on and order by pid
      ob.i$binaryID <- bin.table$binaryID[match(ob.i$rid, bin.table$rid)]
      ob.i <- ob.i[order(ob.i[, "pid"]), ]
      rownames(ob.i) <- ob.i$pid

      ## Create a new copy ordered by locID factor level and save vector for reordering
      ob.i_by_locID <- ob.i[order(ob.i[, "locID"]), ]
      ## ob.i_by_locID$pid <- as.numeric(ob.i_by_locID$pid)

      ## Convert locID factor level (not locID) back to numeric and
      ## reorder by pid
      ob.i_by_locID$locID <- as.numeric(ob.i_by_locID$locID)
      ob.j_reordering <- order(ob.i_by_locID$pid)

      ## Initialise locID.old
      locID.old <- -1

      ## Create a vector of locIDs that are not duplicated
      ind.dup <- !duplicated(ob.i_by_locID$locID)

      ## Calculate distances between observed sites
      for (j in seq_len(NROW(ob.i))) {
        pid.i <- ob.i[j, "pid"]
        locID.i <- ob.i[j, "locID"]
        ind <- ssn$obs$ng.pid == pid.i
        upDist.i <- ssn$obs$DistanceUpstream[ind]

        if (!only_predpts) {
          if (locID.i != locID.old) {
            junk <- get.rid.fc(ob.i_by_locID[ind.dup, "binaryID"], ob.i$binaryID[j])
            ob.j <- getObsRelationshipsDF(
              ssn, pid.i, junk, ind.dup, ob.i,
              ob.i_by_locID, bin.table
            )

            ob.j <- ob.j[ob.j_reordering, ]

            ## obs fills in by column because it's working between obs to obs.
            ind.fc <- ob.j$fc == 1

            dist.obs <- ifelse(ind.fc, upDist.i - ob.j$upDist.j,
              upDist.i - ob.j$juncDist
            )
            current_distance_matrix[, paste(pid.i)] <- ifelse(dist.obs < 0, 0, dist.obs)
          } else {
            current_distance_matrix[, paste(pid.i)] <-
              current_distance_matrix[, paste(pid.old)]
          }
        }

        ## Calculate distance from observed site to pred sites
        if (locID.i != locID.old) {
          ob.j <- getPredRelationshipsDF(ssn, predpts, ind.preds, bin.table, ob.i, j)
          ob.j <- ob.j[order(ob.j[, "pid"]), ]

          ind.fc <- ob.j$fc == 1
          dist.a <- ifelse(ind.fc, ob.j$upDist.j - upDist.i, ob.j$upDist.j - ob.j$juncDist)
          current_distance_matrix_a[paste(pid.i), ] <- ifelse(dist.a < 0, 0, dist.a)

          dist.b <- ifelse(ind.fc, upDist.i - ob.j$upDist.j, upDist.i - ob.j$juncDist)
          current_distance_matrix_b[, paste(pid.i)] <- ifelse(dist.b < 0, 0, dist.b)
        } else {
          ## add column to pred sites
          ## if (!is.null(predpts) && pred.site.no > 0) {
          current_distance_matrix_a[paste(pid.i), ] <-
            current_distance_matrix_a[paste(pid.old), ]
          current_distance_matrix_b[, paste(pid.i)] <-
            current_distance_matrix_b[, paste(pid.old)]
          ## }
        }

        pid.old <- pid.i
        locID.old <- locID.i
      }

      ## Write distance matrices to distance folder --------------------------------------
      ## Observed distance matrix
      if (!only_predpts) {
        file_handle <- file(file.path(ssn$path, "distance", "obs", workspace.name), open = "wb")
        serialize(current_distance_matrix, file_handle, ascii = FALSE)
        close(file_handle)
      }

      ## save obs-pred and pred-obs distance matrices
      file_handle <- file(file.path(ssn$path, "distance", predpts, workspace.name.a), open = "wb")
      serialize(current_distance_matrix_a, file_handle, ascii = FALSE)
      close(file_handle)

      file_handle <- file(file.path(ssn$path, "distance", predpts, workspace.name.b), open = "wb")
      serialize(current_distance_matrix_b, file_handle, ascii = FALSE)
      close(file_handle)
    }

    ## ---------------------------------------------------------------------------------
    ## Calculate distances among observed sites
    ## ---------------------------------------------------------------------------------
    if (site.no > 0 & pred.site.no == 0 & only_predpts == FALSE) {
      ## Create distance matrix full of NAs. If the file already exists
      ## and overwrite == FALSE, stop with an error
      workspace.name <- paste("dist.net", net.num, ".RData", sep = "")
      obs.pids <- sort(as.numeric(ssn$obs$ng.pid[ind.obs]))

      net.name <- paste("net", net.num, sep = "")
      bin.table <- dbReadTable(connect, net.name)

      obs_distance_matrix <- amongSitesDistMat(ssn, obs.pids, name = "obs", bin.table)
      ## Write obs-obs distance matrix to distance folder
      file_handle <- file(file.path(ssn$path, "distance", "obs", workspace.name), open = "wb")
      serialize(obs_distance_matrix, file_handle, ascii = FALSE)
      close(file_handle)
    }


    ## ---------------------------------------------------------------------------------
    ## Calculate distances among prediction sites
    ## ---------------------------------------------------------------------------------
    if (among_predpts & pred.site.no > 0) {
      ## Create distance matrix full of NAs. If the file already exists
      ## and overwrite == FALSE, stop with an error
      workspace.name <- paste("dist.net", net.num, ".RData", sep = "")
      pred.pids <- sort(as.numeric(ssn$preds[[predpts]]$ng.pid[ind.preds]))

      net.name <- paste("net", net.num, sep = "")

      ## This uses bin.table from last network if site.no == 0
      if (site.no == 0 | !exists("bin.table")) {
        bin.table <- dbReadTable(connect, net.name)
      }

      among_distance_matrix <- amongSitesDistMat(ssn, pred.pids, name = predpts, bin.table)
      ## Write preds-preds distance matrix to distance folder
      file_handle <- file(file.path(ssn$path, "distance", predpts, workspace.name),
        open = "wb"
      )
      serialize(among_distance_matrix, file_handle, ascii = FALSE)
      close(file_handle)
    }
  }
}
