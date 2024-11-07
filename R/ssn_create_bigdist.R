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
#' @param no_cores Number of cores to use in computation of the distance matrices. Also, the number of chunks to split the dataset into during computation.
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
#'   predpts = c("pred1km.gpkg", "CapeHorn"),
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
#'   predpts = "CapeHorn", overwrite = TRUE,
#'   among_predpts = TRUE, only_predpts = TRUE
#' )
ssn_create_bigdist <- function(ssn.object, predpts = NULL, overwrite = FALSE,
                               among_predpts = FALSE, only_predpts = FALSE,
                               no_cores = 1) {
  # iterate if predpts is not null
  ## if (length(predpts) > 1) {
  ##   if (only_predpts) {
  ##     x <- lapply(predpts, function(x) {
  ##       ssn_create_distmat(ssn.object,
  ##         predpts = x,
  ##         overwrite,
  ##         among_predpts,
  ##         only_predpts
  ##       )
  ##     })
  ##   } else {
  ##     x1 <- ssn_create_distmat(ssn.object, overwrite = overwrite)
  ##     x2 <- lapply(predpts, function(x) {
  ##       ssn_create_distmat(ssn.object,
  ##         predpts = x,
  ##         overwrite,
  ##         among_predpts,
  ##         only_predpts = TRUE
  ##       )
  ##     })
  ##   }
  ##   return(invisible(NULL))
  ## }


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
    ## Extract netgeom data from predpts and format
    tmp.df <- ssn_get_data(ssn, predpts)
    n.geom <- ssn_get_netgeom(ssn$preds[[predpts]])
    colnames(n.geom)[4:6] <- paste0("ng.", colnames(n.geom[4:6]))
    n.geom$NetworkID <- as.factor(n.geom$NetworkID)
    n.geom$DistanceUpstream <- as.numeric(n.geom$DistanceUpstream)
    tmp.df <- cbind(tmp.df, n.geom)
    ssn <- ssn_put_data(tmp.df, ssn, predpts)
    rm(tmp.df, n.geom)
  }

  ## Extract netgeom and format obs data
  tmp.df <- ssn_get_data(ssn)
  n.geom <- ssn_get_netgeom(ssn$obs)
  colnames(n.geom)[4:6] <- paste0("ng.", colnames(n.geom[4:6]))
  n.geom$NetworkID <- as.factor(n.geom$NetworkID)
  n.geom$DistanceUpstream <- as.numeric(n.geom$DistanceUpstream)
  tmp.df <- cbind(tmp.df, n.geom)
  ssn <- ssn_put_data(tmp.df, ssn)
  rm(tmp.df)

  ## Get netID with observed or predicted sites
  if (!is.null(predpts)) {
    site.nets <- unique(c(
      levels(ssn$obs$NetworkID),
      levels(ssn$preds[[predpts]]$NetworkID)
    ))
  } else {
    site.nets <- unique(c(levels(ssn$obs$NetworkID)))
  }

  net.count <- length(site.nets)
  warned.overwrite <- FALSE

  ## Extract netgeom and format edges data
  ssn$edges <- cbind(ssn$edges, ssn_get_netgeom(ssn$edges))
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

  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  ## close sqlite connection upon function exit
  on.exit({
      dbDisconnect(connect)
      stopCluster(cl)
      registerDoSEQ()
      ##closeAllConnections()
  })


  ## ----------------------------------------------------------------
  ## FOR EACH NETWORK
  ## ----------------------------------------------------------------
  for (i in seq_len(net.count)) {
    ## Set network number and name
    ## net.num <- levels(ssn$edges$NetworkID)[i]
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

      ## Set distance matrix names for o x p
      workspace.name.a <- paste(ssn$path, "/distance/", predpts,
                                "/dist.net", net.num, ".a", sep = "")
      workspace.name.b <- paste(ssn$path, "/distance/", predpts,
                                "/dist.net", net.num, ".b", sep = "")

      ## ## create o x p distance matrix full of NA
      ## current_distance_matrix_a <- matrix(NA,
      ##   nrow = site.no, ncol = pred.site.no,
      ##   dimnames = list(obs.pids, pred.pids)
        ## )
      current_distance_matrix_a <-
            fm.create(filenamebase = workspace.name.a,
              nrow = pred.site.no, ncol = site.no, type = "double")

            rownames(current_distance_matrix_a) <- as.character(pred.pids)
            colnames(current_distance_matrix_a) <- as.character(obs.pids)


      ## ## create p x o distance matrix full of NA
      ## current_distance_matrix_b <- matrix(NA,
      ##   nrow = pred.site.no, ncol = site.no,
      ##   dimnames = list(pred.pids, obs.pids)
      ## )
      current_distance_matrix_b <-
            fm.create(filenamebase = workspace.name.b,
                      nrow = pred.site.no, ncol = site.no,
                      type = "double")
      rownames(current_distance_matrix_b) <- as.character(pred.pids)
      colnames(current_distance_matrix_b) <- as.character(obs.pids)

      ##close(current_distance_matrix_b)
      ##close(current_distance_matrix_a)

      ## Extract binaryID table for the network
      bin.table <- dbReadTable(connect, net.name)

      ## If overwrite == FALSE, check whether distance matrices exist
      if (!overwrite) {

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

        # ## If they are all already there, warn but continue
        # if (all(exists)) {
        #   if (!warned.overwrite) {
        #     warned.overwrite <- TRUE
        #     message("Prediction site distance matrices already exist with overwrite set to FALSE. Not overwriting existing matrices\n")
        #   }
        #   next
        #
        # ## if some exist but some don't that's an error
        # } else if (any(exists) && any(!exists)) {
        #     stop("overwrite was set to FALSE and some (but not all) distance matrices already exist", call. = FALSE)
        # }
      }

      ##browser()
      ## Calculate obs x obs distance matrix
      if (!only_predpts) {

        ## Set o x o distance matrix name
        workspace.name1 <- paste(ssn$path, "/distance/obs/dist.net",
                                net.num, sep = "")

        # ## Check for o x o - look at this later
        # if (!only_predpts) {
        #   obs.exists <- file.exists(file.path(ssn$path, "distance", "obs",
        #                                   workspace.name))
        # } else {
        #   obs.exists <- FALSE
        # }
        #
        # if(obs.exists == TRUE){
        #   message("Observed site distance matrices already exist with overwrite set to FALSE. Not overwriting existing matrices\n")
        #
        # }

        ## Create empty obs distance matrix
        current_distance_matrix <-
            fm.create(filenamebase = workspace.name1,
                      nrow = site.no, ncol = site.no, type = "double")

        rownames(current_distance_matrix) <- as.character(obs.pids)
        colnames(current_distance_matrix) <- as.character(obs.pids)

        ##close(current_distance_matrix)



        ## Calculate obs x obs distances
        if (is.null(getDoParName())) {
            registerDoSEQ() # A little hack to avoid the foreach warning 1st time.
        }

        ## Create vector iterator
        itCol <- isplitVector(obs.pids, chunks = no_cores)

        ans1<- foreach(obs.pids.vec=itCol,
                       .packages = c("SSN2", "filematrix","itertools","iterators"),
                       .errorhandling = "pass") %dopar% {

                         amongSitesBigDistMat(ssn = ssn,
                                              pids = obs.pids.vec,
                                              net.num = net.num,
                                              bin.table = bin.table,
                                              name = "obs",
                                              workspace.name = workspace.name1)
                         finished1 <- TRUE
                         return(finished1)
                       }
        close(current_distance_matrix)

      }

      ## Create vector iterator
      if(only_predpts == TRUE) {
        itCol <- isplitVector(obs.pids, chunks = no_cores)
      }

      ##browser()
      ans2<- foreach(obs.pids.vec=itCol,
                       .packages = c("SSN2", "filematrix","itertools","iterators"),
                       .errorhandling = "stop") %dopar% {

                          amongObsPredsBigDistMat(ssn = ssn,
                                                 obs.pids = obs.pids.vec,
                                                 pred.pids=pred.pids,
                                                 bin.table = bin.table,
                                                 workspace.name.a = workspace.name.a,
                                                 workspace.name.b = workspace.name.b,
                                                 pred.name = predpts)
                           finished2 <- TRUE
                           return(finished2)
                       }
      close(current_distance_matrix_a)
      close(current_distance_matrix_b)
    }
    ## ---------------------------------------------------------------------------------
    ## Case 2: Observed sites, no prediction sites
    ## ---------------------------------------------------------------------------------
    if (site.no > 0 & pred.site.no == 0 & only_predpts == FALSE) {
      ## Create distance matrix full of NAs. If the file already exists
      ## and overwrite == FALSE, stop with an error
      workspace.name1 <- paste(ssn$path, "/distance/obs/dist.net", net.num, sep = "")

      obs.pids <- sort(as.numeric(ssn$obs$ng.pid[ind.obs]))

      net.name <- paste("net", net.num, sep = "")
      bin.table <- dbReadTable(connect, net.name)

      ## obs_distance_matrix <- amongSitesDistMat(ssn, obs.pids, name = "obs", bin.table)
      ## ## Write obs-obs distance matrix to distance folder
      ## file_handle <- file(file.path(ssn$path, "distance", "obs", workspace.name), open = "wb")
      ## serialize(obs_distance_matrix, file_handle, ascii = FALSE)
        ## close(file_handle)

      ## Create vector iterator
      itCol <- isplitVector(obs.pids, chunks = no_cores)

      ans1<- foreach(obs.pids.vec=itCol,.packages = c("SSN2", "filematrix",
                                                        "itertools","iterators"),
                       .errorhandling = "pass") %dopar% {
                           # amongObsBigDistMat(ssn = ssn, net.num = net.num,
                           #                    pids = obs.pids.vec,
                           #                    bin.table = bin.table,
                           #                    workspace.name = workspace.name1)

                          amongSitesBigDistMat(ssn = ssn,
                                              pids = obs.pids.vec,
                                              net.num = net.num,
                                              bin.table = bin.table,
                                              name = "obs",
                                              workspace.name = workspace.name1)
                         finished1 <- TRUE
                         return(finished1)
                       }
      close(current_distance_matrix)
    }


    ## ---------------------------------------------------------------------------------
    ## Case 3: Calculate distances among prediction sites
    ## ---------------------------------------------------------------------------------
    if (among_predpts & pred.site.no > 0) {
      ## Create distance matrix full of NAs. If the file already exists
      ## and overwrite == FALSE, stop with an error
      ##workspace.name <- paste("dist.net", net.num, sep = "")
      pred.pids <- sort(as.numeric(ssn$preds[[predpts]]$ng.pid[ind.preds]))

      preds.path <- paste0(ssn.object$path, "/distance/", predpts,
                           "/dist.net", net.num)

      net.name <- paste("net", net.num, sep = "")

      ## This uses bin.table from last network if site.no == 0
      if (site.no == 0 | !exists("bin.table")) {
        bin.table <- dbReadTable(connect, net.name)
      }

      among_distance_matrix <- fm.create(filenamebase =preds.path,
                                         nrow = pred.site.no,
                                         ncol = pred.site.no,
                                         type = "double")
      rownames(among_distance_matrix) <- as.character(pred.pids)
      colnames(among_distance_matrix) <- as.character(pred.pids)

      ##among_distance_matrix <- amongSitesDistMat(ssn, pred.pids,
        ##name = predpts, bin.table)

      ## Create vector iterator
      itCol <- isplitVector(pred.pids, chunks = no_cores)

        ans3<- foreach(pred.pids.vec=itCol,
                       .packages = c("SSN2", "filematrix","itertools","iterators"),
                       .errorhandling = "pass") %dopar% {
                           amongSitesBigDistMat(ssn = ssn,
                                                pids = pred.pids.vec,
                                                net.num = net.num,
                                                name = predpts,
                                                bin.table = bin.table,
                                                workspace.name = preds.path)
                      finished3 <- TRUE
                      return(finished3)
                       }
        close(among_distance_matrix)
    }
      ## serialize(among_distance_matrix, file_handle, ascii = FALSE)
  }
}
