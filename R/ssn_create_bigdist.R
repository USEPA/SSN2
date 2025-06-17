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
#' @param verbose Logical. If \code{TRUE}, warning messages (if they exist) are printed to the console. Defaults to \code{TRUE}.
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
                               no_cores = 1, verbose = TRUE) {

  # ##iterate if predpts is not null
  # if (length(predpts) > 1) {
  #   if (only_predpts) {
  #     x <- lapply(predpts, function(x) {
  #       ssn_create_bigdist(ssn.object,
  #         predpts = x,
  #         overwrite,
  #         among_predpts,
  #         only_predpts,
  #         no_cores
  #       )
  #     })
  #   } else {
  #     x1 <- ssn_create_bigdist(ssn.object, overwrite = overwrite)
  #     x2 <- lapply(predpts, function(x) {
  #       ssn_create_bigdist(ssn.object,
  #         predpts = x,
  #         overwrite,
  #         among_predpts,
  #         only_predpts = TRUE,
  #         no_cores
  #       )
  #     })
  #   }
  #   return(invisible(NULL))
  # }

  # recreate function argument name
  ssn <- ssn.object

  obs.pids.vec <- NULL
  pred.pids.vec <- NULL

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
    if (!file.exists(file.path(ssn$path, "distance", "obs")) &
        verbose == TRUE) {
      warning("only_predpts == TRUE and distance matrices for observed sites are missing. Set only_predpts = FALSE to calculate distance matrices for observed and prediction sites.\n",
              call. = FALSE)
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
      suppressWarnings(closeAllConnections())
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

    ## Create observed distance matrices if necessary
    ## ----------------------------------------------
    if(site.no > 0) {
      if(only_predpts == FALSE) {

        ## Set o x o distance matrix name
        workspace.name1 <- paste(ssn$path, "/distance/obs/dist.net",
                                 net.num, sep = "")

        obs.pids <- sort(as.numeric(ssn$obs$ng.pid[ind.obs]))

        ## Check for o x o
        obs.exist <- file.exists(paste0(workspace.name1, ".bmat"))

        ## Create empty obs distance matrix
        if((overwrite == TRUE) | obs.exist == FALSE) {

          current_distance_matrix <- fm.create(
            filenamebase = workspace.name1,
            nrow = site.no,
            ncol = site.no,
            type = "double")

          rownames(current_distance_matrix) <- as.character(obs.pids)
          colnames(current_distance_matrix) <- as.character(obs.pids)

          close(current_distance_matrix)
          obs.exist <- FALSE

        } else {
          if(verbose == TRUE){
            message(paste0("observed site distance matrix for netID = ",
                         net.num, " exists and overwrite = FALSE. No changes made to distance matrix.\n"))
          }
        }
      }
    }

    ## -------------------------------------------------------------
    ## CASE 1: IF OBS and PREDS EXIST ON NETWORK
    ## -------------------------------------------------------------
    if (site.no > 0 & pred.site.no > 0) {
      ## get sorted pids to use as dim names
      ##obs.pids <- sort(as.numeric(ssn$obs$ng.pid[ind.obs]))

      ## Get pred pids and sort
      pred.pids <- sort(as.numeric(ssn$preds[[predpts]]$ng.pid[ind.preds]))

      ## Set distance matrix names for o x p
      workspace.name.a <- paste(ssn$path, "/distance/", predpts,
                                "/dist.net", net.num, ".a", sep = "")
      workspace.name.b <- paste(ssn$path, "/distance/", predpts,
                                "/dist.net", net.num, ".b", sep = "")

      ## Check for o x p distance matrics
      preds.exist <- sum(file.exists(paste0(workspace.name.a, ".bmat")),
                         file.exists(paste0(workspace.name.b, ".bmat")))

      ## Return error if only one obs x preds distance matrix is present
      if(preds.exist == 1 & overwrite == FALSE) {
        stop(paste0("Big distance matrix files for ", predpts,
                    " are incomplete. Set overwrite = TRUE to recreate observed versus prediction distance matrices."))
      }

      ## Create empty distance matrices if necessary
      if(preds.exist == 0 | overwrite == TRUE) {

        current_distance_matrix_a <-
          fm.create(filenamebase = workspace.name.a,
                    nrow = pred.site.no,
                    ncol = site.no,
                    type = "double")

        rownames(current_distance_matrix_a) <- as.character(pred.pids)
        colnames(current_distance_matrix_a) <- as.character(obs.pids)
        close(current_distance_matrix_a)

        current_distance_matrix_b <-
          fm.create(filenamebase = workspace.name.b,
                    nrow = pred.site.no, ncol = site.no,
                    type = "double")
        rownames(current_distance_matrix_b) <- as.character(pred.pids)
        colnames(current_distance_matrix_b) <- as.character(obs.pids)
        close(current_distance_matrix_b)

        preds.exist <- FALSE

      } else {
        if(verbose == TRUE){
          message(paste0(predpts, " distance matrices for netID = ",
                       net.num,
                       " exist and overwrite = FALSE. No changes made to distance matrices.\n"))
        preds.exist <- TRUE
        }
      }

      ## Extract binaryID table for the network
      bin.table <- dbReadTable(connect, net.name)

      ## Calculate obs x obs distance matrix
      if (!only_predpts) {

        ## Calculate obs x obs distances
        if(obs.exist == FALSE) {

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

        }
      }

      if(preds.exist == FALSE) {
        # Calculate observed versus prediction distance matrices
        parallel::clusterEvalQ(cl, {
          library(SSN2)
          library(filematrix)
        })
        parLapply(cl, obs.pids,
                  function(x) amongObsPredsBigDistMat(ssn = ssn,
                                                      obs.pids = x,
                                                      pred.pids = pred.pids,
                                                      bin.table = bin.table,
                                                      pred.name = predpts,
                                                      workspace.name.a = workspace.name.a,
                                                      workspace.name.b = workspace.name.b))
      }

    }

    ## ---------------------------------------------------------------------------------
    ## Case 2: Observed sites, no prediction sites
    ## ---------------------------------------------------------------------------------
    if (site.no > 0 & pred.site.no == 0 &
        only_predpts == FALSE & obs.exist == FALSE) {

      ## Open binaryID.db
      bin.table <- dbReadTable(connect, net.name)

      # Create vector iterator
      itCol <- isplitVector(obs.pids, chunks = no_cores)

      ans1<- foreach(obs.pids.vec=itCol,.packages = c("SSN2", "filematrix",
                                                        "itertools","iterators"),
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

    }

    ## ---------------------------------------------------------------------------------
    ## Case 3: Calculate distances among prediction sites
    ## ---------------------------------------------------------------------------------
    if (among_predpts & pred.site.no > 0) {

      pred.pids <- sort(as.numeric(ssn$preds[[predpts]]$ng.pid[ind.preds]))
      preds.path <- paste0(ssn.object$path, "/distance/", predpts,
                           "/dist.net", net.num)

      ## Create preds x preds distance matrix if necessary
      ## -------------------------------------------------
      preds.all.exist <- file.exists(paste0(preds.path, ".bmat"))

      if(overwrite == TRUE | preds.all.exist == FALSE) {

        among_distance_matrix <- fm.create(filenamebase =preds.path,
                                           nrow = pred.site.no,
                                           ncol = pred.site.no,
                                           type = "double")
        rownames(among_distance_matrix) <- as.character(pred.pids)
        colnames(among_distance_matrix) <- as.character(pred.pids)

        close(among_distance_matrix)

        preds.all.exist <- FALSE

      }

      if(preds.all.exist == FALSE) {

        ## Open bin.table if necessary
        if (site.no == 0 | !exists("bin.table")) {
          bin.table <- dbReadTable(connect, net.name)
        }

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
      } else {
        if(verbose == TRUE){
        message(paste0("among_predpts distance matrix exists for ", predpts,
        " netID = ", net.num,
        " and overwrite = FALSE. No changes made to distance matrix.\n"))
        }
      }
    }
    ## --------------------------------------------------------

  } ## Close for each network
}
