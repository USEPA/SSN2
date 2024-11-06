#' Calculate Hydrologic Distances for a large SSN` Object
#'
#' @description Creates a collection of (non-symmetric) matrices containing pairwise downstream hydrologic distances between sites in an 'SSN' object
#'
#' @param ssn 'SSN' object
#' @param predpts a valid predpoints ID from the SSN. Default is \code{NULL}.
#' @param o.write If \code{TRUE}, overwrite existing distance matrices. Defaults to \code{FALSE}.
#' @param amongpreds If \code{TRUE}, compute the distances between the prediction sites. Defaults to \code{FALSE}.
#' @param no.cores Number of cores to use in computation of the distance matrices. Also, the number of chunks to split the dataset into during computation.
#'
#' @details   A distance matrix that contains the hydrologic distance between any two sites in SpatialStreamNetwork object is needed to fit a spatial statistical model using the tail-up and tail-down autocovariance functions described in Ver Hoef and Peterson (2010). These models are implemented in R via \command{glmssn} in the SSN package. The hydrologic distance information needed to model the covariance between flow-connected (i.e. water flows from one location to the other) and flow-unconnected (i.e. water does not flow from one location to the other, but they reside on the same network) locations differs. The total hydrologic distance is a directionless measure; it represents the hydrologic distance between two sites, ignoring flow direction. The hydrologic distance   from each site to a common downstream stream junction is used when creating models for flow-unconnected pairs, which we term downstream hydrologic distance. In contrast, the total hydrologic distance is used for modeling flow-connected pairs, which we term total hydrologic distance. \cr
#' \cr
#' A downstream hydrologic distance matrix provides enough information to meet the data requirements for both the tail-up and tail-down models. When two locations are flow-connected, the downstream hydrologic distance from the upstream location to the downstream location is greater than zero, but it is zero in the other direction. When two locations are flow-unconnected the downstream hydrologic distance will be greater than zero in both directions. A site's downstream hydrologic distance to itself is equal to zero. The format of the downstream hydrologic distance matrix is efficient because distance information needed to fit both the tail-up and tail-down models is only stored once. As an example, a matrix containing the total hydrologic distance between sites is easily calculated by adding the downstream distance matrix to its transpose. \cr
#' \cr
#' The downstream hydrologic distances are calculated based on the binaryIDs and stored as matrices. The matrices are stored in a directory named \sQuote{distance}, which is created by the createBigDistMat function within the .ssn directory. The distance directory will always contain at least one directory named \sQuote{obs}, which contains a number of files for each network that has observed sites residing on it with file extensions .bmat, .desc.txt, nmscol.txt, and .nmsrow.txt. The basefile naming convention is based on the netID number (e.g. dist.net1). Each matrix in the \sQuote{obs} folder contains the information to form a square matrix, which contains the downstream hydrologic distance between each pair of observed sites on the network. Direction is preserved, with columns representing the FROM site and rows representing the TO site. Row and column names correspond to the pid attribute for each site. \cr
#' \cr
#' If the argument \code{predpts} is specified in the call to the function, the downstream hydrologic distances between the observed and prediction sites will also be computed. A new directory is created within the distance directory, with the name corresponding to the predpoints ID (e.g. \dQuote{preds}). A sequence of files is created within this directory, similar to the structure for the observed sites, except that two objects are stored for each network that contains \emph{both} observed and prediction sites. The letters \code{a} and \code{b} are used in the naming convention to distinguish between the two objects (e.g. dist.net1.a and dist.net1.b). The matrices that these objects represent are not necessarily square. In matrices of type \code{a} and \code{b}, rows correspond to prediction locations and columns to observed locations. However, direction is preserved differently in these matrices. For matrices of type \code{a}, columns represent the TO site and rows represent the FROM site. In contrast, the columns represent hte FROM site and the rows represent the TO site in matrices of type \code{b}. Again, row and column names correspond to the pid attribute for each site. \cr
#' \cr
#'   If the argument \code{amongpreds} is set to TRUE, the downstream hydrologic distances will also be computed between prediction sites, for each network. Again these are stored within the distance directory with the name corresponding to the predpoints ID. The naming convention for these prediction to prediction site distance matrices is the same as the distance matrices stored in the \sQuote{obs} directory (e.g. dist.net1). These extra distance matrices are needed to perform block Kriging using the \command{glmssn}. \cr
#' \cr
#' The distance matrices generated using the createBigDistMat function are file-backed, meaning that they are stored as files, rather than in memory. This approach is more computationally efficient when users must calculate a large number of pairwise distances between observed and/or prediction locations. The \code{no.cores} argument designates the number of cores to be used if the code is to be executed in parallel. Note that, the \code{detectCores} function can be used to detect the number of cores available. The \code{no.cores} argument is also used to 'chunk' the dataset into multiple parts for processing, with one part allocated to each core. However, there is an overhead for connecting to each core and so the computational benefits of using the createBigDistMat function over the createDistMat function are only evident when datasets are large.
#'
#' @return The \command{createBigDistMat} function creates a collection of hierarchical directories in the \code{ssn$path} directory, which store the pairwise distances between sites associated with the \link{SpatialStreamNetwork-class} object. See details section for additional information.
#'
#'
#' @author Erin Peterson
#' @export

createBigDistMat <- function(ssn, predpts = NULL, o.write = FALSE, amongpreds = FALSE,
                             no.cores = 1) {

    ## ssn = SpatialStreamNetwork object
    ## predpts = name of prediction points in character format
    ## o.write = should existing distance directories be overwritten
    ## amongpreds = should the pred-pred distance matrix be generated
    ## no.cores = number of cores to use and number of splits in dataset


    ## require(filematrix)
    ## require(foreach)
    ## require(doParallel)
    ## require(itertools)


  #############################################################################
  ## Check that arguments are valid and set up directories
  if(amongpreds && (missing(predpts) || is.null(predpts)))
  {
	stop("A named collection of prediction points must be specified via the predpts option when amongpreds is TRUE")
  }
  ##Check to see whether distance folder exists...
  if (!file.exists(file.path(ssn$path, "distance"))) {
    dir.create(file.path(ssn$path, "distance"))
  }

  ##Check whether an observation folder exists
  obs.path <- file.path(ssn$path, "distance", "obs")
  if (!file.exists(obs.path)) {
    dir.create(obs.path)
  }

  ## And then whether prediction folder exists
  if (!is.null(predpts)) {
    preds.path <- file.path(ssn$path, "distance", predpts)
    if(!file.exists(preds.path)) {
      dir.create(preds.path)
    }

    ## Get ID number from predictions
    count <- 0
    if(length(ssn$preds) > 0) {
	    for (m in 1:length(ssn$preds)) {
	      if (names(ssn$preds)[m] == predpts) {
	           pred.num <- m
             count <- count + 1
        }
	    }
    }

    if (count==0) {
      stop(predpts, " does not exist in SSN")}

    if (count > 1) {
      stop("SSN contains more than one copy of ", predpts)}

    ssn@predpoints@SSNPoints[[pred.num]]@point.data$netID<- as.factor(ssn@predpoints@SSNPoints[[pred.num]]@point.data$netID)
  }

  if (is.null(predpts)) {
      pred.num <- 0}


  ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID<-
      as.factor(ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID)
  net.count <- length(levels(ssn@network.line.coords$NetworkID))
  warned.overwrite <- FALSE

      if (file.exists(file.path(ssn$path,"binaryID.db")) == FALSE)
        stop("binaryID.db is missing from ssn object")

      driver <- RSQLite::SQLite()
      connect.name <- file.path(ssn$path,"binaryID.db")

  connect <- dbConnect(SQLite(), connect.name)


      cl <- makeCluster(no.cores)
      registerDoParallel(cl)

      ## close connections upon function exit
      on.exit({
          dbDisconnect(connect)
          closeAllConnections()
      })
  #######################################################################
  ### FOR EACH NETWORK---------------------------------------------------
  for (i in 1:net.count) {
      net.num <- levels(ssn@network.line.coords$NetworkID)[i]
      ind.obs <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == as.numeric(net.num)

      ## figure out how many observed sites there are in the network
      site.no <- nrow(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])

      if(pred.num > 0) {
          ind.preds <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[,"NetworkID"] %in% net.num
          pred.site.no <- nrow(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])
      } else {
          pred.site.no <- 0
      }

      net.name <- paste("net", net.num, sep = "")


      ##################################################################
      ## If observed sites > 0 #########################################
      if (site.no > 0) {

          workspace.name1 <- paste(obs.path, "/dist.net", net.num, sep = "")

          ## get sorted pids to use as row and column names
          obs.pids<- sort(as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])))

          ## Extract binaryID table
          bin.table <- dbReadTable(connect, net.name)

          exists <- file.exists(paste0(workspace.name1, ".bmat"))

          ###################################################################
          ## If o.write = FALSE and obs file exists
	  if(!o.write & exists) {
                   cat("Distance matrices for observed versus observed sites in network", net.num, "already exists and o.write = FALSE. Not overwriting existing matrices\n\n")

          } else {

          ######################################################################

              ## Create empty obs distance matrix
              current_distance_matrix <-
                    fm.create(filenamebase = workspace.name1,
                              nrow = site.no, ncol = site.no, type = "double")


               rownames(current_distance_matrix) <- as.character(obs.pids)
               colnames(current_distance_matrix) <- as.character(obs.pids)

              ##Calculate distances between observed sites
              if (is.null(getDoParName())) {
                  registerDoSEQ() # A little hack to avoid the foreach warning 1st time.
              }

              ##cl <- makeCluster(no.cores)
              ##registerDoParallel(cl)
              ## registerDoParallel(cores=no.cores)

              ## Create vector iterator
              itCol <- isplitVector(obs.pids, chunks = no.cores)

              ans1<- foreach(obs.pids.vec=itCol,.packages = c("SSN", "filematrix","itertools","iterators"),
                                 .errorhandling = "pass") %dopar% {

                          amongObsBigDistMat(ssn = ssn, net.num = net.num,
                                             pids = obs.pids.vec,
                                             bin.table = bin.table,
                                             workspace.name = workspace.name1)
                          finished1 <- TRUE
                          return(finished1)
                      }

              ## stopCluster(cl)
              ## registerDoSEQ()
          }
      } ## END IF SITE.NO > 0

##################################################################################
      ## Calculate Observed to prediction site distances
      if(site.no>0 & pred.site.no > 0) {
          ## figure out how many prediction sites there are in the network
         ## ind.preds <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[,"NetworkID"] %in% net.num
         ## pred.site.no <- nrow(ssn@obspoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])

         #if(pred.site.no > 0) {

             ## Get pred pids
         pred.pids<- sort(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])))

         ## create empty distance matrices for obs-preds
         workspace.name.a <- paste(preds.path, "/dist.net", net.num, ".a", sep = "")
         workspace.name.b <- paste(preds.path, "/dist.net", net.num, ".b", sep = "")

         exists.a <- file.exists(paste0(workspace.name.a, ".bmat"))
         exists.b <- file.exists(paste0(workspace.name.b, ".bmat"))

###################################################################
         ## If o.write = FALSE and obs file exists
	 if(!o.write & all(exists.a, exists.b)) {
             cat("Distance matrices for observed versus predicted sites in network",
                 net.num, " already exist and o.write = FALSE.Not overwriting existing matrices\n\n")

          } else {
#####################

             current_distance_matrix_a <-
                 fm.create(filenamebase = paste0(workspace.name.a),
                       nrow = pred.site.no, ncol = site.no, type = "double")

             rownames(current_distance_matrix_a) <- as.character(pred.pids)
             colnames(current_distance_matrix_a) <- as.character(obs.pids)

             current_distance_matrix_b <-
                    fm.create(filenamebase = paste0(workspace.name.b),
                              nrow = pred.site.no, ncol = site.no, type = "double")
             rownames(current_distance_matrix_b) <- as.character(pred.pids)
             colnames(current_distance_matrix_b) <- as.character(obs.pids)


             ## Create vector iterator
             itCol <- isplitVector(obs.pids, chunks = no.cores)

             ## cl <- makeCluster(no.cores)
             ## registerDoParallel(cl)

             ans2<- foreach(obs.pids.vec=itCol,.packages = c("SSN", "filematrix","itertools","iterators"),
                             .errorhandling = "pass") %dopar% {

                      amongObsPredsBigDistMat(ssn = ssn, net.num = net.num, pred.num, obs.pids = obs.pids.vec,
                                                    pred.pids=pred.pids, bin.table = bin.table, workspace.name.a = workspace.name.a,
                                                    workspace.name.b = workspace.name.b)
                      finished2 <- TRUE
                      return(finished2)
             }

             ## stopCluster(cl)
             ## registerDoSEQ()
         }
      }## END OBS-PRED SITE DISTANCES

##------------------------------------------------------
## ###############################################################################################

          ## Calculate distances between prediction sites ------------------------------------
      if (amongpreds & pred.num > 0) {

           pred.pids<- sort(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])))
           net.name <- paste("net", net.num, sep = "")

           among.path <- paste0(preds.path, "/dist.net", net.num)

           exists2 <- file.exists(paste0(among.path, ".bmat"))

###################################################################
           ## If o.write = FALSE and preds versus preds file exists
	   if(!o.write & exists2) {
             cat("Distance matrices for predicted versus predicted sites in network",
                 net.num, " already exist and o.write = FALSE. Not overwriting existing matrices\n\n")

           } else {
#####################


               among_distance_matrix <- fm.create(filenamebase =among.path,
                     nrow = pred.site.no, ncol = pred.site.no, type = "double")
               rownames(among_distance_matrix) <- as.character(pred.pids)
               colnames(among_distance_matrix) <- as.character(pred.pids)

               bin.table <- dbReadTable(connect, net.name)
               ##Calculate distances between observed sites
               ## cl <- makeCluster(no.cores)
               ## registerDoParallel(cl)
          ##registerDoParallel(cores=no.cores)

               ## Create vector iterator
               itCol <- isplitVector(pred.pids, chunks = no.cores)

               ans3<- foreach(pred.pids.vec=itCol,.packages = c("SSN", "filematrix","itertools","iterators"),
                              .errorhandling = "pass") %dopar% {

                      amongPredsBigDistMat(ssn = ssn, net.num = net.num, pids = pred.pids.vec,
                            pred.num = pred.num, bin.table = bin.table,
                                           workspace.name = among.path)
                      finished3 <- TRUE
                      return(finished3)
               }
               ## stopCluster(cl)
               ## registerDoSEQ()

           }
      }## End among preds
  }
  stopCluster(cl)
  registerDoSEQ()

}
