#' Helper function to determining distance matrices among sites
#'
#' @param ssn An SSN object.
#' @param pids A list of pid values for prediction sites
#' @param name The network name (obs or prediction name)
#' @param bin.table A binaryID table for the network.
#' @param workspace.name Name of new distance matrix file
#' @param net.num Network number (netID) in character format.
#'
#' @return A distance matrix
#' @export amongSitesBigDistMat
#'
amongSitesBigDistMat <- function(ssn, pids, net.num, name = "obs", bin.table,
                              workspace.name) {
  site.no <- length(pids)
  ##among_distance_matrix <- matrix(NA, nrow = site.no, ncol = site.no)
  among_distance_matrix<- fm.open(filenamebase = workspace.name,
                                  readonly = FALSE)


  # on.exit(
  #     filematrix::close(among_distance_matrix)
  # )
  #diag(among_distance_matrix) <- 0
  # rownames(among_distance_matrix) <- pids
  # colnames(among_distance_matrix) <- pids

  if (name != "obs") {

    all.pids.ind <- ssn$preds[[name]]$NetworkID %in% net.num
    locID.pid.data <- ssn$preds[[name]]$ng.locID[all.pids.ind]
    pid.data <- ssn_get_netgeom(ssn$preds[[name]][all.pids.ind, ], c(
      "pid", "SegmentID", "locID",
      "DistanceUpstream"), reformat = TRUE)

    #pid.data <- as.data.frame(lapply(pid.data, as.numeric))
    colnames(pid.data) <- c("pid", "rid", "locID", "upDist")
  } else {
    all.pids.ind <- ssn$obs$NetworkID %in% net.num
    locID.pid.data <- ssn$obs$locID[all.pids.ind]
    pid.data <- ssn_get_netgeom(ssn$obs[all.pids.ind, ], c(
      "pid", "SegmentID", "locID",
      "DistanceUpstream"), reformat = TRUE)

    ## pid.data <- as.data.frame(sapply(pid.data, as.numeric))
    colnames(pid.data) <- c("pid", "rid", "locID", "upDist")
  }

  pid.data <- pid.data[order(pid.data$pid), ]
  ## New - check if it was already factor
  ##pid.data$locID <- as.factor(pid.data$locID)

  ## Need bin.table
  pid.data$binaryID <- bin.table$binaryID[match(pid.data$rid, bin.table$rid)]
  pid.data <- pid.data[order(pid.data[, "pid"]), ]
  rownames(pid.data) <- pid.data$pid

  ##-----------

  ## New-------
  ob.i_by_locID <- pid.data[order(pid.data[,"locID"]),]
  ob.i_by_locID$pid <- as.numeric(ob.i_by_locID$pid)
  ob.i_by_locID$locID <- as.numeric(ob.i_by_locID$locID)
  ob.j_reordering <- order(ob.i_by_locID$pid)

  ind.dup <- !duplicated(ob.i_by_locID$locID)

  ## locID values can be repeated, in which case they have the same distance data.
  locID.old <- -1

  for (b in seq_len(site.no)) {

    ind.pid <- which(ob.i_by_locID$pid == pids[b])
    pid.b <- ob.i_by_locID[ind.pid, "pid"]
    locID.b <- ob.i_by_locID[ind.pid, "locID"]
    upDist.b <- ob.i_by_locID[ind.pid, "upDist"]

    if (locID.b != locID.old) {

      ## Extract data.frame with columns rid and binaryID
      junk <- get.rid.fc(ob.i_by_locID[ind.dup,"binaryID"],
                               ob.i_by_locID$binaryID[ind.pid])

      ob.j <- getSitesRelationshipsDF(ssn, pid.b,  junk, ind.dup,
                                      name = name,
                                    ob.i_by_locID, bin.table)

      # truncated.binaryIDs <- data.frame(pid = pid.data[, "pid"], junk,
      #                                   stringsAsFactors = FALSE)
      # truncated.binaryIDs$fc <- as.logical(truncated.binaryIDs$fc)
      # truncated.binaryIDs$junc.rid <- bin.table$rid[match(truncated.binaryIDs$binaryID, bin.table$binaryID)]
      #
      # truncated.binaryIDs$juncDist <- ssn$edges$DistanceUpstream[match(
      #   truncated.binaryIDs$junc.rid,
      #   ssn$edges$rid
      # )]
      # truncated.binaryIDs$upDist.j <- pid.data$upDist[match(truncated.binaryIDs$pid, pid.data$pid)]
      ind.fc<-ob.j$fc==1
      ##dist.sites <- ifelse(ind.fc, upDist.b - truncated.binaryIDs$upDist.j,
        ##                   upDist.b - truncated.binaryIDs$juncDist)
      dist.sites <- ifelse(ind.fc, upDist.b - ob.j$upDist.j,
                         upDist.b - ob.j$juncDist)

      col.ind<- colnames(among_distance_matrix) == as.character(pid.b)
      ##among_distance_matrix[, paste(pid.b)] <- ifelse(dist.sites < 0, 0, dist.sites)
      among_distance_matrix[,col.ind] <- ifelse(dist.sites<0, 0, dist.sites)
      locID.old <- locID.b
    } else {
      col.ind <- colnames(among_distance_matrix)== as.character(pid.b)
      among_distance_matrix[, col.ind]<- ifelse(dist.sites<0, 0, dist.sites)
    }
  }

  close(among_distance_matrix)

}
