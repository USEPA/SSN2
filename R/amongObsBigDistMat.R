
amongObsBigDistMat <- function(ssn, net.num, pids, bin.table, workspace.name){

    ##require(filematrix)

    pred.num <- 1

    site.no <- length(pids)
    ##all.pids.ind <- ssnobspoints@SSNPoints[[pred.num]]@point.data$netID %in% net.num
    ##all.pids.ind <- ssn$obs$netID %in% net.num

    ## CREATE DISTANCE MATRIX OUTSIDE OF THIS FUNCTION AND THEN ACCESS IT
    among_distance_matrix <- fm.open(filenamebase = workspace.name, readonly = FALSE)
#
#     locID.pid.data <- attributes(ssn@obspoints@SSNPoints[[pred.num]]@network.point.coords)$locID[all.pids.ind]
#     pid.data <- as.data.frame(cbind(as.numeric(rownames(ssn@obspoints@SSNPoints[[pred.num]]@network.point.coords[all.pids.ind,])),
#         as.numeric(levels(ssn@obspoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[all.pids.ind]))[ssn@obspoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[all.pids.ind]],
#         locID.pid.data, ssn@obspoints@SSNPoints[[pred.num]]@network.point.coords$DistanceUpstream[all.pids.ind]))

    ind.pids <- ssn$obs$ng.pid %in% as.character(pids)
    locID.pid.data <- ssn$obs$locID[ind.pids]
    pid.data <- ssn_get_netgeom(ssn$obs[ind.pids, ], c(
      "pid", "SegmentID", "locID",
      "DistanceUpstream"), reformat = TRUE)

    colnames(pid.data)<- c("pid","rid", "locID", "upDist")
    pid.data <- pid.data[order(pid.data$pid),]
    pid.data$locID <- as.factor(pid.data$locID) ## New

    pid.data$binaryID <- bin.table$binaryID[match(pid.data$rid, bin.table$rid)]
    pid.data <- pid.data[order(pid.data[,"pid"]),]
    rownames(pid.data) <- pid.data$pid

    ## New-------
    ob.i_by_locID <- pid.data[order(pid.data[,"locID"]),]
    ob.i_by_locID$pid <- as.numeric(ob.i_by_locID$pid)
    ob.i_by_locID$locID <- as.numeric(ob.i_by_locID$locID)
    ob.j_reordering <- order(ob.i_by_locID$pid)

    ind.dup <- !duplicated(ob.i_by_locID$locID)
    ##-----------

    ##locID values can be repeated, in which case they have the same distance data.
    locID.old <- -1

    ## for(b in 1:(nrow(ob.i_by_locID))){
    for(b in 1:length(pids)){

       ind.pid <- which(ob.i_by_locID$pid == pids[b])
       pid.b <- ob.i_by_locID[ind.pid, "pid"]
       locID.b <- ob.i_by_locID[ind.pid, "locID"]
       upDist.b <- ob.i_by_locID[ind.pid, "upDist"]


        if(locID.b != locID.old) {
            junk <- get.rid.fc(ob.i_by_locID[ind.dup,"binaryID"], ob.i_by_locID$binaryID[ind.pid])
            ob.j <- getObsRelationshipsDF(ssn, pid.b,  junk, ind.dup, ob.i_by_locID, bin.table)
	    ##upDist.i <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[paste(pid.b),"DistanceUpstream"]
	    upDist.i <- ssn$obs[paste(pid.b),"DistanceUpstream"]
            ob.j <-ob.j[ob.j_reordering,]

            ##obs fills in by column because it's working between obs to obs.
            ind.fc<-ob.j$fc==1
            dist.obs <- ifelse(ind.fc, upDist.b - ob.j$upDist.j,
                                 upDist.b - ob.j$juncDist)
           ## truncated.binaryIDs <- data.frame(pid=pid.data[,"pid"], junk, stringsAsFactors = FALSE)
           ## truncated.binaryIDs$fc <- as.logical(truncated.binaryIDs$fc)
           ## truncated.binaryIDs$junc.rid <- bin.table$rid[match(truncated.binaryIDs$binaryID, bin.table$binaryID)]

           ## truncated.binaryIDs$juncDist <- ssn@data$upDist[match(truncated.binaryIDs$junc.rid, ssn@data$rid)]
           ## truncated.binaryIDs$upDist.j <- pid.data$upDist[match(truncated.binaryIDs$pid, pid.data$pid)]
           ## ind.fc<-truncated.binaryIDs$fc==1
           dist.preds <- ifelse(ind.fc, upDist.b - ob.j$upDist.j, upDist.b - ob.j$juncDist)

            ## Need to expand column again

           col.ind <- colnames(among_distance_matrix) == as.character(pid.b)
           among_distance_matrix[,col.ind] <- ifelse(dist.preds<0, 0,
                                                      dist.preds)

          locID.old <- locID.b
        } else {
            ## Add same column to among_distance__matrix
            col.ind <- colnames(among_distance_matrix) == as.character(pid.b)
            among_distance_matrix[,col.ind] <- ifelse(dist.preds<0, 0,
                                                      dist.preds)
        }
    }
    close(among_distance_matrix)
}
