
## amongPredsBigDistMat(ssn = ssn, net.num = net.num, pids = pred.pids.vec,
##                         pred.num = pred.num, bin.table = bin.table


amongPredsBigDistMat <- function(ssn, net.num, pids, pred.num, bin.table, workspace.name){

    ##require(filematrix)
    pred.site.no <- length(pids)
    all.pids.ind <- ssn@predpoints@SSNPoints[[pred.num]]@point.data$netID %in% net.num

    ## CREATE DISTANCE MATRIX OUTSIDE OF THIS FUNCTION AND THEN ACCESS IT
    among_distance_matrix <- fm.open(filenamebase = workspace.name, readonly = FALSE)

    locID.pid.data <- attributes(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords)$locID[all.pids.ind]
    pid.data <- as.data.frame(cbind(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[all.pids.ind,])),
        as.numeric(levels(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[all.pids.ind]))[ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[all.pids.ind]],
        locID.pid.data, ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$DistanceUpstream[all.pids.ind]))

    colnames(pid.data)<- c("pid","rid", "locID", "upDist")
    pid.data <- pid.data[order(pid.data$pid),]
    pid.data$locID <- as.factor(pid.data$locID)

    pid.data$binaryID <- bin.table$binaryID[match(pid.data$rid, bin.table$rid)]
    pid.data <- pid.data[order(pid.data[,"pid"]),]
    rownames(pid.data) <- pid.data$pid

    ## ------------------------------------------------
    ## New from amongObsBigDistMat---------------------
    ob.i_by_locID <- pid.data[order(pid.data[,"locID"]),]
    ob.i_by_locID$pid <- as.numeric(ob.i_by_locID$pid)
    ob.i_by_locID$locID <- as.numeric(ob.i_by_locID$locID)
    ob.j_reordering <- order(ob.i_by_locID$pid)

    ind.dup <- !duplicated(ob.i_by_locID$locID)
    ##-------------------------------------------------

    ##locID values can be repeated, in which case they have the same distance data.
    locID.old <- -1
    for(b in 1:(length(pids))){

       ## ind.pid <- which(pid.data$pid == pids[b])
       ## locID.b <- pid.data[ind.pid, "locID"]
       ## upDist.b <- pid.data[ind.pid, "upDist"]
       ## pid.b <- pid.data[ind.pid, "pid"]

       ## New from amongObsBigDistMat---------------------
       ind.pid <- which(ob.i_by_locID$pid == pids[b])
       pid.b <- ob.i_by_locID[ind.pid, "pid"]
       locID.b <- ob.i_by_locID[ind.pid, "locID"]
       upDist.b <- ob.i_by_locID[ind.pid, "upDist"]
       ##-------------------------------------------------


        if(locID.b != locID.old) {

           ##junk <- SSN:::get.rid.fc(pid.data[,"binaryID"], pid.data$binaryID[ind.pid])
           junk <- SSN:::get.rid.fc(ob.i_by_locID[ind.dup,"binaryID"], ob.i_by_locID$binaryID[ind.pid])
           ob.j <- getPredsRelationshipsDF(ssn, pid.b,  junk, ind.dup, ob.i_by_locID, bin.table)

           upDist.i <- ssn@predpoints@SSNPoints[[1]]@network.point.coords[paste(pid.b),"DistanceUpstream"]
           ob.j <-ob.j[ob.j_reordering,]

            ##preds fills in by column because it's working between preds to preds.
            ind.fc<-ob.j$fc==1
            dist.obs <- ifelse(ind.fc, upDist.b - ob.j$upDist.j,
                               upDist.b - ob.j$juncDist)
            ##--------------

           ## truncated.binaryIDs <- data.frame(pid=pid.data[,"pid"], junk, stringsAsFactors = FALSE)
           ## truncated.binaryIDs$fc <- as.logical(truncated.binaryIDs$fc)
           ## truncated.binaryIDs$junc.rid <- bin.table$rid[match(truncated.binaryIDs$binaryID, bin.table$binaryID)]

           ## truncated.binaryIDs$juncDist <- ssn@data$upDist[match(truncated.binaryIDs$junc.rid, ssn@data$rid)]
           ## truncated.binaryIDs$upDist.j <- pid.data$upDist[match(truncated.binaryIDs$pid, pid.data$pid)]
           ## ind.fc<-truncated.binaryIDs$fc==1
           ##dist.preds <- ifelse(ind.fc, upDist.b - truncated.binaryIDs$upDist.j,
                           ##upDist.b - truncated.binaryIDs$juncDist)
           dist.preds <- ifelse(ind.fc, upDist.b - ob.j$upDist.j, upDist.b - ob.j$juncDist)

           locID.old <- locID.b

            ## Need to expand column again
            col.ind <- colnames(among_distance_matrix) == as.character(pid.b)
            among_distance_matrix[,col.ind] <- ifelse(dist.preds<0, 0,
                    dist.preds)
        } else {
            ## Add same column to among_distance__matrix
            col.ind <- colnames(among_distance_matrix) == as.character(pid.b)
            among_distance_matrix[,col.ind] <- ifelse(dist.preds<0, 0,
                                                      dist.preds)
        }
    }
    close(among_distance_matrix)
}



