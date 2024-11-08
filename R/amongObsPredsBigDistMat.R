

amongObsPredsBigDistMat <- function(ssn, obs.pids, pred.pids, bin.table,
                                    workspace.name.a, workspace.name.b,
                                    pred.name){

    ##browser()
    site.no <- length(obs.pids)
    pred.site.no <- length(pred.pids)

    ## Extract and add netgeom data to preds
    pred.data<- ssn_get_netgeom(ssn$preds[[pred.name]],
                                c("pid", "SegmentID", "locID",
                                  "DistanceUpstream"),
                                reformat = TRUE)

    pred.pid.ind<- pred.data$pid %in% pred.pids
    pred.data<- pred.data[pred.pid.ind,]
    ##colnames(pred.data)<- c("netID", "rid", "upDist", "ratio", "pid", "locID")
    colnames(pred.data) <- c("pid", "rid", "locID","upDist")


    # tmp.df <- ssn_get_data(ssn, pred.name)
    # n.geom <- ssn_get_netgeom(ssn$preds[[pred.name]], reformat = TRUE)
    # colnames(n.geom)[4:6] <- paste0("ng.", colnames(n.geom[4:6]))
    # # n.geom$NetworkID <- as.factor(n.geom$NetworkID)
    # # n.geom$DistanceUpstream <- as.numeric(n.geom$DistanceUpstream)
    # tmp.df <- cbind(tmp.df, n.geom)
    # ssn <- ssn_put_data(tmp.df, ssn, pred.name)
    # rm(tmp.df, n.geom)

    ## Extract and add netgeom data to obs
    obs.data<- ssn_get_netgeom(ssn$obs,
                               c("pid", "SegmentID", "locID",
                                "DistanceUpstream"),
                                    reformat = TRUE)
    obs.pid.ind<- obs.data$pid %in% obs.pids
    obs.data<- obs.data[obs.pid.ind,]
    colnames(obs.data)<- c("pid", "rid", "locID","upDist")


    ## CREATE DISTANCE MATRIX OUTSIDE OF THIS FUNCTION AND THEN ACCESS IT here
    current_distance_matrix_a <- fm.open(filenamebase = workspace.name.a,
                                         readonly = FALSE)
    current_distance_matrix_b <- fm.open(filenamebase = workspace.name.b,
                                         readonly = FALSE)

    # on.exit(
    #   filematrix::close(current_distance_matrix_a)
    #   filematrix::close(current_distance_matrix_b)
    # )

    ## Create a data.frame of prediction point data
    locID.pred.data <- pred.data$locID

    pred.data <- pred.data[order(pred.data$pid),]
    pred.data$locID <- as.factor(pred.data$locID)

    pred.data$binaryID <- bin.table$binaryID[match(pred.data$rid, bin.table$rid)]
    pred.data <- pred.data[order(pred.data[,"pid"]),]
    rownames(pred.data) <- pred.data$pid

    ##
    pred_by_locID <- pred.data[order(pred.data[,"locID"]),]
    pred_by_locID$pid <- as.numeric(pred_by_locID$pid)
    pred_by_locID$locID <- as.numeric(pred_by_locID$locID)
    pred_reordering <- order(pred_by_locID$pid)

    dup.pred <- !duplicated(pred_by_locID$locID)
    ##-----------

     ## Create a data.frame of observed point data
    locID.obs.data <- obs.data$locID
    obs.data <- obs.data[order(obs.data$pid),]
    obs.data$locID <- as.factor(obs.data$locID)

    obs.data$binaryID <- bin.table$binaryID[match(obs.data$rid, bin.table$rid)]
    obs.data <- obs.data[order(obs.data[,"pid"]),]
    rownames(obs.data) <- obs.data$pid

    ##
    obs_by_locID <- obs.data[order(obs.data[,"locID"]),]
    obs_by_locID$pid <- as.numeric(obs_by_locID$pid)
    obs_by_locID$locID <- as.numeric(obs_by_locID$locID)
    obs_reordering <- order(obs_by_locID$pid)

    ##-----------

    ##locID values can be repeated, in which case they have the same distance data.
    locID.old <- -1
    for(ob in 1:(nrow(obs_by_locID))){

       ind.pid <- which(obs_by_locID$pid == obs.pids[ob])
       pid.ob <- obs_by_locID[ind.pid, "pid"]
       locID.ob <- obs_by_locID[ind.pid, "locID"]
       upDist.ob <- obs_by_locID[ind.pid, "upDist"]


        if(locID.ob != locID.old) {

            junk <- get.rid.fc(pred_by_locID[dup.pred,"binaryID"],
                                     obs_by_locID$binaryID[ind.pid])


            ob.j <- getObsPredsRelationshipsDF(ssn, junk, dup.pred,
                                              pred_by_locID,
                                              obs_by_locID[ind.pid,],
                                              bin.table, pred.name)

            ob.j <-ob.j[pred_reordering,]

            ##obs fills in by column because it's working between obs to obs.
            ind.fc<-ob.j$fc==1

            dist.a <- ifelse(ind.fc, ob.j$upDist.j-upDist.ob,
                             ob.j$upDist.j - ob.j$juncDist)
            ## WE CHANGE THIS SO THAT WE'RE ADDING COLUMNS (TRANSPOSE IT LATER FOR USE)
            col.ind <- colnames(current_distance_matrix_a) == as.character(pid.ob)
            current_distance_matrix_a[,col.ind] <- ifelse(dist.a<0, 0, dist.a)

            dist.b <- ifelse(ind.fc, upDist.ob - ob.j$upDist.j,
                           upDist.ob - ob.j$juncDist)
            current_distance_matrix_b[,col.ind] <- ifelse(dist.b<0, 0, dist.b)

            locID.old <- locID.ob

        } else {
            col.ind <- colnames(current_distance_matrix_a) == as.character(pid.ob)
            current_distance_matrix_a[,col.ind]<- ifelse(dist.a<0, 0, dist.a)

            current_distance_matrix_b[,col.ind] <- ifelse(dist.b<0, 0, dist.b)
            ##locID.old <- locID.ob
        }

    }
    ##browser()
    close(current_distance_matrix_a)
    close(current_distance_matrix_b)
}
