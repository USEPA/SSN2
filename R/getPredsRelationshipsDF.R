

getPredsRelationshipsDF <- function(ssn, pid, junk, ind, ob_by_locID, bin) {

    ## ssn = SpatialStreamNetwork object
    ## pid = integer pid value of interest
    ## junk = ox2 data.frame with columns fc (logical) and binaryID
    ## ind = vector indicator saying whether it is a duplicate
    ## ob = data.frame with pid, rid, locID, and binaryID for sites on network
    ##      ordered by pid
    ## ob_by_locID = data.frame with pid, rid, locID, and binaryID for sites on network
    ##               ordered by locID
    ## bin = binaryID table

    ## Returns a data.frame that relates all sites to pid.i:
    ## pid: numeric
    ## locID: numeric
    ## fc: logical - is the sute fc with the pid of interest
    ## binaryID: binaryID of the common downstream junction
    ## junc.rid: rid for the common downstream junction
    ## upDist.j: upDist for each site

## Create relationships table

      ob.j.r <- data.frame(ob_by_locID[ind, c("pid", "locID")], junk,
                           stringsAsFactors = FALSE)

      ## ob.j.r <- data.frame(ob_by_locID[ind, c("pid", "locID")], junk,
      ##                      stringsAsFactors = FALSE)

      ob.j.r$fc <- as.logical(ob.j.r$fc)
      rownames(ob.j.r)<- ob.j.r$pid

      ## Add column showing rid for common downstream junction (ob.j.r relates all ##sites to pid)
      ## dataframe with pid, locID, fc, binaryID, junc.rid
      ob.j.r$junc.rid <- bin$rid[match(ob.j.r$binaryID, bin$binaryID)]

      ## Expand this data.frame ---
      ## Add extra rows
      ob.j.r[,"rep"] <- aggregate(rep(1, length(ob_by_locID$locID)), by = list(ob_by_locID$locID), sum)[2]
      ob.j <- ob.j.r[rep(row.names(ob.j.r), ob.j.r$rep),1:ncol(ob.j.r)-1]
      ## fix pid values
      ob.j$pid <- ob_by_locID$pid
      rownames(ob.j)<- rownames(ob_by_locID)

      ## Reorder back to same as pid.data

      ## Create some funky rownames, with extension .fc - ADDED OB.J INSTEAD OF OB
      rownames(ob.j) <- paste(rownames(ob.j), ".fc", sep = "")

      ## Don't know why we're doing this...
      ##ob.j$pid <- ob_by_locID$pid[ind]

      ## juncDist is the upstream distance of the common downstream rid junction
      ob.j$juncDist <- ssn@network.line.coords$DistanceUpstream[match(ob.j$junc.rid, ssn@network.line.coords$SegmentID)]

      ## upDist.j is the upDist for each observed site
      ob.j$upDist.j <- ssn@predpoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream[
                    match(ob.j$pid, as.numeric(rownames(ssn@predpoints@SSNPoints[[1]]@network.point.coords)))]

      ob.j

}
