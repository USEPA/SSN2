#' Get spatial relationships for observed points
#'
#' @param ssn SSN object
#' @param pid integer pid value of interest
#' @param junk data.frame with columns fc (logical) and binaryID
#' @param ind vector indicator saying whether it is a duplicate
#' @param ob data.frame with pid, rid, locID, and binaryID for sites on network ordered by pid
#' @param ob_by_locID data.frame with pid, rid, locID, and binaryID for sites on network
#'   ordered by locID factor level (not true locID)
#' @param bin binaryID table
#'
#' @noRd
getObsRelationshipsDF <- function(ssn, pid, junk, ind, ob, ob_by_locID, bin) {
  ## ssn = SpatialStreamNetwork object
  ## pid = integer pid value of interest
  ## junk = data.frame with columns fc (logical) and binaryID
  ## ind = vector indicator saying whether it is a duplicate
  ## ob = data.frame with pid, rid, locID, and binaryID for sites on network
  ##      ordered by pid
  ## ob_by_locID = data.frame with pid, rid, locID, and binaryID for sites on network
  ##               ordered by locID factor level (not true locID)
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
    stringsAsFactors = FALSE
  )

  ob.j.r$fc <- as.logical(ob.j.r$fc)
  rownames(ob.j.r) <- ob.j.r$pid

  ## Add column showing rid for common downstream junction (ob.j.r relates all
  ## sites to pid)
  ob.j.r$junc.rid <- bin$rid[match(ob.j.r$binaryID, bin$binaryID)]

  ## Expand ob.j.r to include repeated measurements
  reps <- as.numeric(ob_by_locID$locID)
  ob.j <- ob.j.r[reps, ]
  ## ob.j <- ob.j.r

  ## Create rownames, with extension .fc
  rownames(ob.j) <- paste(rownames(ob), ".fc", sep = "")

  ## Reassign pid after expanding repeated measurements
  ob.j$pid <- ob_by_locID$pid

  ## juncDist is the upstream distance of the common downstream rid junction
  ob.j$juncDist <- ssn$edges$DistanceUpstream[match(
    ob.j$junc.rid,
    ssn$edges$SegmentID
  )]

  ## upDist.j is the upDist for each observed site
  ob.j$upDist.j <- ssn$obs$DistanceUpstream[
    match(ob.j$pid, ssn$obs$ng.pid)
  ]

  ob.j
}
