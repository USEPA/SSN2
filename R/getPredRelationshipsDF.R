#' Get spatial relationships for prediction points
#'
#' @param ssn SSN object
#' @param predpts Index of predpts
#' @param ind Index for pred values on the network
#' @param bin binaryID table for the network
#' @param ob Data frame  with pid, rid, locID, and binaryID for sites on the
#'   network ordered by pid (ob.i)
#' @param j Row j of data.frame ob
#'
#' @noRd
getPredRelationshipsDF <- function(ssn, predpts, ind, bin, ob, j) {
  ## ssn = SpatialStreamNetwork object
  ## num = index of predpts in SSN
  ## ind = index for pred values that lie on this network
  ## bin = binaryID table for this network
  ## ob =  data.frame with pid, rid, locID, and binaryID for sites on network
  ##         ordered by pid (ob.i)
  ## j = row j of data.frame ob

  pred.tmp <- as.data.frame(cbind(
    ssn$preds[[predpts]]$ng.pid[ind],
    ssn$preds[[predpts]]$SegmentID[ind]
  ))
  colnames(pred.tmp) <- c("pid", "rid")

  pred.tmp$binaryID <- bin$binaryID[match(pred.tmp$rid, bin$rid)]
  pred.tmp <- pred.tmp[order(pred.tmp[, "pid"]), ]
  rownames(pred.tmp) <- pred.tmp$pid

  junk <- get.rid.fc(pred.tmp[, "binaryID"], ob$binaryID[j])
  ob.j <- data.frame(pred.tmp["pid"], junk, stringsAsFactors = FALSE)

  ob.j$pid <- as.numeric(ob.j$pid)
  ob.j$fc <- as.logical(ob.j$fc)

  ob.j$junc.rid <- bin$rid[match(ob.j$binaryID, bin$binaryID)]
  ob.j$juncDist <- ssn$edges$DistanceUpstream[match(ob.j$junc.rid, ssn$edges$SegmentID)]
  ob.j$upDist.j <- ssn$preds[[predpts]]$DistanceUpstream[match(
    ob.j$pid,
    as.numeric(ssn$preds[[predpts]]$ng.pid)
  )]
  ob.j
}
