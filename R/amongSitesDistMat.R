amongSitesDistMat <- function(ssn, pids, name = "obs", bin.table) {
  ## ssn = SSN object
  ## pids = list of pid values for prediction sites
  ## bin.table = binaryID table for the network

  site.no <- length(pids)

  among_distance_matrix <- matrix(NA, nrow = site.no, ncol = site.no)
  diag(among_distance_matrix) <- 0
  rownames(among_distance_matrix) <- pids
  colnames(among_distance_matrix) <- pids

  if (name != "obs") {
    ind.pids <- ssn$preds[[name]]$ng.pid %in% as.character(pids)
    locID.pid.data <- ssn$preds[[name]]$locID[ind.pids]
    pid.data <- ssn_get_netgeometry(ssn$preds[[name]][ind.pids, ], c(
      "pid", "SegmentID", "locID",
      "DistanceUpstream"
    ))
    pid.data <- as.data.frame(lapply(pid.data, as.numeric))
    colnames(pid.data) <- c("pid", "rid", "locID", "upDist")
  } else {
    ind.pids <- ssn$obs$ng.pid %in% as.character(pids)
    locID.pid.data <- ssn$obs$locID[ind.pids]
    pid.data <- ssn_get_netgeometry(ssn$obs[ind.pids, ], c(
      "pid", "SegmentID", "locID",
      "DistanceUpstream"
      ), reformat = TRUE)

    ##pid.data <- as.data.frame(sapply(pid.data, as.numeric))
    colnames(pid.data) <- c("pid", "rid", "locID", "upDist")
  }

  pid.data <- pid.data[order(pid.data$pid), ]

  ## Need bin.table
  pid.data$binaryID <- bin.table$binaryID[match(pid.data$rid, bin.table$rid)]
  pid.data <- pid.data[order(pid.data[, "pid"]), ]
  rownames(pid.data) <- pid.data$pid

  ## locID values can be repeated, in which case they have the same distance data.
  locID.old <- -1
  for (b in seq_len(site.no)) {
    locID.b <- pid.data[b, "locID"]
    upDist.b <- pid.data[b, "upDist"]
    pid.b <- pid.data[b, "pid"]

    if (locID.b != locID.old) {
      junk <- get.rid.fc(pid.data[, "binaryID"], pid.data$binaryID[b])
      truncated.binaryIDs <- data.frame(pid = pid.data[, "pid"], junk, stringsAsFactors = FALSE)
      truncated.binaryIDs$fc <- as.logical(truncated.binaryIDs$fc)
      truncated.binaryIDs$junc.rid <- bin.table$rid[match(truncated.binaryIDs$binaryID, bin.table$binaryID)]

      truncated.binaryIDs$juncDist <- ssn$edges$DistanceUpstream[match(
        truncated.binaryIDs$junc.rid,
        ssn$edges$rid
      )]
      truncated.binaryIDs$upDist.j <- pid.data$upDist[match(truncated.binaryIDs$pid, pid.data$pid)]
      ind.fc <- truncated.binaryIDs$fc == 1
      dist.sites <- ifelse(ind.fc, upDist.b - truncated.binaryIDs$upDist.j,
        upDist.b - truncated.binaryIDs$juncDist
      )
      among_distance_matrix[, paste(pid.b)] <- ifelse(dist.sites < 0, 0, dist.sites)
    }
  }
  among_distance_matrix
}
