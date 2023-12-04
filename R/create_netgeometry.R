## Add netgeom column to sf data.frame sf_data
create_netgeom <- function(sf_data, type = NULL) {
  if (type == "point") {
    sf_data[, "netgeom"] <- paste0("SNETWORK (", paste(
      sf_data$netID, sf_data$rid, sf_data$upDist,
      sf_data$ratio, sf_data$pid, sf_data$locID
    ), ")", sep = "")
  } else {
    sf_data[, "netgeom"] <- paste0("ENETWORK (", paste(
      sf_data$netID,
      sf_data$rid,
      sf_data$upDist
    ),
    ")",
    sep = ""
    )
  }
  return(sf_data) ## Return sf data.frame with netgeom column added
}
