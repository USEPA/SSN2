#' Get initial object for use with Torgegram calculations
#'
#' @param type The Torgegram type (flow-connected, flow-unconnected, Euclidean)
#'
#' @noRd
get_Torgegram_initial_object <- function(type) {
  if (any(!type %in% c("flowcon", "flowuncon", "euclid"))) {
    stop("All elements of type must be \"flowcon\", \"flowuncon\", or \"euclid\".", call. = FALSE)
  }

  if ("flowcon" %in% type || "flowuncon" %in% type) {
    taildown_type <- "exponential" # set this so initial object specified
  } else {
    taildown_type <- "none"
  }

  if ("euclid" %in% type) {
    euclid_type <- "exponential" # set this so initial object specified
  } else {
    euclid_type <- "none"
  }

  initial_object <- get_initial_object(
    tailup_type = "none",
    taildown_type = taildown_type,
    euclid_type = euclid_type,
    nugget_type = "none",
    tailup_initial = NULL,
    taildown_initial = NULL,
    euclid_initial = NULL,
    nugget_initial = NULL
  )
}
