#' @title Extract netgeometry column
#'
#' @description Extract topological information from netgeometry column
#'
#' @param x An sf data.frame found in an \code{SSN} object or the
#'   netgeometry column as a vector
#'
#' @param netvars Network coordinate variables to return. Default is
#'   "all". For edges, valid column names include: "NetworkID",
#'   "SegmentID", and "DistanceUpstream". For point datasets, valid column
#'   names include "NetworkID", "SegmentID", "DistanceUpstream", "ratio", "pid",
#'   and "locID".
#' @param reformat Convert network coordinate variables from character to numeric.
#'
#' @details When an \code{SSN} object is generated using the
#'   \code{importSSN} function, a text column named "netgeometry" is added
#'   to the edges, observed sites, and prediction sites (if they
#'   exist) data.frames. The netgeometry column contains data used to
#'   describe how edge and site features relate to one another in
#'   topological space. For edges, netgeometry values contain the
#'   "ENETWORK" prefix, with 3 space delimited values in parentheses:
#'   "ENETWORK (NetworkID SegmentID DistanceUpstream)". For point
#'   datasets (observed and prediction sites), the values contain the
#'   "SNETWORK" prefix, followed by 6 space delimited values in parentheses:
#'   "SNETWORK (NetworkID SegmentID DistanceUpstream ratio pid locID)". The
#'   \code{ssn_get_netgeometry} function extracts and converts these
#'   values from text to numeric, returning either a data.frame
#'   (default) or vector containing the variables requested via
#'   \code{netvars}.
#'
#' @return If more than one column is requested using netvars, the
#'   function returns a data.frame (default). If only one column is
#'   requested, the result is a vector.
#'
#' @name ssn_get_netgeometry
#' @export
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' ssn_get_netgeometry(mf04p$obs)
#' ssn_get_netgeometry(mf04p$edges, "DistanceUpstream")
ssn_get_netgeometry <- function(x, netvars = "all", reformat = FALSE) {
  # I think this should be an SSN obejct and we should have another column
  # for "type" which can be "edges", "obs", or a prediction name
  if (inherits(x, "SSN")) {
    stop("An object of class SSN is not a valid input", call. = FALSE)
  }


  if (inherits(x, "data.frame")) {
    x <- x$netgeometry
  }

  ## delete "network"
  x <- gsub("[ES]{1}NETWORK ", "", x) # include the space to avoid issues later
  ## delete brackets
  x <- gsub("[()]{1}", "", x)
  ## split string at spaces
  x_split <- strsplit(x, " ")
  ## Get a matrix
  x_df <- as.data.frame(do.call(rbind, x_split))
  ## Detect number of columns
  ncols <- NCOL(x_df)
  ## Then attach names
  if (ncols == 3) {
    colnames(x_df) <- c("NetworkID", "SegmentID", "DistanceUpstream")
  } else if (ncols == 6) {
    colnames(x_df) <- c("NetworkID", "SegmentID", "DistanceUpstream", "ratio", "pid", "locID")
  } else {
    stop("Invalid number of columns", call. = FALSE)
  }

  if (reformat == TRUE) {
    x_df <- as.data.frame(apply(x_df, 2, as.numeric, simplify = FALSE))
  }

  ## Return result
  if (any(netvars == "all")) {
    return(x_df)
  } else {
    if (sum(netvars %in% colnames(x_df)) != length(netvars)) {
      stop("Invalid netvars requested. Check spelling and data type.", call. = FALSE)
    }
    return(x_df[, c(netvars), drop = FALSE])
  }
}
