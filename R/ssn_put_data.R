#' Put an sf data.frame in an SSN object
#'
#' @description The \command{ssn_put_data} function puts an sf
#'   data.frame representing observation or prediction
#'   data into an SSN, ssn_lm, or ssn_glm object.
#'
#' @param data sf data.frame with point geometry.
#' @param x An object of class SSN, ssn_lm, or ssn_glm.
#' @param name the internal name of the data set in the object
#'   \code{x}. For observed data, this will always be "obs", the
#'   default.
#' @param resize_data Logical. Indicates whether sf_df can have a
#'   different number of features than the current data.frame in the
#'   object. Default is FALSE.
#'
#' @details The internal \code{name} for observed data in objects of
#'   class \code{SSN}, \code{ssn_lm}, and \code{ssn_glm} is "obs" and it is the
#'   default. If another \code{name} is specified, it must represent a
#'   prediction dataset in the object. For \code{SSN} objects,
#'   these names are obtained using the call \code{names(x$preds)}. For
#'   all other object classes, the names are obtained using the call
#'   names(x$ssn.object$preds).
#'
#'   The \code{resize_sf_data} argument specifies whether sf_data can have a
#'   different number of features (i.e., rows) than the sf data.frame
#'   it is replacing. Care should be taken when resize_df is set to
#'   TRUE, especially if the new sf_data has more features than the
#'   existing sf data.frame. In these cases, the user is responsible
#'   for ensuring that the additional features have the correct
#'   spatial, topological, and attribute data to accurately represent
#'   spatial relationships in the SSN object.
#'
#' @return Returns an object of the same class as x, which contains
#'   the sf data.frame sf_data.
#'
#' @name ssn_put_data
#' @export
#'
#' @seealso [ssn_get_data()]
#'
#' @examples
#' data(mf04p)
#' ## Extract observation data.frame from SSN object
#' obs.df <- ssn_get_data(mf04p)
#' ## Create a new column for summer mean temperature and set Value in
## first row to NA
#' obs.df$Value <- obs.df$Summer_mn
#' obs.df$Value[1] <- NA
#'
#' ## Put the modified sf data.frame into the SSN object
#' mf04p <- ssn_put_data(obs.df, mf04p)
#' head(ssn_get_data(mf04p)[, c("Summer_mn", "Value")])
ssn_put_data <-
  function(data, x, name = "obs", resize_data = FALSE) {
    # changing argument names
    sf_data <- data
    resize_sf_data <- resize_data

    ## Check input object types
    if (sum(class(sf_data)[1] != "sf" |
      st_geometry_type(sf_data) != "POINT") > 0) {
      stop(paste0("sf_data must be an sf data.frame with POINT geometry"))
    }

    if (!"netgeom" %in% names(sf_data)) {
      stop("'netgeom' column is missing in sf_data\n")
    }

    if (!inherits(x, c("SSN", "ssn_lm", "ssn_glm"))) {
      stop(paste0("Error: an object of class ", class(x), " is not a valid input."))
    }

    ## Set x to SSN object for ssn_lm and ssn_glm objects
    if (inherits(x, c("ssn_lm", "ssn_glm"))) {
      x1 <- x$ssn.object
    }
    if (inherits(x, "SSN")) {
      x1 <- x
    }

    ## Replace obs data.frame
    if (name == "obs") {
      ## Compare number of rows in sf_data vs x1$obs
      if ((nrow(sf_data) != nrow(x1$obs)) & resize_sf_data == FALSE) {
        stop("The current sf data.frame in x has a different number of point features than sf_data and resize_sf_data == FALSE. If this is not an error, set resize_sf_data = TRUE.")
      } else {
        x1$obs <- sf_data
      }
    }

    ## Replace preds data.frame
    if (name != "obs") {
      np <- length(x1$preds)
      if (length(np) == 0) stop("No prediction data sets in object")

      if (!(name %in% attributes(x1$preds)$names)) {
        stop(paste(
          "name does not match any prediction datasets, try one of these: ",
          names(x1$preds)
        ))
      }

      ## Compare number of rows in sf_data vs x1$preds[[name]]
      if ((nrow(sf_data) != nrow(x1$preds[[name]])) & resize_sf_data == FALSE) {
        stop("The current sf data.frame in x has a different number of point features than sf_data and resize_sf_data == FALSE. If this is not an error, set resize_sf_data = TRUE.")
      } else {
        x1$preds[[name]] <- sf_data
      }
    }
    if (inherits(x, c("ssn_lm", "ssn_glm"))) x$ssn.object <- x1
    if (inherits(x, "SSN")) x <- x1

    return(x)
  }
