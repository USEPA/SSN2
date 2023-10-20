#' Get a data.frame from an SSN, ssn_lm, or ssn_glm object
#'
#' @description The \command{ssn_get_data} function extracts an
#'   sf data.frame for the observation or prediction data from
#'   an \code{SSN}, \code{ssn_lm}, or \code{ssn_glm} object.
#'
#' @param x An object of class \code{SSN}, \code{ssn_lm}, or \code{ssn_glm}.
#' @param name the internal name of the dataset in the object
#'   \code{x}. For observed values, this will always be "obs", the
#'   default.
#'
#' @details The internal \code{name} for observed data in objects of
#'   class \code{SSN} is "obs" and it is the
#'   default. If another \code{name} is specified, it must represent a
#'   prediction data set in the \code{SSN},
#'   \code{ssn_lm}, or \code{ssn_glm} object. For \code{SSN} objects,
#'   these names are obtained using the call \code{names(x$preds)}. For
#'   all other object classes, the names are obtained using the call
#'   \code{names(x$ssn.object$preds)}.
#'
#' @return An sf data.frame
#'
#' @name ssn_get_data
#' @export
#'
#' @seealso [ssn_put_data()]
#'
#' @examples
#' ## Extract observed data from an SSN object
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, predpts = "pred1km", overwrite = TRUE)
#'
#' obs.df <- ssn_get_data(mf04p)
#' dim(obs.df)
#'
#' ## Extract prediction data from an SSN object
#' names(mf04p$preds)
#' pred1km.df <- ssn_get_data(mf04p, name = "pred1km")
#' names(pred1km.df)
#'
#' ## extract observed data from an ssn_lm object
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' obs.mod.df <- ssn_get_data(ssn_mod)
#' summary(obs.mod.df)
ssn_get_data <- function(x, name = "obs") {
  ## Check input object type
  if (!class(x) %in% c("SSN", "ssn_lm", "ssn_glm")) {
    stop(paste0("Error: an object of class ", class(x), " is not a valid input."))
  }

  if (class(x) %in% c("ssn_lm", "ssn_glm")) x <- x$ssn.object


  if (name == "obs") {
    return(x$obs)
  }

  if (name != "obs") {
    np <- length(x$preds)
    if (length(np) == 0) stop("No prediction data sets in object")
    if (!name %in% attributes(x$preds)$names) {
      tmp <- names(x$preds)[1]
      for (k in 2:np) {
        tmp <- paste0(tmp, ", ", names(x$preds)[k])
      }
      stop(paste("name does not match any prediction data sets, try one of these: ", tmp))
    }
    return(x$preds[[name]])
  }
}
