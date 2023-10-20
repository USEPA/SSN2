#' Glance at many fitted model objects
#'
#' @description \code{glances()} repeatedly calls \code{glance()} on several
#'   fitted model objects and binds the output together, sorted by a column of interest.
#'
#' @param object Fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param ... Additional fitted model objects from [ssn_lm()] or [ssn_glm()].
#' @param sort_by Sort by a \code{glance} statistic (i.e., the name of a column
#'   output from \code{glance()} or the order of model input (\code{sort_by = "order"}).
#'   The default is \code{"AICc"}.
#' @param decreasing Should \code{sort_by} be decreasing or not? The default is \code{FALSE}.
#'
#' @return A tibble where each row represents the output of \code{glance()} for
#'   each fitted model object.
#'
#' @name glances.SSN2
#' @method glances ssn_lm
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
#' # tailup only
#' ssn_mod1 <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' # taildown only
#' ssn_mod2 <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   taildown_type = "exponential"
#' )
#' glances(ssn_mod1, ssn_mod2)
glances.ssn_lm <- function(object, ..., sort_by = "AICc", decreasing = FALSE) {
  model_list <- c(list(object), list(...))
  model_list_names <- c(as.character(as.list(substitute(list(object)))[-1]), as.character(as.list(substitute(list(...)))[-1]))
  model_glance <- lapply(model_list, function(x) glance(x))
  model_bind <- do.call(rbind, model_glance)
  model_bind <- cbind(data.frame(model = model_list_names), model_bind)
  if (sort_by == "order") {
    model_bind <- model_bind[order(seq_len(NROW(model_bind)), decreasing = decreasing), , drop = FALSE]
  } else {
    model_bind <- model_bind[order(model_bind[[substitute(sort_by)]], decreasing = decreasing), , drop = FALSE]
  }
  tibble::as_tibble(model_bind)
}

#' @rdname glances.SSN2
#' @method glances ssn_glm
#' @export
glances.ssn_glm <- glances.ssn_lm
