#' @rdname predict.SSN2
#' @method predict ssn_lmRF
#' @export
predict.ssn_lmRF <- function(object, newdata, local, se.fit = FALSE,
                              interval = "none", type = "response", block = FALSE, ...) {
  if (missing(local)) local <- NULL
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before predicting from an ssn_lmRF object.", call. = FALSE)
  }
  if (!identical(se.fit, FALSE)) {
    stop("se.fit is not supported for ssn_lmRF point predictions.", call. = FALSE)
  }
  if (!identical(interval, "none")) {
    stop("interval must be \"none\" for ssn_lmRF point predictions.", call. = FALSE)
  }
  if (!identical(type, "response")) {
    stop("type must be \"response\" for ssn_lmRF point predictions.", call. = FALSE)
  }
  if (!identical(block, FALSE)) {
    stop("block prediction is not supported for ssn_lmRF.", call. = FALSE)
  }

  ranger_args <- get_ssn_lmRF_predict_args(list(...))
  if (missing(newdata)) {
    newdata <- if (is.null(object$missing_prediction_name)) "all" else object$missing_prediction_name
  }
  newdata_names <- resolve_newdata_name(object$ssn_lm, newdata)
  available <- names(object$ssn_lm$ssn.object$preds)
  if (!length(newdata_names) || any(!newdata_names %in% available)) {
    stop("newdata must name one or more prediction sets in the fitted SSN object.", call. = FALSE)
  }

  predict_one <- function(newdata_name) {
    prediction_data <- get_prediction_newdata(object$ssn_lm, newdata_name)$newdata
    if (!NROW(prediction_data)) return(NULL)
    ranger_prediction <- do.call(
      stats::predict,
      c(
        list(object = object$ranger, data = as.data.frame(sf::st_drop_geometry(prediction_data)), type = "response"),
        ranger_args
      )
    )$predictions
    residual_prediction <- predict(object$ssn_lm, newdata = newdata_name, local = local)
    if (length(ranger_prediction) != length(residual_prediction) || any(!is.finite(ranger_prediction))) {
      stop("ranger prediction did not return one finite value per SSN prediction row.", call. = FALSE)
    }
    prediction <- as.numeric(ranger_prediction) + as.numeric(residual_prediction)
    names(prediction) <- names(residual_prediction)
    prediction
  }

  if (length(newdata_names) == 1L) return(predict_one(newdata_names))
  predictions <- lapply(newdata_names, predict_one)
  names(predictions) <- newdata_names
  predictions
}

#' Validate and pass through \code{ranger::predict.ranger()} arguments for \code{predict.ssn_lmRF()}
#'
#' @param dots A named list of \code{...} arguments from
#'   \code{predict.ssn_lmRF()}.
#'
#' @return \code{dots}, if every argument is named and one of the supported
#'   \code{ranger::predict.ranger()} controls (\code{num.trees}, \code{seed},
#'   \code{num.threads}, \code{verbose}); otherwise an error.
#'
#' @noRd
get_ssn_lmRF_predict_args <- function(dots) {
  dot_names <- names(dots)
  if (length(dots) && (is.null(dot_names) || any(!nzchar(dot_names)))) {
    stop("All ssn_lmRF prediction arguments must be named.", call. = FALSE)
  }
  unsupported <- intersect(dot_names, c("scale", "df", "level", "terms", "na.action"))
  if (length(unsupported)) {
    stop(
      "Unsupported ssn_lmRF prediction argument(s): ", paste(unsupported, collapse = ", "), ".",
      call. = FALSE
    )
  }
  allowed <- c("num.trees", "seed", "num.threads", "verbose")
  unknown <- setdiff(dot_names, allowed)
  if (length(unknown)) {
    stop(
      "Unsupported ssn_lmRF prediction argument(s): ", paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  dots
}
