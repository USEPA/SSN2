#' @method AUROC ssn_glm
#' @export
AUROC.ssn_glm <- function(object, ...) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Install the pROC package before using AUROC().", call. = FALSE)
  } else {

    if (object$family != "binomial") {
      stop("AUROC() only available when family is \"binomial\".", call. = FALSE)
    }

    if (any(object$size != 1)) {
      stop("AUROC() only available for binary models (i.e., models whose response indicates a single success or failure).", call. = FALSE)
    }
    dotlist <- list(...)
    if (!("quiet" %in% names(dotlist))) {
      dotlist$quiet <- TRUE
    }
    as.numeric(do.call(pROC::auc, c(list(response = as.vector(object$y), predictor = fitted(object)), dotlist)))
  }
}
