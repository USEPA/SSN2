#' @method AUROC ssn_glm
#' @export
AUROC.ssn_glm <- function(object, ...) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Install the pROC package before using AUROC().", call. = FALSE)
  } else {

    if (object$family != "binomial") {
      stop("AUROC() only available when family is \"binomial\".", call. = FALSE)
    }

    # AUROC requires binary (single trial) outcomes; a binomial model with
    # size > 1 represents aggregated successes/failures, not a single
    # classification label, so ROC curves don't apply
    if (any(object$size != 1)) {
      stop("AUROC() only available for binary models (i.e., models whose response indicates a single success or failure).", call. = FALSE)
    }
    dotlist <- list(...)
    if (!("quiet" %in% names(dotlist))) {
      # suppress pROC's default console message unless the user asked for it
      dotlist$quiet <- TRUE
    }
    # compare observed binary outcomes to fitted probabilities across all
    # classification thresholds to get the area under the ROC curve
    as.numeric(do.call(pROC::auc, c(list(response = as.vector(object$y), predictor = fitted(object)), dotlist)))
  }
}
