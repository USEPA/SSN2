#' Fit a random forest stream-network residual model
#'
#' @description Fit random forest spatial residual models for stream network data
#'   using random forest to fit the mean and a spatial linear model to fit the residuals.
#'   The spatial linear model fit to the residuals can incorporate variety of estimation
#'   methods, allowing for random effects, anisotropy, partition factors, and big data methods.
#'
#' @param formula A two-sided linear formula describing the fixed effect structure
#'   of the model, with the response to the left of the \code{~} operator and
#'   the terms on the right, separated by \code{+} operators.
#' @param ssn.object A spatial stream network object with class \code{SSN}.
#' @param ... Additional named arguments to [ranger::ranger()] or [ssn_lm()].
#'
#' @details The random forest residual spatial linear model is described by
#'   Fox et al. (2020). A random forest model is fit to the mean portion of the
#'   model specified by \code{formula} using \code{ranger::ranger()}. Residuals
#'   are computed and used as the response variable in an intercept-only spatial
#'   linear model fit using [ssn_lm()]. This model object is intended for use with
#'   \code{predict()} to perform prediction, also called random forest
#'   regression Kriging.
#' 
#' @return A list with several elements to be used with \code{predict()}. These
#'   elements include the function call (named \code{call}), the random forest object
#'   fit to the mean (named \code{ranger}),
#'   the spatial stream network model object fit to the residuals
#'   (named \code{ssn_lm}), and an object can contain data for
#'   locations at which to predict (called \code{newdata}). The \code{newdata}
#'   object contains the set of
#'   observations in \code{data} whose response variable is \code{NA}.
#'
#' @note This function does not perform any internal scaling. If optimization is not
#'   stable due to extremely large variances, scale relevant variables
#'   so they have variance 1 before optimization.
#'
#' @references
#' Fox, E.W., Ver Hoef, J.M., & Olsen, A.R. (2020). Comparing spatial
#' regression to random forests for large environmental data sets.
#' \emph{PLoS ONE}, 15(3), e0229509.
#'
#' @export
#'
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf <- ssn_import(temp_path,
#'   predpts = "CapeHorn", overwrite = TRUE
#' )
#' ssn_create_distmat(mf, predpts = "CapeHorn", overwrite = TRUE)
#' rf_fit <- ssn_lmRF(
#'   Summer_mn ~ ELEV_DEM, mf, tailup_type = "exponential",
#'   additive = "afvArea", num.trees = 50, num.threads = 1, seed = 1
#' )
#' print(rf_fit)
#' head(predict(rf_fit, "CapeHorn"))
#' }
ssn_lmRF <- function(formula, ssn.object, ...) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install the ranger package before using ssn_lmRF().", call. = FALSE)
  }
  if (length(attr(stats::terms(formula), "offset"))) {
    stop("Offsets are not supported by ssn_lmRF(); ranger has no compatible offset interface.", call. = FALSE)
  }

  dot_expressions <- as.list(match.call(expand.dots = FALSE)$...)
  if (is.null(dot_expressions)) dot_expressions <- list()
  if ("additive" %in% names(dot_expressions) && is.symbol(dot_expressions$additive)) {
    dot_expressions$additive <- deparse1(dot_expressions$additive)
  }
  dots <- lapply(dot_expressions, eval, envir = parent.frame())
  rf_args <- get_ssn_lmRF_fit_args(dots)
  rf_args$ranger <- check_ssn_lmRF_ranger_args(rf_args$ranger)

  obdata <- ssn.object$obs
  forest_data <- sf::st_drop_geometry(obdata)
  model_frame <- stats::model.frame(formula, forest_data, na.action = stats::na.pass)
  response <- stats::model.response(model_frame)
  if (!is.numeric(response) || is.matrix(response)) {
    stop("The ssn_lmRF response must be a single numeric variable.", call. = FALSE)
  }
  missing_response <- is.na(response)
  if (all(missing_response)) {
    stop("The ssn_lmRF response cannot be entirely missing.", call. = FALSE)
  }
  if (any(!is.finite(response[!missing_response]))) {
    stop("The ssn_lmRF response must be finite where observed.", call. = FALSE)
  }

  ranger_fit <- do.call(
    ranger::ranger,
    c(list(formula = formula, data = forest_data[!missing_response, , drop = FALSE]), rf_args$ranger)
  )
  check_ssn_lmRF_ranger_fit(ranger_fit, sum(!missing_response))

  residual_name <- make.unique(c(names(obdata), ".ssn_lmRF_residual"))[length(names(obdata)) + 1L]
  residual <- rep(NA_real_, NROW(obdata))
  residual[!missing_response] <- response[!missing_response] - ranger_fit$predictions
  residual_ssn <- ssn.object
  residual_ssn$obs[[residual_name]] <- residual
  residual_formula <- stats::as.formula(paste(residual_name, "~ 1"), env = environment(formula))
  residual_fit <- do.call(
    ssn_lm,
    c(list(formula = residual_formula, ssn.object = residual_ssn), rf_args$ssn_lm)
  )

  ranger_fit$call <- NA
  residual_fit$call <- NA
  structure(
    list(
      call = match.call(),
      ranger = ranger_fit,
      ssn_lm = residual_fit,
      missing_prediction_name = if (any(missing_response)) ".missing" else NULL,
      residual_name = residual_name
    ),
    class = "ssn_lmRF"
  )
}

#' Split \code{ssn_lmRF()}'s \code{...} into \code{ranger::ranger()} and \code{ssn_lm()} argument groups
#'
#' @param dots A named list of \code{...} arguments from \code{ssn_lmRF()}.
#'
#' @return A list with \code{ranger} and \code{ssn_lm}, each a named list of
#'   the arguments belonging to that function; errors if any argument is
#'   unnamed or unrecognized by either function (or the optimizer/model
#'   arguments \code{method}/\code{control}/\code{hessian}/\code{contrasts}).
#'
#' @noRd
get_ssn_lmRF_fit_args <- function(dots) {
  dot_names <- names(dots)
  if (length(dots) && (is.null(dot_names) || any(!nzchar(dot_names)))) {
    stop("All arguments supplied to ssn_lmRF() must be named.", call. = FALSE)
  }
  ranger_names <- setdiff(names(formals(ranger::ranger)), c("formula", "data", "..."))
  ssn_lm_names <- setdiff(names(formals(ssn_lm)), c("formula", "ssn.object", "..."))
  optim_names <- c("method", "control", "hessian")
  model_names <- "contrasts"
  allowed <- union(ranger_names, union(ssn_lm_names, c(optim_names, model_names)))
  unknown <- setdiff(dot_names, allowed)
  if (length(unknown)) {
    stop(
      "Unsupported ssn_lmRF argument(s): ", paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  ranger_index <- dot_names %in% ranger_names
  ssn_lm_index <- dot_names %in% c(ssn_lm_names, optim_names, model_names)
  list(ranger = dots[ranger_index], ssn_lm = dots[ssn_lm_index])
}

#' Validate and default \code{ssn_lmRF()}'s forwarded \code{ranger::ranger()} arguments
#'
#' Rejects an explicit \code{write.forest = FALSE}, \code{oob.error = FALSE},
#' or \code{probability = TRUE} (each incompatible with random forest
#' regression kriging), and defaults \code{write.forest}/\code{oob.error} to
#' \code{TRUE} when not supplied.
#'
#' @param ranger_args A named list of arguments to forward to
#'   \code{ranger::ranger()}.
#'
#' @return \code{ranger_args}, with \code{write.forest}/\code{oob.error}
#'   defaulted to \code{TRUE}.
#'
#' @noRd
check_ssn_lmRF_ranger_args <- function(ranger_args) {
  if ("write.forest" %in% names(ranger_args) && !isTRUE(ranger_args$write.forest)) {
    stop("ssn_lmRF() requires write.forest = TRUE for prediction.", call. = FALSE)
  }
  if ("oob.error" %in% names(ranger_args) && !isTRUE(ranger_args$oob.error)) {
    stop("ssn_lmRF() requires oob.error = TRUE to form residuals.", call. = FALSE)
  }
  if ("probability" %in% names(ranger_args) && isTRUE(ranger_args$probability)) {
    stop("ssn_lmRF() requires a regression forest (probability = FALSE).", call. = FALSE)
  }
  if (!("write.forest" %in% names(ranger_args))) ranger_args$write.forest <- TRUE
  if (!("oob.error" %in% names(ranger_args))) ranger_args$oob.error <- TRUE
  ranger_args
}

#' Validate a fitted \code{ranger} forest is usable for \code{ssn_lmRF()}
#'
#' @param ranger_fit A fitted \code{ranger::ranger()} model object.
#' @param n The expected number of (non-missing-response) out-of-bag
#'   predictions.
#'
#' @return \code{ranger_fit}, invisibly, if it is a regression forest with a
#'   saved forest and \code{n} finite out-of-bag predictions; otherwise an
#'   error.
#'
#' @noRd
check_ssn_lmRF_ranger_fit <- function(ranger_fit, n) {
  if (!identical(ranger_fit$treetype, "Regression")) {
    stop("ssn_lmRF() requires ranger to fit a regression forest.", call. = FALSE)
  }
  if (is.null(ranger_fit$forest)) {
    stop("ssn_lmRF() requires a saved ranger forest.", call. = FALSE)
  }
  if (length(ranger_fit$predictions) != n || any(!is.finite(ranger_fit$predictions))) {
    stop("ssn_lmRF() requires finite ranger out-of-bag predictions for every observed response.", call. = FALSE)
  }
  invisible(ranger_fit)
}
