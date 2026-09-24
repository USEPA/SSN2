#' Create a covariance matrix
#'
#' Create a covariance matrix from a fitted model object.
#'
#' @param object A fitted model object (e.g., [ssn_lm()] or [ssn_glm()]).
#' @param newdata If omitted, the covariance matrix of
#'   the observed data is returned. If provided, \code{newdata} is
#'   a character string naming the prediction data set (accessible via
#'   \code{object$ssn.object$preds}) for which the covariance is desired.
#'   Note that the prediction data must be in the original SSN object used
#'   to fit \code{object}.
#' @param cov_type The type of covariance matrix returned. If \code{newdata}
#'   is omitted, the \eqn{n \times n} covariance matrix of the observed
#'   data is returned, where \eqn{n} is the sample size used to fit \code{object}.
#'   If \code{newdata} is provided and \code{cov_type} is \code{"pred.obs"} (the default),
#'   the \eqn{m \times n} covariance matrix of the predicted and observed data is returned,
#'   where \eqn{m} is the number of observations in the prediction data.
#'   If \code{newdata} is provided and \code{cov_type} is \code{"obs.pred"},
#'   the \eqn{n \times m} covariance matrix of the observed and prediction data is returned.
#'   If \code{newdata} is provided and \code{cov_type} is \code{"pred.pred"},
#'   the \eqn{m \times m} covariance matrix of the prediction data is returned.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return If \code{newdata} is omitted, the covariance matrix of the observed
#'   data, which has dimension n x n, where n is the sample size used to fit \code{object}.
#'   If \code{newdata} is provided, the covariance matrix between the unobserved (new)
#'   data and the observed data, which has dimension m x n, where m is the number of
#'   new observations and n is the sample size used to fit \code{object}.
#'
#' @name covmatrix.SSN2
#' @method covmatrix ssn_lm
#' @order 1
#' @export
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
#'
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' covmatrix(ssn_mod)
#' covmatrix(ssn_mod, "CapeHorn")
covmatrix.ssn_lm <- function(object, newdata, cov_type, ...) {
  params_object <- object$coefficients$params_object

  initial_object <- get_initial_object_from_coef(object)

  if (missing(newdata)) {
    cov_type <- "obs.obs"
  } else if (missing(cov_type)) {
    cov_type <- "pred.obs"
  }

  if (cov_type != "obs.obs" && is.null(newdata)) {
    stop("No prediction data for which to create a covariance matrix.", call. = FALSE)
  }

  if (cov_type == "obs.obs") {
    de_scale <- sum(params_object$tailup[["de"]], params_object$taildown[["de"]], params_object$euclid[["de"]])
    randcov_names <- get_randcov_names(object$random)
    randcov_Zs <- get_randcov_Zs(object$ssn.object$obs, randcov_names)
    partition_matrix_val <- partition_matrix(object$partition_factor, object$ssn.object$obs)
    if ((is.logical(object$missing_index) && sum(object$missing_index) > 0) || (!is.logical(object$missing_index) && length(object$missing_index) > 0)) { # this (and list format below) is for "putting stuff back together" when there is missingness in observed data
      object$ssn.object$obs <- rbind(object$ssn.object$obs, object$ssn.object$preds$.missing)
      reorder_val <- order(c(object$observed_index, object$missing_index))
      object$ssn.object$obs <- object$ssn.object$obs[reorder_val, , drop = FALSE]
    }
    tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
    taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")
    backend <- select_square_dist_backend(object$ssn.object, "obs", tailup_none, taildown_none)
    dist_object <- get_dist_object(object$ssn.object, initial_object, object$additive, object$anisotropy, backend = backend)
    # this is to subset the data by observed index
    dist_object <- get_dist_object_oblist(dist_object, object$observed_index, local_index = rep(1, object$n))
    dist_object <- dist_object[[1]] # unlist
    cov_val <- get_cov_matrix(params_object, dist_object, randcov_Zs, partition_matrix_val,
      object$anisotropy,
      de_scale = de_scale, diagtol = object$diagtol
    )
  } else if (cov_type == "obs.pred" || cov_type == "pred.obs") {
    newdata_name <- newdata
    newdata <- object$ssn.object$preds[[newdata_name]]
    dist_pred_object <- get_dist_pred_object(object, newdata_name, initial_object)
    cov_val <- get_cov_vector(params_object, dist_pred_object, object$ssn.object$obs, newdata, object$partition_factor, object$anisotropy)
    if (cov_type == "obs.pred") {
      cov_val <- t(cov_val)
    }
  } else if (cov_type == "pred.pred") {
    newdata_name <- newdata
    newdata <- object$ssn.object$preds[[newdata_name]]
    if (identical(newdata_name, ".missing")) {
      rows <- seq_len(NROW(newdata))
      return(get_block_pred_covariance(object, newdata_name, rows, rows))
    }
    de_scale <- sum(params_object$tailup[["de"]], params_object$taildown[["de"]], params_object$euclid[["de"]])
    randcov_names <- get_randcov_names(object$random)
    extended_random_xlev <- extend_randcov_xlev(object$random_xlev, newdata, randcov_names)
    randcov_Zs <- get_randcov_Zs(newdata, randcov_names, xlev_list = extended_random_xlev)
    partition_matrix_val <- partition_matrix(object$partition_factor, newdata)
    tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
    taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")
    backend <- select_square_dist_backend(object$ssn.object, newdata_name, tailup_none, taildown_none)
    dist_predbk_object <- get_dist_predbk_object(object, newdata_name, initial_object, backend = backend)
    cov_val <- get_cov_matrix(params_object, dist_predbk_object, randcov_Zs, partition_matrix_val,
      object$anisotropy,
      de_scale = de_scale, diagtol = object$diagtol
    )
  } else {
    stop('cov_type must be "obs.obs", "obs.pred", "pred.obs", "pred.pred"', call. = FALSE)
  }

  # return covariance value as a base R matrix (not a Matrix matrix)
  as.matrix(cov_val)
}

#' @rdname covmatrix.SSN2
#' @method covmatrix ssn_glm
#' @order 2
#' @export
covmatrix.ssn_glm <- covmatrix.ssn_lm
