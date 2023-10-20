#' Create a covariance matrix
#'
#' Create a covariance matrix from a fitted model object.
#'
#' @param object A fitted model object (e.g., [ssn_lm()] or [ssn_glm()]).
#' @param newdata If omitted, the covariance matrix of
#'   the observed data is returned. If provided, \code{newdata} is
#'   a data frame or \code{sf} object that contains coordinate information
#'   required to construct the covariance between \code{newdata} and
#'   the observed data. If a data frame, \code{newdata}
#'   must contain variables that represent coordinates having the same name as
#'   the coordinates from the observed data used to fit \code{object}. If an
#'   \code{sf} object, coordinates are obtained from the geometry of \code{newdata}.
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
#' @return A covariance matrix (see \code{cov_type}).
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


  tailup_type <- remove_covtype(class(coef(object, type = "tailup")))
  taildown_type <- remove_covtype(class(coef(object, type = "taildown")))
  euclid_type <- remove_covtype(class(coef(object, type = "euclid")))
  nugget_type <- remove_covtype(class(coef(object, type = "nugget")))

  initial_object <- get_initial_object(
    tailup_type = tailup_type,
    taildown_type = taildown_type,
    euclid_type = euclid_type,
    nugget_type = nugget_type,
    tailup_initial = NULL,
    taildown_initial = NULL,
    euclid_initial = NULL,
    nugget_initial = NULL
  )

  if (missing(newdata)) {
    cov_type <- "obs.obs"
  } else if (missing(cov_type)) {
    cov_type <- "pred.obs"
  }


  if (cov_type == "obs.obs") {
    de_scale <- sum(params_object$tailup[["de"]], params_object$taildown[["de"]], params_object$euclid[["de"]])
    randcov_names <- get_randcov_names(object$random)
    randcov_Zs <- get_randcov_Zs(object$ssn.object$obs, randcov_names)
    partition_matrix_val <- partition_matrix(object$partition_factor, object$ssn.object$obs)
    if (length(object$missing_index) > 0) { # this (and list format below) is for "putting stuff back together" when there is missingness in observed data
      object$ssn.object$obs <- rbind(object$ssn.object$obs, object$ssn.object$preds$.missing)
      reorder_val <- order(c(object$observed_index, object$missing_index))
      object$ssn.object$obs <- object$ssn.object$obs[reorder_val, , drop = FALSE]
    }
    dist_object <- get_dist_object(object$ssn.object, initial_object, object$additive, object$anisotropy)
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
    de_scale <- sum(params_object$tailup[["de"]], params_object$taildown[["de"]], params_object$euclid[["de"]])
    randcov_names <- get_randcov_names(object$random)
    randcov_Zs <- get_randcov_Zs(newdata, randcov_names)
    partition_matrix_val <- partition_matrix(object$partition_factor, newdata)
    dist_predbk_object <- get_dist_predbk_object(object, newdata_name, initial_object)
    cov_val <- get_cov_matrix(params_object, dist_predbk_object, randcov_Zs, partition_matrix_val,
      object$anisotropy,
      de_scale = de_scale, diagtol = object$diagtol
    )
  } else {
    stop("Invalid \"cov_type\" argument.", call. = FALSE)
  }

  # return covariance value as a base R matrix (not a Matrix matrix)
  as.matrix(cov_val)
}

#' @rdname covmatrix.SSN2
#' @method covmatrix ssn_glm
#' @order 2
#' @export
covmatrix.ssn_glm <- covmatrix.ssn_lm
