#' Augment data with information from fitted model objects
#'
#' @description Augment accepts a fitted model object and a data set and adds
#'   information about each observation in the data set. New columns always
#'   begin with a \code{.} prefix to avoid overwriting columns in the original
#'   data set.
#'
#'   Augment behaves differently depending on whether the original data or new data
#'   requires augmenting. Typically, when augmenting the original data, only the fitted
#'   model object is specified, and when augmenting new data, the fitted model object
#'   and \code{newdata} are specified. When augmenting the original data, diagnostic
#'   statistics are augmented to each row in the data set. When augmenting new data,
#'   predictions and optional intervals (confidence or prediction) or standard errors are augmented to each
#'   row in the new data set.
#'
#' @param x A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param drop A logical indicating whether to drop extra variables in the
#'   fitted model object \code{x} when augmenting. The default for \code{drop} is \code{TRUE}.
#'   \code{drop} is ignored if augmenting \code{newdata}.
#' @param newdata A vector that contains the names of the prediction \code{sf}
#'   objects from the original \code{ssn.object} requiring prediction.
#'   All of the original explanatory variables used to create the fitted model object \code{x}
#'   must be present in each prediction \code{sf} object represented by \code{newdata}.
#'   Defaults to \code{NULL}, which indicates
#'   that nothing has been passed to \code{newdata} and augmenting occurs
#'   for the original data. The value \code{"ssn"} is shorthand for specifying
#'   all prediction \code{sf} objects.
#' @param se_fit Logical indicating whether or not a \code{.se.fit} column should
#'   be added to augmented output. Passed to \code{predict()} and
#'   defaults to \code{FALSE}.
#' @param interval Character indicating the type of confidence interval columns to
#'   add to the augmented \code{newdata} output. Passed to \code{predict()} and defaults
#'   to \code{"none"}.
#' @param level Tolerance/confidence level. The default is \code{0.95}.
#' @param local A list or logical. If a list, specific list elements described
#'   in [predict.ssn_lm()] or [predict.ssn_glm()] control the big data approximation behavior.
#'   If a logical, \code{TRUE} chooses default list elements for the list version
#'   of \code{local} as specified in [predict.ssn_lm()] or [predict.ssn_glm()]. Defaults to \code{FALSE},
#'   which performs exact computations.
#' @param ... Additional arguments to \code{predict()} when augmenting \code{newdata}.
#'
#' @details \code{augment()} returns a tibble as an \code{sf} object.
#'
#'   Missing response values from the original data can be augmented as if
#'   they were a \code{newdata} object by providing \code{".missing"} to the
#'   \code{newdata} argument.
#'
#' @return When augmenting the original data set, a tibble with additional columns
#'   \itemize{
#'     \item \code{.fitted}: Fitted value
#'     \item \code{.resid}: Response residual (the difference between observed and fitted values)
#'     \item \code{.hat}: Leverage (diagonal of the hat matrix)
#'     \item \code{.cooksd}: Cook's distance
#'     \item \code{.std.resid}: Standardized residuals
#'     \item \code{.se.fit}: Standard error of the fitted value.
#'   }
#'
#'   When augmenting a new data set, a tibble with additional columns
#'   \itemize{
#'     \item \code{.fitted}: Predicted (or fitted) value
#'     \item \code{.lower}: Lower bound on interval
#'     \item \code{.upper}: Upper bound on interval
#'     \item \code{.se.fit}: Standard error of the predicted (or fitted) value
#'   }
#'
#'   When predictions for all prediction objects are desired, the output is a list
#'   where each element has a name that matches the prediction objects and values
#'   that are the predictions.
#'
#' @name augment.SSN2
#' @method augment ssn_lm
#' @order 1
#' @export
#'
#' @seealso [tidy.SSN2()] [glance.SSN2()]
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
#' augment(ssn_mod)
#' augment(ssn_mod, newdata = "CapeHorn")
augment.ssn_lm <- function(x, drop = TRUE, newdata = NULL, se_fit = FALSE,
                           interval = c("none", "confidence", "prediction"),
                           level = 0.95, local, ...) {



  interval <- match.arg(interval)

  # set data and newdata
  if (is.null(newdata)) {
    if (drop) {
      data <- model.frame(x)
    } else {
      data <- sf::st_drop_geometry(x$ssn.object$obs)
    }
  } else {
    data <- model.frame(x)
  }

  if (is.null(newdata)) {
    augment_data <- tibble::tibble(.fitted = fitted(x))
    if (se_fit) {
      preds_data <- predict(x, newdata = data, se.fit = se_fit, interval = "confidence", ...)
      augment_data$.se.fit <- preds_data$se.fit
    }
    tibble_out <- tibble::tibble(cbind(data, augment_data, influence(x)))
    tibble_out$pid <- ssn_get_netgeom(x$ssn.object$obs, netvars = "pid")$pid
    coords <- sf::st_coordinates(x$ssn.object$obs)
    tibble_out$.xcoord <- coords[, 1, drop = TRUE]
    tibble_out$.ycoord <- coords[, 2, drop = TRUE]
    tibble_out <- sf::st_as_sf(tibble_out, coords = c(".xcoord", ".ycoord"), crs = x$crs)
  } else {
    if (missing(local)) local <- NULL
    newdata_name <- newdata
    if (newdata_name == "all") {
      newdata_name <- names(x$ssn.object$preds)
    }
    tibble_out <- lapply(newdata_name, function(y) {
      newdata <- x$ssn.object$preds[[y]]
      if (NROW(newdata) == 0) {
        return(NULL)
      }
      preds_newdata <- predict(x,
        newdata = y, se.fit = se_fit, interval = interval,
        level = level, local = local, ...
      )
      if (se_fit) {
        if (interval %in% c("confidence", "prediction")) {
          augment_newdata <- tibble::tibble(
            .fitted = preds_newdata$fit[, "fit"]
          )
          augment_newdata$.lower <- preds_newdata$fit[, "lwr"]
          augment_newdata$.upper <- preds_newdata$fit[, "upr"]
        } else {
          augment_newdata <- tibble::tibble(
            .fitted = preds_newdata$fit
          )
        }
        augment_newdata$.se.fit <- preds_newdata$se.fit
      } else {
        if (interval %in% c("confidence", "prediction")) {
          augment_newdata <- tibble::tibble(
            .fitted = preds_newdata[, "fit"]
          )
          augment_newdata$.lower <- preds_newdata[, "lwr"]
          augment_newdata$.upper <- preds_newdata[, "upr"]
        } else {
          augment_newdata <- tibble::tibble(
            .fitted = preds_newdata
          )
        }
      }
      coords <- sf::st_coordinates(newdata)
      tibble_out <- tibble::tibble(cbind(sf::st_drop_geometry(newdata), augment_newdata))
      augment_newdata$pid <- ssn_get_netgeom(newdata, netvars = "pid")$pid
      tibble_out$.xcoord <- coords[, 1, drop = TRUE]
      tibble_out$.ycoord <- coords[, 2, drop = TRUE]
      tibble_out <- sf::st_as_sf(tibble_out, coords = c(".xcoord", ".ycoord"), crs = x$crs)
    })

    names(tibble_out) <- newdata_name
    if (length(tibble_out) == 1) {
      tibble_out <- tibble_out[[1]] # unlist if only one element
    }
  }
  tibble_out
}

#' @param type.predict The scale (\code{response} or \code{link}) of fitted
#'   values and predictions obtained using \code{ssn_glm()} objects.
#' @param type.residuals The residual type (\code{deviance}, \code{pearson}, or \code{response})
#'   of fitted models from \code{ssn_glm()} objects. Ignored if
#'   \code{newdata} is specified.
#' @param newdata_size The \code{size} value for each observation in \code{newdata}
#'   used when predicting for the binomial family.
#' @param var_correct A logical indicating whether to return the corrected prediction
#'   variances when predicting via models fit using \code{ssn_glm}. The default is
#'   \code{TRUE}.
#' @rdname augment.SSN2
#' @method augment ssn_glm
#' @export
augment.ssn_glm <- function(x, drop = TRUE, newdata = NULL, type.predict = c("link", "response"),
                            type.residuals = c("deviance", "pearson", "response"), se_fit = FALSE,
                            interval = c("none", "confidence", "prediction"),
                            newdata_size, level = 0.95, local = local, var_correct = TRUE, ...) {

  type.predict <- match.arg(type.predict)
  type.residuals <- match.arg(type.residuals)
  interval <- match.arg(interval)

  # set data and newdata
  if (is.null(newdata)) {
    if (drop) {
      data <- model.frame(x)
    } else {
      data <- sf::st_drop_geometry(x$ssn.object$obs)
    }
  } else {
    data <- model.frame(x)
  }

  if (is.null(newdata)) {
    augment_data <- tibble::tibble(.fitted = fitted(x, type = type.predict))
    if (se_fit) {
      preds_data <- predict(x, newdata = data, type = type.predict, se.fit = se_fit, interval = "confidence", ...)
      augment_data$.se.fit <- preds_data$se.fit
    }
    tibble_out <- tibble::tibble(cbind(data, augment_data, influence(x, type = type.residuals)))
    coords <- sf::st_coordinates(x$ssn.object$obs)
    tibble_out$.xcoord <- coords[, 1, drop = TRUE]
    tibble_out$.ycoord <- coords[, 2, drop = TRUE]
    tibble_out <- sf::st_as_sf(tibble_out, coords = c(".xcoord", ".ycoord"), crs = x$crs)
  } else {
    if (missing(newdata_size)) newdata_size <- NULL
    if (missing(local)) local <- NULL
    newdata_name <- newdata
    if (newdata_name == "all") {
      newdata_name <- names(x$ssn.object$preds)
    }
    if (missing(newdata_size)) newdata_size <- NULL
    tibble_out <- lapply(newdata_name, function(y) {
      newdata <- x$ssn.object$preds[[y]]
      if (NROW(newdata) == 0) {
        return(NULL)
      }
      # need to update newdata size in case there is one for each prediction object
      preds_newdata <- predict(x,
        newdata = y, type = type.predict, se.fit = se_fit, interval = interval,
        newdata_size = newdata_size, level = level,
        var_correct = FALSE, local = local, ...
      )
      if (se_fit) {
        if (interval %in% c("confidence", "prediction")) {
          augment_newdata <- tibble::tibble(
            .fitted = preds_newdata$fit[, "fit"]
          )
          augment_newdata$.lower <- preds_newdata$fit[, "lwr"]
          augment_newdata$.upper <- preds_newdata$fit[, "upr"]
        } else {
          augment_newdata <- tibble::tibble(
            .fitted = preds_newdata$fit
          )
        }
        augment_newdata$.se.fit <- preds_newdata$se.fit
      } else {
        if (interval %in% c("confidence", "prediction")) {
          augment_newdata <- tibble::tibble(
            .fitted = preds_newdata[, "fit"]
          )
          augment_newdata$.lower <- preds_newdata[, "lwr"]
          augment_newdata$.upper <- preds_newdata[, "upr"]
        } else {
          augment_newdata <- tibble::tibble(
            .fitted = preds_newdata
          )
        }
      }
      coords <- sf::st_coordinates(newdata)
      tibble_out <- tibble::tibble(cbind(sf::st_drop_geometry(newdata), augment_newdata))
      tibble_out$.xcoord <- coords[, 1, drop = TRUE]
      tibble_out$.ycoord <- coords[, 2, drop = TRUE]
      tibble_out <- sf::st_as_sf(tibble_out, coords = c(".xcoord", ".ycoord"), crs = x$crs)
    })

    names(tibble_out) <- newdata_name
    if (length(tibble_out) == 1) {
      tibble_out <- tibble_out[[1]] # unlist if only one element
    }
  }
  tibble_out
}
