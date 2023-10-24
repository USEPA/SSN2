#' Tidy a fitted model object
#'
#' @description Tidy a fitted model object into a summarized tibble.
#'
#' @param x A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param conf.int Logical indicating whether or not to include a confidence interval
#'   in the tidied output. The default is \code{FALSE}.
#' @param conf.level The confidence level to use for the confidence interval if
#'   \code{conf.int} is \code{TRUE}. Must be strictly greater than 0 and less than 1.
#'   The default is 0.95, which corresponds to a 95 percent confidence interval.
#' @param effects The type of effects to tidy. Available options are \code{"fixed"}
#'   (fixed effects), \code{"tailup"} (tailup covariance parameters),
#'   \code{"taildown"} (taildown covariance parameters), \code{"euclid"} (Euclidean
#'   covariance parameters), \code{"nugget"} (nugget covariance parameter),
#'   \code{"dispersion"} (dispersion parameter if relevant), \code{"ssn"} for all
#'   of \code{"tailup"}, \code{"taildown"}, \code{"euclid"}, \code{"nugget"}, and
#'   \code{"dispersion"}, and \code{"randcov"} (random effect variances). The default is \code{"fixed"}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @return A tidy tibble of summary information \code{effects}.
#'
#' @name tidy.SSN2
#' @method tidy ssn_lm
#' @export
#'
#' @seealso [glance.SSN2()] [augment.SSN2()]
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' tidy(ssn_mod)
tidy.ssn_lm <- function(x, conf.int = FALSE,
                        conf.level = 0.95, effects = "fixed", ...) {

  if (conf.int && (conf.level < 0 || conf.level > 1)) {
    stop("conf.level must be between 0 and 1.", call. = FALSE)
  }

  if (effects == "fixed") {
    result <- tibble::as_tibble(summary(x)$coefficients$fixed,
      rownames = "term", .name_repair = "minimal"
    )
    colnames(result) <- c(
      "term", "estimate", "std.error",
      "statistic", "p.value"
    )

    if (conf.int) {
      ci <- tibble::as_tibble(
        confint(x,
          level = conf.level,
          type = "fixed"
        ),
        rownames = "term", .name_repair = "minimal"
      )
      colnames(ci) <- c("term", "conf.low", "conf.high")
      result <- tibble::as_tibble(base::merge(result, ci, by = "term"), .name_repair = "minimal")
    }
  } else if (effects == "ssn") {
    tailup_coef <- unclass(coefficients(x, type = "tailup"))
    tailup_result <- tibble::tibble(
      effect = "tailup",
      term = names(tailup_coef),
      estimate = unname(tailup_coef),
      is_known = unname(x$is_known$tailup)
    )

    taildown_coef <- unclass(coefficients(x, type = "taildown"))
    taildown_result <- tibble::tibble(
      effect = "taildown",
      term = names(taildown_coef),
      estimate = unname(taildown_coef),
      is_known = unname(x$is_known$taildown)
    )

    euclid_coef <- unclass(coefficients(x, type = "euclid"))
    euclid_result <- tibble::tibble(
      effect = "euclid",
      term = names(euclid_coef),
      estimate = unname(euclid_coef),
      is_known = unname(x$is_known$euclid)
    )

    if (!x$anisotropy) {
      which_rotate <- which(euclid_result$term == "rotate")
      which_scale <- which(euclid_result$term == "scale")
      euclid_result <- euclid_result[-c(which_rotate, which_scale), , drop = FALSE]
    }

    nugget_coef <- unclass(coefficients(x, type = "nugget"))
    nugget_result <- tibble::tibble(
      effect = "nugget",
      term = names(nugget_coef),
      estimate = unname(nugget_coef),
      is_known = unname(x$is_known$nugget)
    )

    result <- tibble::as_tibble(rbind(tailup_result, taildown_result, euclid_result, nugget_result), .name_repair = "minimal")
  } else if (effects == "tailup") {
    tailup_coef <- unclass(coefficients(x, type = "tailup"))
    result <- tibble::tibble(
      term = names(tailup_coef),
      estimate = unname(tailup_coef),
      is_known = unname(x$is_known$tailup)
    )
  } else if (effects == "taildown") {
    taildown_coef <- unclass(coefficients(x, type = "taildown"))
    result <- tibble::tibble(
      term = names(taildown_coef),
      estimate = unname(taildown_coef),
      is_known = unname(x$is_known$taildown)
    )
  } else if (effects == "euclid") {
    euclid_coef <- unclass(coefficients(x, type = "euclid"))
    result <- tibble::tibble(
      term = names(euclid_coef),
      estimate = unname(euclid_coef),
      is_known = unname(x$is_known$euclid)
    )

    if (!x$anisotropy) {
      which_rotate <- which(result$term == "rotate")
      which_scale <- which(result$term == "scale")
      result <- result[-c(which_rotate, which_scale), , drop = FALSE]
    }
  } else if (effects == "nugget") {
    nugget_coef <- unclass(coefficients(x, type = "nugget"))
    result <- tibble::tibble(
      term = names(nugget_coef),
      estimate = unname(nugget_coef),
      is_known = unname(x$is_known$nugget)
    )
  } else if (effects == "randcov") {
    randcov_coef <- unclass(coefficients(x, type = "randcov")) # not needed now
    if (is.null(randcov_coef)) {
      result <- NULL
    } else {
      result <- tibble::tibble(
        term = names(randcov_coef),
        estimate = unname(randcov_coef),
        is_known = unname(x$is_known$randcov)
      )
    }
  }
  result
}





#' @rdname tidy.SSN2
#' @method tidy ssn_glm
#' @export
tidy.ssn_glm <- function(x, conf.int = FALSE,
                         conf.level = 0.95, effects = "fixed", ...) {

  if (conf.int && (conf.level < 0 || conf.level > 1)) {
    stop("conf.level must be between 0 and 1.", call. = FALSE)
  }

  if (effects == "fixed") {
    result <- tibble::as_tibble(summary(x)$coefficients$fixed,
      rownames = "term", .name_repair = "minimal"
    )
    colnames(result) <- c(
      "term", "estimate", "std.error",
      "statistic", "p.value"
    )

    if (conf.int) {
      ci <- tibble::as_tibble(
        confint(x,
          level = conf.level,
          type = "fixed"
        ),
        rownames = "term", .name_repair = "minimal"
      )
      colnames(ci) <- c("term", "conf.low", "conf.high")
      result <- tibble::as_tibble(base::merge(result, ci, by = "term"), .name_repair = "minimal")
    }
  } else if (effects == "ssn") {
    tailup_coef <- unclass(coefficients(x, type = "tailup"))
    tailup_result <- tibble::tibble(
      effect = "tailup",
      term = names(tailup_coef),
      estimate = unname(tailup_coef),
      is_known = unname(x$is_known$tailup)
    )

    taildown_coef <- unclass(coefficients(x, type = "taildown"))
    taildown_result <- tibble::tibble(
      effect = "taildown",
      term = names(taildown_coef),
      estimate = unname(taildown_coef),
      is_known = unname(x$is_known$taildown)
    )

    euclid_coef <- unclass(coefficients(x, type = "euclid"))
    euclid_result <- tibble::tibble(
      effect = "euclid",
      term = names(euclid_coef),
      estimate = unname(euclid_coef),
      is_known = unname(x$is_known$euclid)
    )

    if (!x$anisotropy) {
      which_rotate <- which(euclid_result$term == "rotate")
      which_scale <- which(euclid_result$term == "scale")
      euclid_result <- euclid_result[-c(which_rotate, which_scale), , drop = FALSE]
    }

    nugget_coef <- unclass(coefficients(x, type = "nugget"))
    nugget_result <- tibble::tibble(
      effect = "nugget",
      term = names(nugget_coef),
      estimate = unname(nugget_coef),
      is_known = unname(x$is_known$nugget)
    )

    dispersion_coef <- unclass(coefficients(x, type = "dispersion"))
    dispersion_result <- tibble::tibble(
      effect = "dispersion",
      term = names(dispersion_coef),
      estimate = unname(dispersion_coef),
      is_known = unname(x$is_known$dispersion)
    )

    result <- tibble::as_tibble(rbind(
      tailup_result, taildown_result, euclid_result, nugget_result,
      dispersion_result
    ), .name_repair = "minimal")
  } else if (effects == "tailup") {
    tailup_coef <- unclass(coefficients(x, type = "tailup"))
    result <- tibble::tibble(
      term = names(tailup_coef),
      estimate = unname(tailup_coef),
      is_known = unname(x$is_known$tailup)
    )
  } else if (effects == "taildown") {
    taildown_coef <- unclass(coefficients(x, type = "taildown"))
    result <- tibble::tibble(
      term = names(taildown_coef),
      estimate = unname(taildown_coef),
      is_known = unname(x$is_known$taildown)
    )
  } else if (effects == "euclid") {
    euclid_coef <- unclass(coefficients(x, type = "euclid"))
    result <- tibble::tibble(
      term = names(euclid_coef),
      estimate = unname(euclid_coef),
      is_known = unname(x$is_known$euclid)
    )

    if (!x$anisotropy) {
      which_rotate <- which(result$term == "rotate")
      which_scale <- which(result$term == "scale")
      result <- result[-c(which_rotate, which_scale), , drop = FALSE]
    }
  } else if (effects == "nugget") {
    nugget_coef <- unclass(coefficients(x, type = "nugget"))
    result <- tibble::tibble(
      term = names(nugget_coef),
      estimate = unname(nugget_coef),
      is_known = unname(x$is_known$nugget)
    )
  } else if (effects == "dispersion") {
    dispersion_coef <- unclass(coefficients(x, type = "dispersion"))
    result <- tibble::tibble(
      term = names(dispersion_coef),
      estimate = unname(dispersion_coef),
      is_known = unname(x$is_known$dispersion)
    )
  } else if (effects == "randcov") {
    randcov_coef <- unclass(coefficients(x, type = "randcov")) # not needed now
    if (is.null(randcov_coef)) {
      result <- NULL
    } else {
      result <- tibble::tibble(
        term = names(randcov_coef),
        estimate = unname(randcov_coef),
        is_known = unname(x$is_known$randcov)
      )
    }
  }
  result
}
