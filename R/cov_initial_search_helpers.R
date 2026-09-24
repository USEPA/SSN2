#' Cross a spatial covariance-starting-value grid with random-effect allocations
#'
#' When random effects are present, expands \code{cov_grid} into three
#' regimes -- spatial-dominant (90% spatial variance, random effects share
#' 10%), random-dominant (10% spatial, random effects share 90%, one row per
#' effect where that effect alone dominates when there is more than one), and
#' an even 50/50 split -- matching spmodel's random-effect starting-value
#' normalization.
#'
#' @param cov_grid A spatial covariance-starting-value grid from
#'   \code{\link{cov_initial_design}()}/\code{\link{build_cov_initial_grid}()}.
#' @param initial_NA_object A joint covariance initial-value object; only
#'   used to check whether \code{randcov_initial} is present.
#' @param data_object A model data object with \code{randcov_names}.
#' @param ns2 The overall variance anchor (1.2 times the OLS residual
#'   variance).
#'
#' @return \code{cov_grid} unchanged if there are no random effects;
#'   otherwise the three-regime expanded grid, with spatial variance columns
#'   rescaled and one column added per random effect.
#'
#' @noRd
add_randcov_grid <- function(cov_grid, initial_NA_object, data_object, ns2) {
  if (is.null(initial_NA_object$randcov_initial)) {
    return(cov_grid)
  }

  # find the number of random effects
  nvar_randcov <- length(data_object$randcov_names)
  if (nvar_randcov == 1) {
    randcov_grid <- matrix(1, nrow = 1, ncol = 1)
    max_row <- 1
  } else {
    # create a grid of random effects
    randcov_grid <- matrix(0.1 / (nvar_randcov - 1), nrow = nvar_randcov, ncol = nvar_randcov)
    diag(randcov_grid) <- 0.9
    randcov_grid <- rbind(randcov_grid, matrix(1 / nvar_randcov, nrow = 1, ncol = nvar_randcov))
    max_row <- nvar_randcov + 1
  }
  randcov_grid <- as.data.frame(ns2 * randcov_grid)
  colnames(randcov_grid) <- data_object$randcov_names

  # cov_grid1 focuses on spatial parameters
  cov_grid1 <- merge(cov_grid, 0.1 * randcov_grid[max_row, , drop = FALSE], by = NULL)
  cov_grid1[, c("tailup_de", "taildown_de", "euclid_de", "nugget")] <- 0.9 * cov_grid1[, c("tailup_de", "taildown_de", "euclid_de", "nugget")]

  # cov_grid2 focuses on random effects
  cov_grid2 <- merge(cov_grid, 0.9 * randcov_grid, by = NULL)
  cov_grid2[, c("tailup_de", "taildown_de", "euclid_de", "nugget")] <- 0.1 * cov_grid2[, c("tailup_de", "taildown_de", "euclid_de", "nugget")]

  # cov_grid3 is an even spread between spatial and random effects
  cov_grid3 <- merge(cov_grid, 0.5 * randcov_grid, by = NULL)
  cov_grid3[, c("tailup_de", "taildown_de", "euclid_de", "nugget")] <- 0.5 * cov_grid3[, c("tailup_de", "taildown_de", "euclid_de", "nugget")]

  # bind them together and rename as cov_grid to match case without random effects
  rbind(cov_grid1, cov_grid2, cov_grid3)
}

#' Build Euclidean range/shape starting values, expanding over shape when unknown
#'
#' @param euclid_initial The Euclidean covariance initial-value object.
#' @param euclid_range A distance (or vector of distances) to convert into a
#'   range starting value via \code{\link{get_cov_initial_range}()}.
#' @param euclid_max The Euclidean distance anchor (currently unused; kept
#'   for a call signature consistent with related helpers).
#' @param extra An optional vector of shape (\code{"extra"}) starting values
#'   to use directly, recycled to \code{length(euclid_range)}; if
#'   \code{NULL}, family-specific defaults are used (or the supplied known
#'   value, if any).
#'
#' @return A list with \code{range} and \code{extra} (\code{NULL} for
#'   Euclidean families without a shape parameter).
#'
#' @noRd
get_euclid_start_values <- function(euclid_initial, euclid_range, euclid_max, extra = NULL) {
  euclid_type <- remove_covtype(class(euclid_initial)[1])
  if (!euclid_has_extra(euclid_type)) {
    return(list(range = get_cov_initial_range(euclid_initial, euclid_range), extra = NULL))
  }

  if (is.null(extra)) {
    extra <- rep(switch(euclid_type, matern = c(1, 4), cauchy = c(0.5, 2),
                        pexponential = c(0.4, 1.6)), each = 2)
  }
  extra <- rep(extra, length.out = length(euclid_range))
  if ("extra" %in% names(euclid_initial$initial) &&
      !is.na(euclid_initial$initial[["extra"]])) {
    extra <- rep(euclid_initial$initial[["extra"]], length(euclid_range))
  }
  list(range = get_cov_initial_range(euclid_initial, euclid_range), extra = extra)
}

#' Extract the Euclidean starting values from the best grid candidate
#'
#' @param euclid_initial The Euclidean covariance initial-value object (used
#'   only for its type, to decide whether \code{extra} applies).
#' @param min_params A named vector/list (one row of a covariance-starting-value
#'   grid) with \code{euclid_de}, \code{euclid_range}, optionally
#'   \code{euclid_extra}, \code{rotate}, and \code{scale}.
#'
#' @return A named vector with \code{de}, \code{range}, (if applicable)
#'   \code{extra}, \code{rotate}, and \code{scale}.
#'
#' @noRd
update_euclid_grid_initial <- function(euclid_initial, min_params) {
  euclid_type <- remove_covtype(class(euclid_initial)[1])
  values <- c(de = min_params[["euclid_de"]], range = min_params[["euclid_range"]])
  if (euclid_has_extra(euclid_type)) {
    values <- c(values, extra = min_params[["euclid_extra"]])
  }
  c(values, rotate = min_params[["rotate"]], scale = min_params[["scale"]])
}

#' Pin already-known covariance-parameter values across every grid row
#'
#' For each covariance/random-effect field with a known (non-\code{NA})
#' value in \code{initial_object}, overwrites that field in every row of
#' \code{cov_grid} rather than letting the grid search vary it; used for
#' fields not specific to either the Gaussian or GLM starting-value search
#' (family-specific fields, such as GLM dispersion, are pinned by each
#' caller's own \code{cov_grid_replace}/\code{cov_grid_replace_glm}).
#'
#' @param cov_grid A covariance-starting-value grid.
#' @param initial_object A joint covariance initial-value object, with any
#'   known fields already resolved (non-\code{NA}).
#' @param data_object A model data object with \code{randcov_names}.
#'
#' @return \code{cov_grid} with known fields pinned.
#'
#' @noRd
cov_grid_replace_shared <- function(cov_grid, initial_object, data_object) {
  if (!is.na(initial_object$tailup_initial$initial[["de"]])) {
    cov_grid[, "tailup_de"] <- initial_object$tailup_initial$initial[["de"]]
  }

  if (!is.na(initial_object$tailup_initial$initial[["range"]])) {
    cov_grid[, "tailup_range"] <- initial_object$tailup_initial$initial[["range"]]
  }

  if (!is.na(initial_object$taildown_initial$initial[["de"]])) {
    cov_grid[, "taildown_de"] <- initial_object$taildown_initial$initial[["de"]]
  }

  if (!is.na(initial_object$taildown_initial$initial[["range"]])) {
    cov_grid[, "taildown_range"] <- initial_object$taildown_initial$initial[["range"]]
  }

  if (!is.na(initial_object$euclid_initial$initial[["de"]])) {
    cov_grid[, "euclid_de"] <- initial_object$euclid_initial$initial[["de"]]
  }

  if (!is.na(initial_object$euclid_initial$initial[["range"]])) {
    cov_grid[, "euclid_range"] <- initial_object$euclid_initial$initial[["range"]]
  }

  if ("euclid_extra" %in% names(cov_grid) &&
      "extra" %in% names(initial_object$euclid_initial$initial) &&
      !is.na(initial_object$euclid_initial$initial[["extra"]])) {
    cov_grid[, "euclid_extra"] <- initial_object$euclid_initial$initial[["extra"]]
  }

  if (!is.na(initial_object$euclid_initial$initial[["rotate"]])) {
    cov_grid[, "rotate"] <- initial_object$euclid_initial$initial[["rotate"]]
  }

  if (!is.na(initial_object$euclid_initial$initial[["scale"]])) {
    cov_grid[, "scale"] <- initial_object$euclid_initial$initial[["scale"]]
  }

  if (!is.na(initial_object$nugget_initial$initial[["nugget"]])) {
    cov_grid[, "nugget"] <- initial_object$nugget_initial$initial[["nugget"]]
  }

  for (x in data_object$randcov_names) {
    if (!is.na(initial_object$randcov_initial$initial[[x]])) {
      cov_grid[, x] <- initial_object$randcov_initial$initial[[x]]
    }
  }
  cov_grid
}
