#' Find initial values for optimization.
#'
#' @param initial_NA_object Initial object with relevant NAs to indicate which require estimation.
#' @param ssn.object SSN object.
#' @param data_object Data object.
#' @param estmethod Estimation method.
#'
#' @return Initial values
#' @noRd
cov_initial_search_glm <- function(initial_NA_object, ssn.object, data_object, estmethod) {
  # find the initial "new" sample variance
  s2 <- data_object$s2
  ns2 <- 1.2 * s2

  # create a grid of initial values for the covariance parameters
  var_min <- 0.05
  tailup_de <- pmax(ns2 * c(rep(9 / 10, 2), rep(1 / 3 * 1 / 10, 2), rep(1 / 3 * 1 / 10, 2), 1 / 3 * 1 / 10, rep(1 / 4, 2)), var_min)
  taildown_de <- pmax(ns2 * c(rep(1 / 3 * 1 / 10, 2), rep(9 / 10, 2), rep(1 / 3 * 1 / 10, 2), 1 / 3 * 1 / 10, rep(1 / 4, 2)), var_min)
  euclid_de <- pmax(ns2 * c(rep(1 / 3 * 1 / 10, 2), rep(1 / 3 * 1 / 10, 2), rep(9 / 10, 2), 1 / 3 * 1 / 10, rep(1 / 4, 2)), var_min)
  nugget <- pmax(ns2 * c(rep(1 / 3 * 1 / 10, 2), rep(1 / 3 * 1 / 10, 2), rep(1 / 3 * 1 / 10, 2), 9 / 10, rep(1 / 4, 2)), var_min)

  # find the maximum tail and euclidean distances to consider
  tail_max <- data_object$tail_max
  euclid_max <- data_object$euclid_max

  # find the initial ranges based on the maximum distances
  tailup_range <- c(1 / 4 * tail_max, 3 / 4 * tail_max, rep(1 / 4, 6) * tail_max, 3 / 4 * tail_max) # can be custom for each cov function
  taildown_range <- c(rep(1 / 4, 3) * tail_max, 3 / 4 * tail_max, rep(1 / 4, 4) * tail_max, 3 / 4 * tail_max)
  euclid_range <- c(rep(1 / 4, 5) * euclid_max, 3 / 4 * euclid_max, rep(1 / 4, 2) * euclid_max, 3 / 4 * euclid_max)


  # create a grid of initial covariance parameter values
  cov_grid <- data.frame(
    tailup_de = tailup_de, taildown_de = taildown_de, euclid_de = euclid_de,
    nugget = nugget, tailup_range = tailup_range, taildown_range = taildown_range,
    euclid_range = euclid_range, rotate = 0, scale = 1, dispersion = 1
  )

  cov_grid2 <- cov_grid
  cov_grid2$dispersion <- 100
  cov_grid <- rbind(cov_grid, cov_grid2)

  # if there is anisotropy test midpoints
  if (data_object$anisotropy) {
    cov_grid2 <- cov_grid
    cov_grid2$rotate <- pi / 2
    cov_grid2$scale <- 0.5
    cov_grid <- rbind(cov_grid, cov_grid2)
  }

  # if there are random effects add to the grid
  if (!is.null(initial_NA_object$randcov)) {
    # find the number of random effects
    nvar_randcov <- length(data_object$randcov_names)
    if (nvar_randcov == 1) {
      randcov_grid <- matrix(1, nrow = 1, ncol = 1)
      max_row <- 1
    } else {
      # create a grid of random effects
      randcov_grid <- matrix(0.1 / nvar_randcov, nrow = nvar_randcov, ncol = nvar_randcov)
      diag(randcov_grid) <- 0.9
      randcov_grid <- rbind(randcov_grid, matrix(1 / nvar_randcov, nrow = 1, ncol = nvar_randcov))
      max_row <- nvar_randcov + 1
    }
    randcov_grid <- as.data.frame(ns2 * randcov_grid)
    colnames(randcov_grid) <- data_object$randcov_names

    # cov_grid1 focuses on spatial parameters
    cov_grid1 <- merge(cov_grid, 0.1 * randcov_grid[max_row, , drop = FALSE], by = NULL)
    cov_grid1[, c("tailup_de", "taildown_de", "euclid_de", "nugget")] <- 0.9 * cov_grid1[, c("tailup_de", "taildown_de", "euclid_de", "nugget")]
    cov_grid2_index <- which(vapply(seq_len(NROW(cov_grid)), function(x) {
      length(unique(unlist(cov_grid[x, c("tailup_de", "taildown_de", "euclid_de", "nugget")]))) == 1 &
        cov_grid[x, "tailup_range"] == 1 / 4 * tail_max &
        cov_grid[x, "rotate"] == 0 &
        cov_grid[x, "scale"] == 1
    }, logical(1)))

    # cov_grid2 focuses on random effects
    cov_grid2 <- merge(cov_grid[cov_grid2_index, , drop = FALSE], 0.9 * randcov_grid[seq(1, max_row, by = 1), , drop = FALSE], by = NULL)
    cov_grid2[, c("tailup_de", "taildown_de", "euclid_de", "nugget")] <- 0.1 * cov_grid2[, c("tailup_de", "taildown_de", "euclid_de", "nugget")]

    # cov_grid3 is an even spread between spatial and random effects
    cov_grid3_index <- which(vapply(seq_len(NROW(cov_grid)), function(x) {
      length(unique(unlist(cov_grid[x, c("tailup_de", "taildown_de", "euclid_de", "nugget")]))) == 1 &
        cov_grid[x, "rotate"] == 0 &
        cov_grid[x, "scale"] == 1
    }, logical(1)))
    cov_grid3 <- merge(cov_grid[cov_grid3_index, , drop = FALSE], 0.5 * randcov_grid[seq(1, max_row, by = 1), , drop = FALSE], by = NULL)
    cov_grid3[, c("tailup_de", "taildown_de", "euclid_de", "nugget")] <- 0.5 * cov_grid3[, c("tailup_de", "taildown_de", "euclid_de", "nugget")]

    # bind them together and rename as cov_grid to match case without random effects
    cov_grid <- rbind(cov_grid1, cov_grid2, cov_grid3)
  }

  # find the unique rows of the grid
  cov_grid <- unique(cov_grid_replace_glm(cov_grid, initial_NA_object, data_object))
  # split into list
  cov_grid_splits <- split(cov_grid, seq_len(NROW(cov_grid)))
  # iterate through list
  objvals <- vapply(cov_grid_splits, function(x) eval_grid_glm(x, initial_NA_object, ssn.object, data_object, estmethod), numeric(1))
  # find parameters that yield the minimum -2ll
  min_params <- unlist(cov_grid_splits[[which.min(objvals)]])
  # store this as new NA object
  updated_NA_object <- initial_NA_object
  updated_NA_object$tailup_initial$initial <- c(de = min_params[["tailup_de"]], range = min_params[["tailup_range"]])
  updated_NA_object$taildown_initial$initial <- c(de = min_params[["taildown_de"]], range = min_params[["taildown_range"]])
  updated_NA_object$euclid_initial$initial <- c(
    de = min_params[["euclid_de"]], range = min_params[["euclid_range"]],
    rotate = min_params[["rotate"]], scale = min_params[["scale"]]
  )
  updated_NA_object$nugget_initial$initial <- c(nugget = min_params[["nugget"]])
  updated_NA_object$dispersion_initial$initial <- c(dispersion = min_params[["dispersion"]])

  if (!is.null(updated_NA_object$randcov_initial)) {
    updated_NA_object$randcov_initial$initial <- min_params[data_object$randcov_names]
  }


  # return best parameters
  best_params <- list(initial_object = updated_NA_object)
}

#' Evaluate initial values
#'
#' @param cov_grid_split A set of initial values.
#' @param initial_NA_object Initial object with relevant NAs to indicate which require estimation.
#' @param ssn.object SSN object.
#' @param data_object Data object.
#' @param estmethod Estimation method.
#'
#' @return The minus twice log likelihood
#' @noRd
eval_grid_glm <- function(cov_grid_split, initial_NA_object, ssn.object, data_object, estmethod) {
  cov_grid <- unlist(cov_grid_split)
  # params object
  params_object <- get_params_object_grid_glm(cov_grid, initial_NA_object)

  # gloglik products
  lapll_prods <- laploglik_products(params_object, data_object, estmethod)

  # minus two gloglik
  get_minustwolaploglik(lapll_prods, data_object, estmethod)
}

#' Replace initial values with known values.
#'
#' @param cov_grid A grid of initial values.
#' @param initial_object An initial object.
#' @param data_object Data object.
#'
#' @return Initial values with known values replaced.
#' @noRd
cov_grid_replace_glm <- function(cov_grid, initial_object, data_object) {
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

  if (!is.na(initial_object$euclid_initial$initial[["rotate"]])) {
    cov_grid[, "rotate"] <- initial_object$euclid_initial$initial[["rotate"]]
  }

  if (!is.na(initial_object$euclid_initial$initial[["scale"]])) {
    cov_grid[, "scale"] <- initial_object$euclid_initial$initial[["scale"]]
  }

  if (!is.na(initial_object$nugget_initial$initial[["nugget"]])) {
    cov_grid[, "nugget"] <- initial_object$nugget_initial$initial[["nugget"]]
  }

  if (!is.na(initial_object$dispersion_initial$initial[["dispersion"]])) {
    cov_grid[, "dispersion"] <- initial_object$dispersion_initial$initial[["dispersion"]]
  }

  for (x in data_object$randcov_names) {
    if (!is.na(initial_object$randcov_initial$initial[[x]])) {
      cov_grid[, x] <- initial_object$randcov_initial$initial[[x]]
    }
  }

  cov_grid
}
