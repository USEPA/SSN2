get_randcov_list <- function(index, randcov_Zs = NULL, randcov_names = NULL) {
  if (is.null(randcov_names)) {
    randcov_list <- NULL
  } else {
    index_vals <- split(seq_along(index), index)
    randcov_list <- lapply(index_vals, function(x) get_randcov_index_list(x, index, randcov_Zs, randcov_names))
    names(randcov_list) <- names(index_vals)
  }
  randcov_list
}

get_randcov_index_list <- function(index_val, index, randcov_Zs, randcov_names) {
  Z_lists <- lapply(randcov_names, function(x) get_randcov_var_list(x, index_val, index, randcov_Zs))
  names(Z_lists) <- randcov_names
  Z_lists
}

get_randcov_var_list <- function(randcov_name, index_val, index, randcov_Zs) {
  Z_list <- randcov_Zs[[randcov_name]][["Z"]][index_val, , drop = FALSE]
  if (is.null(randcov_Zs[[randcov_name]][["ZZt"]])) {
    ZZt_list <- NULL
  } else {
    ZZt_list <- randcov_Zs[[randcov_name]][["ZZt"]][index_val, index_val, drop = FALSE]
  }
  if (is.null(randcov_Zs[[randcov_name]][["ZtZ"]])) {
    ZtZ_list <- NULL
  } else {
    ZtZ_list <- randcov_Zs[[randcov_name]][["ZtZ"]][index_val, index_val, drop = FALSE]
  }
  list(Z = Z_list, ZZt = ZZt_list, ZtZ = ZtZ_list)
}

#' Get random effect design matrices
#'
#' @param data Data
#' @param randcov_names Random effect names
#' @param ZZt Product of random effect design matrix and its transpose
#' @param ZtZ Transpose of ZZt
#' @param xlev_list Levels of random effects
#'
#' @noRd
get_randcov_Zs <- function(data, randcov_names = NULL, ZZt = TRUE, ZtZ = FALSE, xlev_list = NULL) {
  if (is.null(randcov_names)) {
    randcov_Zs <- NULL
  } else {
    randcov_Zs <- lapply(randcov_names, get_randcov_Z, data, ZZt, ZtZ, xlev_list)
    names(randcov_Zs) <- randcov_names
  }
  randcov_Zs
}

#' Get a random effect design matrix
#'
#' @param randcov_name Random effect name
#' @param randcov_names Random effect names
#' @param ZZt Product of random effect design matrix and its transpose
#' @param ZtZ Transpose of ZZt
#' @param xlev_list Levels of random effects
#'
#' @noRd
get_randcov_Z <- function(randcov_name, data, ZZt = TRUE, ZtZ = FALSE, xlev_list = NULL) {
  bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
  Z_reform <- reformulate(bar_split[[2]], intercept = FALSE)
  if (is.null(xlev_list)) {
    Z_frame <- model.frame(Z_reform, data = data, drop.unused.levels = FALSE)
  } else {
    Z_frame <- model.frame(Z_reform, data = data, drop.unused.levels = FALSE, xlev = xlev_list[[randcov_name]])
  }
  if (any(!attr(terms(Z_frame), "dataClasses") %in% c("character", "factor", "ordered"))) {
    stop("Random effect grouping variables must be categorical or factor.", call. = FALSE)
  }
  # Compound grouping terms need every model-frame column, not just the first.
  if (NCOL(Z_frame) == 1) {
    Z_factor <- Z_frame[[1]]
  } else {
    Z_factor <- interaction(Z_frame, drop = FALSE, sep = ":")
  }
  Z_index <- t(fac2sparse(Z_factor, drop.unused.levels = FALSE))
  if (bar_split[[1]] == "1") {
    Z <- Z_index
  } else {
    Z_mod_reform <- reformulate(bar_split[[1]], intercept = FALSE)
    Z_mod_frame <- model.frame(Z_mod_reform, data = data, drop.unused.levels = FALSE)
    Z_mod <- model.matrix(Z_mod_reform, Z_mod_frame)
    if (NCOL(Z_mod) > 1) {
      stop("All variable names to the left of | in random must be numeric.", call. = FALSE)
    }
    Z <- as.vector(Z_mod) * Z_index
  }
  # find and replace values not observed
  Z_levels_observed <- which(colSums(abs(Z)) > 0)
  Z <- Z[, Z_levels_observed, drop = FALSE]
  if (ZZt) {
    ZZt <- tcrossprod(Z, Z)
  } else {
    ZZt <- NULL
  }

  if (ZtZ) {
    ZtZ <- crossprod(Z, Z)
  } else {
    ZtZ <- NULL
  }

  list(Z = Z, ZZt = ZZt, ZtZ = ZtZ)
}

#' Find the names of random effects and coerce them (if needed) to consistent
#'   structure (1 | ranef) or (x | ranef)
#'
#' @param random  A random effects formula (one sided e.g., ~ random effects)
#'   where random intercepts are specified by ~ group or ~ (1 | group) and
#'   random slosep are specified by ~ (x | group)
#'
#' @return Names of random effects
#'
#' @noRd
get_randcov_names <- function(random = NULL) {
  if (is.null(random)) {
    new_labels <- NULL
  } else {
    # get the random formula and turn it into a named vector here
    labels_initial <- labels(terms(random))
    new_labels <- unlist(lapply(labels_initial, get_randcov_name))
  }
  new_labels
}

#' Adjust a single name of a random effect
#'
#' @param label Name of a random effect
#'
#' @noRd
get_randcov_name <- function(label) {
  if (grepl("|", label, fixed = TRUE)) {
    new_label <- label
  } else {
    new_label <- paste("1", label, sep = " | ")
  }
  if (grepl("/", new_label, fixed = TRUE)) {
    bar_split <- unlist(strsplit(new_label, " | ", fixed = TRUE))
    dash_split <- unlist(strsplit(bar_split[[2]], "/", fixed = TRUE))
    front <- bar_split[[1]]
    backs <- dash_split
    new_label <- lapply(seq_along(backs), function(x) paste(front, paste(backs[seq(from = 1, to = x, by = 1)], collapse = ":"), sep = " | "))
  }
  new_label <- unlist(lapply(new_label, function(x) get_randcov_label(x)))
  new_label
}

get_randcov_label <- function(label) {
  strsplits <- strsplit(label, " | ", fixed = TRUE)
  terms_fronts <- terms(reformulate(strsplits[[1]][[1]]))
  labels_fronts <- labels(terms_fronts)
  if (attr(terms_fronts, "intercept") == 1) {
    labels_fronts <- c("1", labels_fronts)
  }
  form_fronts <- lapply(labels_fronts, function(x) paste(x, strsplits[[1]][[2]], sep = " | "))
}

extend_randcov_xlev <- function(random_xlev, newdata, randcov_names) {
  if (is.null(random_xlev) || is.null(randcov_names)) {
    return(random_xlev)
  }
  extended <- lapply(randcov_names, function(randcov_name) {
    old_xlev <- random_xlev[[randcov_name]]
    if (is.null(old_xlev)) {
      return(NULL)
    }
    bar_split <- unlist(strsplit(randcov_name, " | ", fixed = TRUE))
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
    newdata_mf <- model.frame(reform_bar2, newdata, na.action = na.pass)
    newdata_xlev <- .getXlevels(terms(newdata_mf), newdata_mf)
    mapply(function(old, var) union(old, newdata_xlev[[var]]), old_xlev, names(old_xlev), SIMPLIFY = FALSE)
  })
  names(extended) <- randcov_names
  extended
}

get_randcov_xlev <- function(randcov_names, obdata) {
  randcov_xlev <- lapply(randcov_names, function(x) {
    bar_split <- unlist(strsplit(x, " | ", fixed = TRUE))
    reform_bar2 <- reformulate(bar_split[[2]], intercept = FALSE)
    mf <- model.frame(reform_bar2, obdata)
    .getXlevels(terms(mf), mf)
  })
  names(randcov_xlev) <- randcov_names
  randcov_xlev
}

randcov_matrix <- function(randcov_params = NULL, randcov_list) {
  if (is.null(randcov_params)) {
    randcov_matrix_val <- NULL
  } else {
    randcov_names <- names(randcov_params)
    # var times ZZt
    randcov_matrices <- lapply(randcov_names, function(x) randcov_params[[x]] * randcov_list[[x]][["ZZt"]])
    randcov_matrix_val <- Reduce("+", randcov_matrices)
  }
  randcov_matrix_val
}
