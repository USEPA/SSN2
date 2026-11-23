NULL

#' Build a model matrix and frame, dropping rows with missing predictors
#'
#' Resolves \code{formula}'s data (including a bare \code{.}, via
#' \code{\link{get_dot_formula_data}()}), builds the model matrix, drops any
#' rows with a missing predictor value (erroring is handled elsewhere; this
#' just excludes them so the response's own missingness can still be used to
#' define prediction rows), then rebuilds the model frame/matrix on the
#' cleaned data and checks for perfect collinearity.
#'
#' @param formula A model formula.
#' @param obdata A data frame (or \code{sf} object) to build the model
#'   matrix from.
#' @param sf_column_name The geometry column name, or \code{NULL}.
#' @param ... Additional arguments, such as \code{contrasts}.
#'
#' @return A list with \code{obdata} (rows with missing predictors removed),
#'   \code{obdata_model_frame}, \code{terms_val}, \code{X} (the model
#'   matrix), \code{dots} (with \code{contrasts} resolved), \code{xlevels},
#'   \code{p} (matrix rank), \code{n} (row count), and \code{formula} (with a
#'   bare \code{.} expanded).
#'
#' @noRd
get_model_matrix_object <- function(formula, obdata, sf_column_name = NULL, ...) {
  # data used to resolve a bare "." in formula, excluding SSN geometry/topology
  # metadata columns unless they are explicitly named elsewhere in formula
  dot_data <- get_dot_formula_data(formula, obdata, sf_column_name)

  # finding model frame
  obdata_model_frame <- model.frame(formula, dot_data, drop.unused.levels = TRUE, na.action = na.pass)
  # finding contrasts as ...
  dots <- list(...)
  if (!"contrasts" %in% names(dots)) dots$contrasts <- NULL

  # model matrix with potential NA
  X <- model.matrix(formula, obdata_model_frame, contrasts = dots$contrasts)
  # finding rows w/out NA
  ob_predictors <- complete.cases(X)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors.", call. = FALSE)
  }
  # subset obdata (and its dot_data counterpart) by nonNA predictors
  obdata <- obdata[ob_predictors, , drop = FALSE]
  dot_data <- dot_data[ob_predictors, , drop = FALSE]

  # new model frame
  obdata_model_frame <- model.frame(formula, dot_data, drop.unused.levels = TRUE, na.action = na.omit)
  # find terms
  terms_val <- terms(obdata_model_frame)
  if ("." %in% all.vars(formula)) {
    formula <- formula(terms_val)
  }
  # find X
  X <- model.matrix(formula, obdata_model_frame, contrasts = dots$contrasts)
  # find induced contrasts and xlevels
  dots$contrasts <- attr(X, "contrasts")
  xlevels <- .getXlevels(terms_val, obdata_model_frame)
  # find p
  p <- as.numeric(Matrix::rankMatrix(X, method = "qr"))
  if (p < NCOL(X)) {
    warning("There are perfect collinearities detected in X (the matrix of explanatory variables). This may make the model fit unreliable or may cause an error while model fitting. Consider removing redundant explanatory variables and refitting the model.", call. = FALSE)
  }
  # find sample size
  n <- NROW(X)

  list(
    obdata = obdata, obdata_model_frame = obdata_model_frame, terms_val = terms_val,
    X = X, dots = dots, xlevels = xlevels, p = p, n = n, formula = formula
  )
}

#' Resolve the data a bare \code{.} in a formula should expand against
#'
#' When \code{formula} contains no bare \code{.}, returns \code{obdata}
#' unchanged. Otherwise drops the geometry column (unless it is itself named
#' elsewhere in \code{formula}) and any SSN topology/metadata columns
#' (\code{netgeom}, \code{netID}, \code{rid}, \code{upDist}, \code{ratio},
#' \code{pid}, \code{locID}) not explicitly named elsewhere in
#' \code{formula}, so \code{.} expands only to genuine candidate predictors.
#'
#' @param formula A model formula.
#' @param obdata A data frame (or \code{sf} object).
#' @param sf_column_name The geometry column name, or \code{NULL}.
#'
#' @return \code{obdata}, with irrelevant columns dropped when \code{formula}
#'   contains a bare \code{.}.
#'
#' @noRd
get_dot_formula_data <- function(formula, obdata, sf_column_name = NULL) {
  if (!"." %in% all.vars(formula)) {
    return(obdata)
  }
  referenced <- setdiff(all.vars(formula), ".")

  if (!is.null(sf_column_name) && !(sf_column_name %in% referenced) &&
    sf_column_name %in% colnames(obdata)) {
    obdata <- sf::st_drop_geometry(obdata)
  }

  topology_cols <- c("netgeom", "netID", "rid", "upDist", "ratio", "pid", "locID")
  topology_cols <- intersect(topology_cols, colnames(obdata))
  drop_cols <- setdiff(topology_cols, referenced)
  if (length(drop_cols) > 0) {
    obdata <- obdata[, setdiff(colnames(obdata), drop_cols), drop = FALSE]
  }
  obdata
}

#' Check that a supplied \code{local$index} matches the fitted sample size
#'
#' @param index The supplied \code{local$index} grouping vector, or
#'   \code{NULL}.
#' @param n The number of non-missing (fitted) response observations.
#'
#' @return \code{NULL}, invisibly, if \code{index} is \code{NULL} or has
#'   length \code{n}; otherwise an error.
#'
#' @noRd
check_local_index_length <- function(index, n) {
  if (!is.null(index) && length(index) != n) {
    stop(
      "local$index must have the same length as the non-missing (non-NA) response vector (",
      n, "), but has length ", length(index), ". ",
      "Observations with a missing response are excluded from fitting (they become prediction locations instead), so they must also be excluded from local$index.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Check that a set of model input columns has no missing values
#'
#' @param vars A character vector of column names to check (columns not
#'   present in \code{data} are ignored).
#' @param data A data frame.
#' @param label A short label for the error message (e.g.
#'   \code{"partition_factor"}, \code{"random effect"}).
#'
#' @return \code{NULL}, invisibly, if no checked column has a missing value;
#'   otherwise an error naming the offending column(s).
#'
#' @noRd
check_no_na_columns <- function(vars, data, label) {
  vars <- intersect(vars, colnames(data))
  if (length(vars) == 0) {
    return(invisible(NULL))
  }
  na_counts <- vapply(vars, function(x) sum(is.na(data[[x]])), numeric(1))
  bad <- names(na_counts)[na_counts > 0]
  if (length(bad) > 0) {
    stop(
      "Missing values found in ", label, " variable(s): ", paste0("'", bad, "'", collapse = ", "), ". ",
      "Rows with a missing value in a modeling input (other than the response) are not allowed; remove or impute them before fitting.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Check that a response vector is numeric and has nonzero variance
#'
#' @param y The response vector.
#'
#' @return \code{NULL}, invisibly, if \code{y} is a numeric vector with
#'   variance greater than zero; otherwise an error.
#'
#' @noRd
check_response_numeric_and_variable <- function(y) {
  # see if response is numeric
  if (!is.numeric(y)) {
    stop("Response variable must be numeric", call. = FALSE)
  }

  # error if no variance
  if (var(y) == 0) {
    stop("The response has no variability. Model fit unreliable.", call. = FALSE)
  }
}

#' Check that the number of fixed effects is smaller than the sample size
#'
#' @param p The number of fixed effects (model matrix rank).
#' @param n The sample size.
#' @param refit_fun_text The name of the fitting function to name in the error
#'   message (e.g. \code{"ssn_lm"} or \code{"ssn_glm"}), matching whichever
#'   function's data object is being built.
#'
#' @return \code{NULL}, invisibly, if \code{p < n}; otherwise an error.
#'
#' @noRd
check_p_lt_n <- function(p, n, refit_fun_text = "ssn_lm") {
  # error if p >= n
  if (p >= n) {
    stop("The number of fixed effects is at least as large as the number of observations (p >= n). Consider reducing the number of fixed effects and rerunning ", refit_fun_text, "().", call. = FALSE)
  }
}

#' Validate and normalize a \code{partition_factor} formula
#'
#' Checks that \code{partition_factor} names exactly one categorical/factor
#' variable with no missing values, then rebuilds it as a one-sided,
#' no-intercept formula (its canonical stored form).
#'
#' @param partition_factor A one-sided formula, or \code{NULL}.
#' @param obdata A data frame to validate \code{partition_factor} against.
#'
#' @return The normalized \code{partition_factor} formula, or \code{NULL}.
#'
#' @noRd
coerce_partition_factor <- function(partition_factor, obdata) {
  # coerce to factor
  if (!is.null(partition_factor)) {
    partition_factor_labels <- labels(terms(partition_factor))
    if (length(partition_factor_labels) > 1) {
      stop("Only one variable can be specified in partition_factor.", call. = FALSE)
    }
    check_no_na_columns(partition_factor_labels, obdata, "partition_factor")
    partition_mf <- model.frame(partition_factor, obdata)
    # ordered accepted alongside character/factor: partitioning only needs
    # group membership (category equality), never the ordinal spacing
    if (any(!attr(terms(partition_mf), "dataClasses") %in% c("character", "factor", "ordered"))) {
      stop("Partition factor variable must be categorical or factor.", call. = FALSE)
    }
    partition_factor <- reformulate(partition_factor_labels, intercept = FALSE)
    # partition_factor <- reformulate(paste0("as.character(", partition_factor_labels, ")"), intercept = FALSE)
  }
  partition_factor
}

#' Build the random-effect grouping/initial-value pieces shared across fitting paths
#'
#' Validates \code{random} has no missing values, builds its grouping
#' matrices/labels, and resolves \code{randcov_initial} -- either defaulting
#' it, or (if supplied) validating and canonicalizing its term names.
#'
#' @param random A one- or two-sided random effect formula, or \code{NULL}.
#' @param randcov_initial A random-effect variance initial-value object, or
#'   \code{NULL}.
#' @param obdata A data frame containing \code{random}'s grouping variables.
#' @param local_index A grouping vector (all-\code{1L} for a single, ungrouped
#'   fit) used to partition observations before building group matrices.
#'
#' @return A list with \code{randcov_initial}, \code{randcov_list}, \code{randcov_names},
#'   and \code{randcov_xlev} (all \code{NULL} if \code{random} is \code{NULL}).
#'
#' @noRd
build_randcov_list <- function(random, randcov_initial, obdata, local_index) {
  # store random effects list
  if (is.null(random)) {
    randcov_initial <- NULL
    randcov_list <- NULL
    randcov_names <- NULL
    randcov_xlev <- NULL
  } else {
    check_no_na_columns(all.vars(random), obdata, "random effect")
    randcov_names <- get_randcov_names(random)
    randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
    randcov_list <- get_randcov_list(local_index, randcov_Zs, randcov_names)
    randcov_xlev <- get_randcov_xlev(randcov_names, obdata)
    if (is.null(randcov_initial)) {
      randcov_initial <- randcov_initial()
    } else {
      randcov_given_names <- unlist(lapply(
        names(randcov_initial$initial),
        function(x) labels(terms(reformulate(x)))
      ))
      randcov_initial_names <- unique(unlist(lapply(randcov_given_names, get_randcov_name)))
      if (length(randcov_initial_names) != length(names(randcov_initial$initial))) {
        stop("No / can be specified in randcov_initial(). Please specify starting
             values for each variable (e.g., a/b = a + a:b)", call. = FALSE)
      }
      names(randcov_initial$initial) <- randcov_initial_names
      names(randcov_initial$is_known) <- randcov_initial_names
    }
  }
  list(randcov_initial = randcov_initial, randcov_list = randcov_list, randcov_names = randcov_names, randcov_xlev = randcov_xlev)
}

#' Build per-group partition-factor indicator matrices
#'
#' @param partition_factor A one-sided partition factor formula, or
#'   \code{NULL}.
#' @param obdata_list A list of (grouped) data frames.
#'
#' @return A list of partition indicator matrices, one per element of
#'   \code{obdata_list}, or \code{NULL} if \code{partition_factor} is
#'   \code{NULL}.
#'
#' @noRd
build_partition_list <- function(partition_factor, obdata_list) {
  # store partition matrix list
  if (!is.null(partition_factor)) {
    partition_list <- lapply(obdata_list, function(x) partition_matrix(partition_factor, x))
  } else {
    partition_list <- NULL
  }
  partition_list
}

#' Resolve range_constrain into a shared bound and a per-component decision
#'
#' Matches spmodel's \code{range_constrain}: when requested, each active range
#' parameter (tailup, taildown, euclid) is optimized on a bounded logit-odds
#' scale instead of an unconstrained log scale, capped at \code{4 *} the
#' bounding-box diagonal of the observed coordinates. Unlike spmodel (which
#' has a single range parameter), SSN2 has three range parameters spanning
#' two distance metrics (tailup/taildown act on hydrologic distance, euclid on
#' Euclidean distance); a single shared bound based on the (Euclidean)
#' bounding-box diagonal is used for all three, for simplicity and
#' consistency with the bounding box's use elsewhere in the package (e.g. the
#' big-data distance anchors in \code{get_data_object_bigdata()}).
#' Constraining is skipped per-parameter when that range is already known
#' (fixed), when its covariance type is \code{"none"}, or when its own
#' initial value already exceeds the bound (mirroring spmodel's own
#' auto-disable rules).
#'
#' @param obdata The observed data (used only for its coordinates).
#' @param initial_object A joint covariance initial-value object.
#' @param range_constrain A logical indicating whether constraining was
#'   requested.
#'
#' @return A list with \code{range_constrain_value} (the shared bound, or
#'   \code{NULL} if nothing ends up constrained) and a logical for each of
#'   \code{tailup_range_constrain}, \code{taildown_range_constrain}, and
#'   \code{euclid_range_constrain}.
#'
#' @noRd
get_range_constrain_setup <- function(obdata, initial_object, range_constrain) {
  if (!is.logical(range_constrain) || length(range_constrain) != 1 || is.na(range_constrain)) {
    stop("range_constrain must be TRUE or FALSE.", call. = FALSE)
  }

  bbox <- sf::st_bbox(obdata)
  bbox_dist <- sqrt((bbox[["xmax"]] - bbox[["xmin"]])^2 + (bbox[["ymax"]] - bbox[["ymin"]])^2)
  max_range_scale <- 4
  range_constrain_value <- max_range_scale * bbox_dist

  resolve_component <- function(component_initial, none_class) {
    if (!range_constrain || inherits(component_initial, none_class)) {
      return(FALSE)
    }
    # no initial range supplied (common case, later filled in by the
    # covariance starting-value grid search) -- nothing to compare, so
    # constrain freely
    if (!("range" %in% names(component_initial$initial))) {
      return(TRUE)
    }
    if (isTRUE(component_initial$is_known[["range"]])) {
      return(FALSE)
    }
    component_range <- component_initial$initial[["range"]]
    if (!is.na(component_range) && component_range > range_constrain_value) {
      return(FALSE)
    }
    TRUE
  }

  tailup_range_constrain <- resolve_component(initial_object$tailup_initial, "tailup_none")
  taildown_range_constrain <- resolve_component(initial_object$taildown_initial, "taildown_none")
  euclid_range_constrain <- resolve_component(initial_object$euclid_initial, "euclid_none")

  if (!any(tailup_range_constrain, taildown_range_constrain, euclid_range_constrain)) {
    range_constrain_value <- NULL
  }

  list(
    range_constrain_value = range_constrain_value,
    tailup_range_constrain = tailup_range_constrain,
    taildown_range_constrain = taildown_range_constrain,
    euclid_range_constrain = euclid_range_constrain
  )
}
