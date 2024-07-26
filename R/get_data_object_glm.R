#' Get the glm data object for use with many other functions.
#'
#' @param formula Model formula.
#' @param ssn.object SSN object.
#' @param additive Name of the additive function value column.
#' @param anisotropy Whether there is anisotropy.
#' @param initial_object The initial value object.
#' @param random Random effect formula.
#' @param randcov_initial The initial random effect object.
#' @param partition_factor Partition factor formula.
#' @param local Spatial indexing argument (not yet implemented)
#' @param ... Additional arguments
#'
#' @return The glm data object that contains various pieces of important information
#'   required for modeling.
#' @noRd
get_data_object_glm <- function(formula, ssn.object, family, additive, anisotropy,
                                initial_object, random, randcov_initial, partition_factor, local, ...) {
  sf_column_name <- attributes(ssn.object$obs)$sf_column
  crs <- attributes(ssn.object$obs[[sf_column_name]])$crs

  ## get response value in pid (data) order
  na_index <- is.na(sf::st_drop_geometry(ssn.object$obs)[[all.vars(formula)[1]]])
  ## get index in pid (data) order
  observed_index <- which(!na_index)
  missing_index <- which(na_index)

  # get ob data and frame objects
  obdata <- ssn.object$obs[observed_index, , drop = FALSE]
  # finding model frame
  obdata_model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.pass)
  # finding contrasts as ...
  dots <- list(...)
  if (!"contrasts" %in% names(dots)) dots$contrasts <- NULL

  # model matrix with potential NA
  X <- model.matrix(formula, obdata_model_frame, contrasts = dots$contrasts)
  # finding rows w/out NA
  ob_predictors <- complete.cases(X)
  if (any(!ob_predictors)) {
    stop("Cannot have NA values in predictors", call. = FALSE)
  }
  # subset obdata by nonNA predictors
  obdata <- obdata[ob_predictors, , drop = FALSE]

  # new model frame
  obdata_model_frame <- model.frame(formula, obdata, drop.unused.levels = TRUE, na.action = na.omit)
  # find terms
  terms_val <- terms(obdata_model_frame)
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
  # find response
  y_modr <- model.response(obdata_model_frame)
  if (NCOL(y_modr) == 2) {
    y <- y_modr[, 1, drop = FALSE]
    size <- rowSums(y_modr)
  } else {
    if (family == "binomial") {
      if (is.factor(y_modr)) {
        if (length(levels(y_modr)) != 2) {
          stop("When family is binomial, a factor response must have exactly two levels.", call. = FALSE)
        }
        y_modr <- ifelse(y_modr == levels(y_modr)[1], 0, 1)
      }
      if (is.logical(y_modr)) {
        y_modr <- ifelse(y_modr, 1, 0) # or as.numeric()
      }
      size <- rep(1, n)
    } else {
      size <- NULL
    }
    y <- as.matrix(y_modr, ncol = 1)
  }

  # handle offset
  offset <- model.offset(obdata_model_frame)
  if (!is.null(offset)) {
    offset <- as.matrix(offset, ncol = 1)
  }

  # see if response is numeric
  if (!is.numeric(y)) {
    stop("Response variable must be numeric", call. = FALSE)
  }

  # error if no variance
  if (var(y) == 0) {
    stop("The response has no variability. Model fit unreliable.", call. = FALSE)
  }

  # checks on y
  response_checks_glm(family, y, size)

  # error if p >= n
  if (p >= n) {
    stop("The number of fixed effects is at least as large as the number of observations (p >= n). Consider reducing the number of fixed effects and rerunning ssn_lm().", call. = FALSE)
  }

  # find s2 for initial values
  y_trans <- log(y + 1)
  qr_val <- qr(X)
  R_val <- qr.R(qr_val)
  betahat <- backsolve(R_val, qr.qty(qr_val, y_trans))
  resid <- y_trans - X %*% betahat
  s2 <- sum(resid^2) / (n - p)
  # diagtol <- 1e-4
  diagtol <- min(1e-4, 1e-4 * s2)

  # correct anisotropy
  anisotropy <- get_anisotropy_corrected(anisotropy, initial_object)

  # coerce to factor
  if (!is.null(partition_factor)) {
    partition_factor_labels <- labels(terms(partition_factor))
    if (length(partition_factor_labels) > 1) {
      stop("Only one variable can be specified in partition_factor.", call. = FALSE)
    }
    partition_mf <- model.frame(partition_factor, obdata)
    if (any(!attr(terms(partition_mf), "dataClasses") %in% c("character", "factor"))) {
      stop("Partition factor variable must be categorical or factor.", call. = FALSE)
    }
    partition_factor <- reformulate(partition_factor_labels, intercept = FALSE)
    # partition_factor <- reformulate(paste0("as.character(", partition_factor_labels, ")"), intercept = FALSE)
  }

  # find index
  # if (is.null(local)) {
  #   if (n > 3000) {
  #     local <- TRUE
  #     message("Because the sample size exceeds 3000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun splm() with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
  #   } else {
  #     local <- FALSE
  #   }
  # }

  local <- list(index = rep(1, n))
  local <- get_local_list_estimation(local, obdata, n, partition_factor)

  # store data list
  obdata_list <- split.data.frame(obdata, local$index)

  # store X and y
  X_list <- split.data.frame(X, local$index)
  y_list <- split.data.frame(y, local$index)
  ones_list <- lapply(obdata_list, function(x) matrix(rep(1, nrow(x)), ncol = 1))
  if (!is.null(size)) {
    size_list <- split(size, local$index) # just split because vector not matrix
    size <- as.vector(do.call("c", size_list)) # rearranging size by y list
  }

  # organize offset (as a one col matrix)
  if (!is.null(offset)) {
    offset <- do.call("rbind", (split.data.frame(offset, local$index)))
  }

  # store random effects list
  if (is.null(random)) {
    randcov_initial <- NULL
    randcov_list <- NULL
    randcov_names <- NULL
  } else {
    randcov_names <- get_randcov_names(random)
    randcov_Zs <- get_randcov_Zs(obdata, randcov_names)
    randcov_list <- get_randcov_list(local$index, randcov_Zs, randcov_names)
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

  # store partition matrix list
  if (!is.null(local$partition_factor)) {
    partition_list <- lapply(obdata_list, function(x) partition_matrix(local$partition_factor, x))
  } else {
    partition_list <- NULL
  }

  # find dist object
  dist_object <- get_dist_object(ssn.object, initial_object, additive, anisotropy)

  # find maxes
  tailup_none <- inherits(initial_object$tailup_initial, "tailup_none")
  taildown_none <- inherits(initial_object$taildown_initial, "taildown_none")
  if (tailup_none && taildown_none) {
    tail_max <- Inf
  } else {
    tail_max <- max(dist_object$hydro_mat * dist_object$mask_mat)
  }

  euclid_none <- inherits(initial_object$euclid_initial, "euclid_none")
  if (euclid_none) {
    euclid_max <- Inf
  } else {
    if (anisotropy) {
      euclid_max <- max(as.matrix(dist(cbind(dist_object$.xcoord, dist_object$.ycoord))))
    } else {
      euclid_max <- max(dist_object$euclid_mat) # no anisotropy
    }
  }

  # find dist observed object
  dist_object <- get_dist_object_oblist(dist_object, observed_index, local$index)
  # rename as oblist to not store two sets
  dist_object_oblist <- dist_object

  # store order
  order <- unlist(split(seq_len(n), local$index), use.names = FALSE)

  # store global pid
  pid <- ssn_get_netgeom(ssn.object$obs, "pid")$pid

  # restructure ssn
  ssn.object <- restruct_ssn_missing(ssn.object, observed_index, missing_index)

  list(
    anisotropy = anisotropy,
    additive = additive,
    contrasts = dots$contrasts,
    crs = crs,
    diagtol = diagtol,
    # dist_object = dist_object,
    dist_object_oblist = dist_object_oblist,
    euclid_max = euclid_max,
    family = family,
    formula = formula,
    local_index = local$index,
    missing_index = missing_index,
    n = n,
    ncores = local$ncores,
    observed_index = observed_index,
    offset = offset,
    ones_list = ones_list,
    order = order,
    p = p,
    parallel = local$parallel,
    partition_factor_initial = partition_factor,
    partition_factor = local$partition_factor,
    partition_list = partition_list,
    pid = pid,
    randcov_initial = randcov_initial,
    randcov_list = randcov_list,
    randcov_names = randcov_names,
    sf_column_name = sf_column_name,
    size = size,
    ssn.object = ssn.object,
    s2 = s2,
    tail_max = tail_max,
    terms = terms_val,
    var_adjust = local$var_adjust,
    X_list = X_list,
    xlevels = xlevels,
    y_list = y_list
  )
}
