get_data_object <- function(formula, ssn.object, additive, anisotropy,
                            initial_object, random, randcov_initial, partition_factor, local,
                            range_constrain = FALSE, ...) {
  sf_column_name <- attributes(ssn.object$obs)$sf_column
  crs <- attributes(ssn.object$obs[[sf_column_name]])$crs

  ## get response value in pid (data) order
  na_index <- is.na(sf::st_drop_geometry(ssn.object$obs)[[all.vars(formula)[1]]])
  ## get index in pid (data) order
  observed_index <- which(!na_index)
  missing_index <- which(na_index)
  # observed_index <- !na_index
  # missing_index <- na_index

  # get ob data and frame objects
  obdata <- ssn.object$obs[observed_index, , drop = FALSE]
  mm <- get_model_matrix_object(formula, obdata, sf_column_name = sf_column_name, ...)
  obdata <- mm$obdata
  formula <- mm$formula
  obdata_model_frame <- mm$obdata_model_frame
  terms_val <- mm$terms_val
  X <- mm$X
  dots <- mm$dots
  xlevels <- mm$xlevels
  p <- mm$p
  n <- mm$n
  # find response
  y <- as.matrix(model.response(obdata_model_frame), ncol = 1)
  # adjust response to reflect offset
  offset <- model.offset(obdata_model_frame)
  if (!is.null(offset)) {
    offset <- as.matrix(offset, ncol = 1)
    y <- y - offset
  }

  check_response_numeric_and_variable(y)

  check_p_lt_n(p, n)

  # find s2 for initial values
  qr_val <- qr(X)
  R_val <- qr.R(qr_val)
  betahat <- backsolve(R_val, qr.qty(qr_val, y))
  resid <- y - X %*% betahat
  s2 <- sum(resid^2) / (n - p)
  diagtol <- 0

  # correct anisotropy
  anisotropy <- get_anisotropy_corrected(anisotropy, initial_object)

  partition_factor <- coerce_partition_factor(partition_factor, obdata)

  local <- list(index = rep(1, n))
  local <- get_local_list_estimation(local, obdata, n, partition_factor)

  # store data list
  obdata_list <- split.data.frame(obdata, local$index)

  # store X and y
  X_list <- split.data.frame(X, local$index)
  y_list <- split.data.frame(y, local$index)
  ones_list <- lapply(obdata_list, function(x) matrix(rep(1, nrow(x)), ncol = 1))

  # organize offset (as a one col matrix)
  if (!is.null(offset)) {
    offset <- do.call("rbind", (split.data.frame(offset, local$index)))
  }

  rc <- build_randcov_list(random, randcov_initial, obdata, local$index)
  randcov_initial <- rc$randcov_initial
  randcov_list <- rc$randcov_list
  randcov_names <- rc$randcov_names
  randcov_xlev <- rc$randcov_xlev

  partition_list <- build_partition_list(local$partition_factor, obdata_list)
  partition_xlev <- get_partition_xlev(local$partition_factor, obdata)

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

  range_setup <- get_range_constrain_setup(obdata, initial_object, range_constrain)
  range_constrain_value <- range_setup$range_constrain_value
  tailup_range_constrain <- range_setup$tailup_range_constrain
  taildown_range_constrain <- range_setup$taildown_range_constrain
  euclid_range_constrain <- range_setup$euclid_range_constrain

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
    formula = formula,
    local_index = local$index,
    missing_index = missing_index,
    n = n,
    ncores = local$ncores,
    observed_index = observed_index,
    obdata_list = NULL,
    offset = offset,
    ones_list = ones_list,
    order = order,
    order_bigdata = NULL,
    p = p,
    parallel = local$parallel,
    partition_factor_initial = partition_factor,
    partition_factor = local$partition_factor,
    partition_list = partition_list,
    pid = pid,
    randcov_initial = randcov_initial,
    randcov_list = randcov_list,
    randcov_names = randcov_names,
    randcov_xlev = randcov_xlev,
    partition_xlev = partition_xlev,
    sf_column_name = sf_column_name,
    ssn.object = ssn.object,
    s2 = s2,
    tail_max = tail_max,
    range_constrain_value = range_constrain_value,
    tailup_range_constrain = tailup_range_constrain,
    taildown_range_constrain = taildown_range_constrain,
    euclid_range_constrain = euclid_range_constrain,
    terms = terms_val,
    var_adjust = local$var_adjust,
    X_list = X_list,
    xlevels = xlevels,
    y_list = y_list
  )
}
