#' Specify spatial indexing covariance matrix adjustment. Not currently
#'   relevant as spatial indexing has not yet been implemented.
#'
#' @param invcov_betahat_list Placeholder.
#' @param betahat_list Placeholder.
#' @param betahat Placeholder.
#' @param eigenprods_list Placeholder.
#' @param data_object Placeholder.
#' @param params_object Placeholder.
#' @param cov_betahat_noadjust Placeholder.
#' @param var_adjust Placeholder.
#'
#' @noRd
cov_betahat_adjust <- function(invcov_betahat_list, betahat_list,
                               betahat, eigenprods_list, data_object, params_object,
                               cov_betahat_noadjust, var_adjust) {

  P <- length(betahat_list)
  # reset var_adjust if only one partition
  # var_adjust == "theoretical" ||
  if (P == 1 || (inherits(params_object$tailup, "tailup_none") && inherits(params_object$taildown, "taildown_none") && inherits(params_object$euclid, "euclid_none") && inherits(params_object$nugget, "tailup_nugget"))) {
    var_adjust <- "none"
  }

  if (var_adjust == "empirical") {
    cov_betahat_adjust_list <- lapply(betahat_list, function(x) {
      betahat_diff <- x - betahat
      tcrossprod(betahat_diff, betahat_diff)
    })
    cov_betahat_adjust_val <- Reduce("+", cov_betahat_adjust_list) / (P * (P - 1))
  } else if (var_adjust == "pooled") {
    cov_betahat_list <- tryCatch(
      error = function(cnd) NULL,
      lapply(invcov_betahat_list, function(x) chol2inv(chol(forceSymmetric(x))))
    )
    if (is.null(cov_betahat_list)) {
      var_adjust <- "none"
      warning(
        "At least one partition's inverse covariance matrix is singular. Readjusting using var_adjust = \"none\".",
        call. = FALSE
      )
    } else {
      cov_betahat_adjust_val <- Reduce("+", cov_betahat_list) / P^2
    }
  }

  if (var_adjust == "theoretical") {
    index_grid <- expand.grid(d1 = seq_len(P), d2 = seq_len(P))
    index_grid <- index_grid[index_grid$d2 > index_grid$d1, , drop = FALSE]
    index_list <- split(index_grid, seq_len(NROW(index_grid)))
    if (data_object$parallel) {
      W_adjust_list <- parallel::parLapply(data_object$cl, index_list, get_W_ij_parallel,
                                           eigenprods_list = eigenprods_list,
                                           params_object = params_object, randcov_params = randcov_params,
                                           randcov_names = names(randcov_params), data_object = data_object
      )
    } else {
      W_adjust_list <- lapply(index_list, function(x) {
        get_W_ij(
          x$d1, x$d2, eigenprods_list,
          params_object, randcov_params,
          names(randcov_params), data_object
        )
      })
    }
    W_adjust <- Reduce("+", W_adjust_list)
    cov_betahat_adjust_val <- cov_betahat_noadjust + cov_betahat_noadjust %*% W_adjust %*% cov_betahat_noadjust
  } else if (var_adjust == "none") {
    cov_betahat_adjust_val <- cov_betahat_noadjust
  }

  cov_betahat_adjust_val
}

get_W_ij <- function(d1_index, d2_index, eigenprods_list, params_object,
                     randcov_params, randcov_names, data_object) {


  d1 <- data_object$obdata_list[[d1_index]]
  d2 <- data_object$obdata_list[[d2_index]]


  dist_object_bigdata_cross <- get_dist_object_bigdata_cross(d1, d2, params_object, data_object)


  # random effects
  if (!is.null(data_object$randcov_list)) {
    Zs_cross <- lapply(randcov_names, function(x) {
      list(
        ZZt = Matrix::tcrossprod(
          data_object$randcov_list[[d1_index]][[x]]$Z,
          data_object$randcov_list[[d2_index]][[x]]$Z
        )
      )
    })
    names(Zs_cross) <- randcov_names
  } else {
    randcov_params <- NULL
    Zs_cross <- NULL
  }


  # partition matrix
  if (!is.null(data_object$partition_factor)) {
    # finding the formula
    partition_formula <- reformulate(labels(terms(data_object$partition_factor)), intercept = FALSE)
    #
    d1_partition <- Matrix::Matrix(model.matrix(partition_formula, d1), sparse = TRUE)
    d2_partition <- Matrix::Matrix(model.matrix(partition_formula, d2), sparse = TRUE)
    partition_matrix_cross_val <- Matrix::tcrossprod(d1_partition, d2_partition)
  } else {
    partition_matrix_cross_val <- NULL
  }

  cov_d1d2_cross <- get_cov_matrix_cross(params_object, dist_object_bigdata_cross, randcov_params, Zs_cross, partition_matrix_cross_val, data_object)

  part1 <- eigenprods_list[[d1_index]]$SigInv_X
  part2 <- eigenprods_list[[d2_index]]$SigInv_X
  W_ij_half <- Matrix::crossprod(part1, cov_d1d2_cross) %*% part2
  W_ij <- W_ij_half + t(W_ij_half)
}


get_W_ij_parallel <- function(index_list, eigenprods_list, params_object,
                              randcov_params, randcov_names, data_object) {
  d1_index <- index_list$d1
  d2_index <- index_list$d2
  get_W_ij(
    d1_index, d2_index, eigenprods_list,
    params_object, randcov_params,
    randcov_params, data_object
  )
}
