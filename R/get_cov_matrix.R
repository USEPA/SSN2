get_cov_matrix <- function(params_object, dist_object_oblist,
                           randcov_list = NULL, partition_list = NULL, anisotropy,
                           de_scale, diagtol) {
  # THIS IS ELEMENT BY ELEMENT

  # could use Reduce if arguments the same

  # compute spatial covariance matrix
  cov_matrix <- cov_matrix(params_object$tailup, dist_object_oblist) + # tailup covariance
    cov_matrix(params_object$taildown, dist_object_oblist) + # taildown covariance
    cov_matrix(params_object$euclid, dist_object_oblist, anisotropy) + # euclid covariance
    cov_matrix(params_object$nugget, dist_object_oblist, de_scale, diagtol) # nugget covariance

  # add random effects if necessary
  if (!is.null(randcov_list)) {
    cov_matrix <- cov_matrix + randcov_matrix(params_object$randcov, randcov_list)
  }
  # apply partition factor if necessary
  if (!is.null(partition_list)) {
    cov_matrix <- cov_matrix * partition_list
  }

  # return matrix
  cov_matrix
}

get_cov_matrix_list <- function(params_object, data_object) {
  # this is a list of elements

  de_scale <- sum(params_object$tailup[["de"]], params_object$taildown[["de"]], params_object$euclid[["de"]])

  if (is.null(params_object$randcov) & is.null(data_object$partition_list)) {
    # no random effects or partition factor
    cov_matrix_list <- mapply(
      d = data_object$dist_object_oblist,
      function(d) {
        get_cov_matrix(params_object, d,
          anisotropy = data_object$anisotropy, de_scale = de_scale, diagtol = data_object$diagtol
        )
      },
      SIMPLIFY = FALSE
    )
  } else if (!is.null(params_object$randcov) & is.null(data_object$partition_list)) {
    # random effects no partition factor
    cov_matrix_list <- mapply(
      d = data_object$dist_object_oblist,
      r = data_object$randcov_list,
      function(d, r) {
        get_cov_matrix(params_object, d, r,
          anisotropy = data_object$anisotropy, de_scale = de_scale, diagtol = data_object$diagtol
        )
      },
      SIMPLIFY = FALSE
    )
  } else if (is.null(params_object$randcov) & !is.null(data_object$partition_list)) {
    # no random effects partition factor
    cov_matrix_list <- mapply(
      d = data_object$dist_object_oblist,
      p = data_object$partition_list,
      function(d, p) {
        get_cov_matrix(params_object, d,
          partition_list = p, anisotropy = data_object$anisotropy, de_scale = de_scale, diagtol = data_object$diagtol
        )
      },
      SIMPLIFY = FALSE
    )
  } else {
    # random effects and partition factor
    cov_matrix_list <- mapply(
      d = data_object$dist_object_oblist,
      r = data_object$randcov_list,
      p = data_object$partition_list,
      function(d, r, p) {
        get_cov_matrix(params_object, d, r, p,
          anisotropy = data_object$anisotropy, de_scale = de_scale, diagtol = data_object$diagtol
        )
      },
      SIMPLIFY = FALSE
    )
  }
  cov_matrix_list
}
