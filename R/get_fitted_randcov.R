get_fitted_randcov <- function(data_object, randcov_params, SigInv_r_list) {
  if (is.null(names(randcov_params))) return(NULL)
  partition <- if ("partition_factor_initial" %in% names(data_object)) {
    data_object$partition_factor_initial
  } else {
    data_object$partition_factor
  }
  if (is.null(partition)) {
    out <- lapply(names(randcov_params), function(term) {
      Z <- data_object$randcov_list[[1L]][[term]]$Z
      total <- numeric(NCOL(Z))
      count <- integer(NCOL(Z))
      for (k in seq_along(SigInv_r_list)) {
        value <- randcov_params[[term]] * as.numeric(crossprod(
          data_object$randcov_list[[k]][[term]]$Z, SigInv_r_list[[k]]))
        total <- total + value
        count <- count + (value != 0)
      }
      value <- stats::setNames(total / pmax(count, 1L), colnames(Z))
      value[order(names(value))]
    })
  } else {
    observed <- data_object$ssn.object$obs
    # Local fitting reorders rows before building its precision-residual blocks.
    if (!is.null(data_object$order_bigdata)) {
      observed <- observed[data_object$order_bigdata, , drop = FALSE]
    }
    design <- get_randcov_design(observed, randcov_params, partition)
    index <- if (is.null(data_object$local_index)) list(seq_len(NROW(observed))) else
      unname(split(seq_len(NROW(observed)), data_object$local_index))
    total <- numeric(length(design$variance))
    count <- integer(length(total))
    for (k in seq_along(index)) {
      value <- design$variance * as.numeric(crossprod(
        design$observed[index[[k]], , drop = FALSE], SigInv_r_list[[k]]))
      total <- total + value
      count <- count + (value != 0)
    }
    value <- total / pmax(count, 1L)
    out <- lapply(names(randcov_params), function(term) {
      rows <- which(design$info$term == term)
      stats::setNames(value[rows], make.unique(design$info$level[rows]))
    })
  }
  stats::setNames(out, names(randcov_params))
}
