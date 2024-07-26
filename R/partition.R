#' Create a partition matrix
#'
#' @param partition_factor A partition factor (formula)
#' @param data Data
#'
#' @return A partition matrix
#'
#' @noRd
partition_matrix <- function(partition_factor = NULL, data) {
  if (is.null(partition_factor)) {
    partition_matrix_val <- NULL
  } else {
    # finding the formula
    partition_formula <- reformulate(labels(terms(partition_factor)), intercept = FALSE)
    # use regular contrasts here so matrix all zeros and ones
    partition_model_frame <- model.frame(partition_formula, data)
    if (length(unique(as.character(unlist(partition_model_frame)))) == 1) {
      partition_model_val <- Matrix::Matrix(matrix(1, nrow = NROW(partition_model_frame), ncol = 1), sparse = TRUE)
      # partition_model_val <- Matrix::Matrix(as.matrix(partition_model_frame), sparse = TRUE)
    } else {
      partition_model_val <- Matrix::Matrix(model.matrix(partition_formula, data), sparse = TRUE)
    }
    partition_matrix_val <- tcrossprod(partition_model_val, partition_model_val)
  }
  partition_matrix_val
}



#' Get names of partition factors
#'
#' @param partition_factor A partition factor
#'
#' @noRd
get_partition_names <- function(partition_factor) {
  # could use old version for partition names
  # get the partition formula and turn it into a named vector here
  labels_initial <- labels(terms(partition_factor))
  new_labels <- unlist(lapply(labels_initial, get_partition_name))
}

#' Get adjusted name of a partition factor
#'
#' @param label A partition factor name
#'
#' @noRd
get_partition_name <- function(label) {
  if (grepl("|", label, fixed = TRUE)) {
    new_label <- label
  } else {
    new_label <- paste("1", label, sep = " | ")
  }
  new_label
}
