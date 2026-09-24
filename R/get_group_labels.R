#' Row group labels derived from a one-sided formula's model.matrix()
#'
#' @param reform A one-sided formula for the grouping variable(s).
#' @param data The data.
#' @param xlev Optional factor levels to enforce (from \code{.getXlevels()}),
#'   so this call spans the same levels as a different (e.g. full) data set.
#' @param na_pass If \code{TRUE}, use \code{na.action = na.pass} (for
#'   newdata-side calls, so a row is not dropped/rejected merely because a
#'   different row has a missing value elsewhere); if \code{FALSE}, use the
#'   default \code{na.action}. With \code{na_pass = TRUE}, a missing grouping
#'   value errors because it does not identify exactly one indicator column.
#'
#' @details \code{model.matrix()} one-hot encodes the (possibly
#'   multi-variable) grouping formula into dummy columns named with its own
#'   "varname + level" convention; for each row, \code{which()} finds the
#'   single column that is 1, and its column name becomes that row's group
#'   label. This is the shared implementation behind random-effect grouping
#'   labels and the partition factor's own labels, at both
#'   fitting/observed-data time and newdata/prediction time. Expects a single
#'   combined grouping term (e.g. an interaction like \code{g:h}), not an
#'   additive multi-term formula (\code{g + h}), since only a single combined
#'   term one-hot encodes to exactly one column per row.
#'
#' @return A character vector, one group label per row of \code{data}.
#'
#' @noRd
model_matrix_group_labels <- function(reform, data, xlev = NULL, na_pass = FALSE) {
  mf <- if (na_pass) {
    model.frame(reform, data, na.action = na.pass, xlev = xlev)
  } else if (is.null(xlev)) {
    model.frame(reform, data)
  } else {
    model.frame(reform, data, xlev = xlev)
  }
  mx <- model.matrix(reform, mf)
  names_mx <- colnames(mx)
  split_rows <- split(mx, seq_len(NROW(mx)))
  names_mx[vapply(split_rows, function(y) which(as.logical(y)), numeric(1))]
}
