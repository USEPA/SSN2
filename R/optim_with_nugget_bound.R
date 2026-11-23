#' Call optim with an optional nugget lower bound
#'
#' @noRd
optim_with_nugget_bound <- function(optim_args, optim_dotlist, lower = NULL) {
  if (is.null(lower)) {
    return(do.call("optim", c(optim_args, optim_dotlist)))
  }
  parameter <- names(optim_args$par)
  upper <- optim_dotlist$upper
  if (!identical(parameter, "nugget_log") ||
    length(lower) != 1L || !is.finite(lower) || length(upper) != 1L || !is.finite(upper) || lower >= upper) {
    stop("Invalid bounds for the nugget variance search.", call. = FALSE)
  }
  optim_dotlist$lower <- lower
  do.call("optim", c(optim_args, optim_dotlist))
}

#' Lower search bound for a sole unknown nugget variance
#'
#' @noRd
nugget_only_lower <- function(par, initial_object, diagtol) {
  if (!identical(names(par), "nugget_log") || inherits(initial_object$nugget_initial, "nugget_none")) {
    return(NULL)
  }
  params <- get_params_object_known(initial_object)
  params$nugget[["nugget"]] <- 0
  params <- floor_estimated_nugget(params, list(nugget = FALSE), diagtol)
  floor <- params$nugget[["nugget"]]
  if (isTRUE(floor == 0)) return(NULL)
  log(floor)
}
