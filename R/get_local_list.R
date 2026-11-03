validate_simulation_local <- function(local, options) {
  if (is.logical(local) && length(local) == 1L && !is.na(local)) return(as.logical(local))
  if (!is.list(local)) stop("local must be FALSE, TRUE, or a list.", call. = FALSE)
  if (length(local)) {
    keys <- names(local)
    if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
      stop("local options must have nonempty, unique, nonmissing names.", call. = FALSE)
    }
    unknown <- setdiff(keys, options)
    if (length(unknown)) {
      stop("Unknown local option: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
  }
  local
}

get_local_list_conditional <- function(local, object, newdata = NULL) {
  if (is.null(newdata)) newdata <- object$ssn.object$obs[FALSE, , drop = FALSE]
  if (is.null(local)) {
    if (object$n > 5000 || NROW(newdata) > 5000) {
      local <- TRUE
      message("Because the data size or number of conditional simulations exceeds 5,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else local <- FALSE
  }
  random_options <- c("size_random", "groups_random", "method_random")
  local <- validate_simulation_local(local, c("method", "size", "ordering", "chunk_size", random_options))
  if (identical(local, FALSE)) return(list(exact = TRUE, method = "all"))
  if (identical(local, TRUE)) local <- list()
  random_settings <- Filter(Negate(is.null), local[intersect(names(local), random_options)])
  for (name in intersect(names(random_settings), c("size_random", "groups_random"))) {
    value <- random_settings[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value)) {
      stop(paste0("local$", name, " must be a positive integer."), call. = FALSE)
    }
  }
  if (!is.null(random_settings[["groups_random"]]) && random_settings[["groups_random"]] > object$n) {
    stop("local$groups_random cannot exceed the number of observed sites.", call. = FALSE)
  }
  method_random <- random_settings[["method_random"]]
  if (!is.null(method_random) && (!is.character(method_random) || length(method_random) != 1L ||
      is.na(method_random) || !method_random %in% c("kmeans", "random"))) {
    stop('local$method_random must be "kmeans" or "random".', call. = FALSE)
  }
  if (length(random_settings) && (is.null(object$random) || identical(local[["method"]], "all"))) {
    stop("local$size_random, local$groups_random, and local$method_random are only used for Vecchia simulation with random effects.", call. = FALSE)
  }
  if (!is.null(object$random)) {
    if (!is.null(local[["size_random"]])) {
      local$groups_random <- ceiling(object$n / local[["size_random"]])
    } else if (!is.null(local[["groups_random"]])) {
      local$size_random <- ceiling(object$n / local[["groups_random"]])
    } else {
      local$size_random <- 1000L
      local$groups_random <- ceiling(object$n / 1000L)
    }
    if (is.null(local[["method_random"]])) local$method_random <- "kmeans"
  }
  local <- get_local_list_conditional_vecchia(local, object, newdata)
  if (local[["method"]] == "all") return(list(exact = TRUE, method = "all"))
  local$exact <- FALSE
  local$chunk_size <- get_conditional_chunk_size(local)
  local
}

get_local_list_conditional_vecchia <- function(local, object, newdata) {
  if (is.null(local[["size"]])) local$size <- if (inherits(object, "ssn_glm")) 60L else 30L
  if (is.null(local[["method"]])) local$method <- "covariance"
  method <- local[["method"]]
  if (length(method) != 1L || is.na(method) || !method %in% c("all", "covariance")) {
    stop('local$method must be "all" or "covariance".', call. = FALSE)
  }
  size <- local[["size"]]
  if (method != "all" && (!is.numeric(size) || length(size) != 1L || is.na(size) ||
      !is.finite(size) || size < 1 || size != floor(size))) {
    stop("local$size must be a positive integer.", call. = FALSE)
  }
  local$ordering <- get_decorrelate_ordering(local[["ordering"]])
  if (method != "all" && !is.null(object$random)) {
    local$random_blocks <- conditional_random_blocks(object, local[["size_random"]],
      local[["groups_random"]], local[["method_random"]])
  }
  if (method != "all" && inherits(object, "ssn_glm")) local$ordering_observed <- local$ordering
  local$order <- if (method == "all" || NROW(newdata) <= 1L) seq_len(NROW(newdata)) else
    get_decorrelate_order(get_decorrelate_rows(newdata), local$ordering, newdata)
  local$inv_order <- order(local$order)
  local
}
