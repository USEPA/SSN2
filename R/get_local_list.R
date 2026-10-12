get_local_list_conditional_lowrank <- function(local, object, newdata, n, n_pred) {
  names_local <- names(local)

  if (!"method_base" %in% names_local) local$method_base <- "base"
  if (length(local$method_base) != 1L || !local$method_base %in% c("base", "all")) {
    stop("local$method_base must be \"base\" or \"all\".", call. = FALSE)
  }
  if (!"method_new" %in% names_local) local$method_new <- "base"
  if (length(local$method_new) != 1L || !local$method_new %in% c("base", "all")) {
    stop("local$method_new must be \"base\" or \"all\".", call. = FALSE)
  }
  if (!"size_base" %in% names_local) local$size_base <- 5000L
  if (!"size_new" %in% names_local) local$size_new <- 1000L
  if (!"reorder_base" %in% names_local) local$reorder_base <- "grts"
  if (length(local$reorder_base) != 1L || !local$reorder_base %in% c("none", "random", "grts")) {
    stop("local$reorder_base must be \"none\", \"random\", or \"grts\".", call. = FALSE)
  }
  if (!"reorder_new" %in% names_local) local$reorder_new <- "random"
  if (length(local$reorder_new) != 1L || !local$reorder_new %in% c("none", "random")) {
    stop("local$reorder_new must be \"none\" or \"random\".", call. = FALSE)
  }
  if (!"kmeans_new" %in% names_local) {
    local$kmeans_new <- !identical(local$reorder_new, "none")
  }
  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$size_base >= n) local$method_base <- "all"
  if (local$size_new >= n_pred) local$method_new <- "all"

  index_base <- seq_len(n)
  index_new <- seq_len(n_pred)

  if (!identical(local$method_base, "all")) {
    observed <- object$ssn.object$obs
    rows_base <- get_decorrelate_rows(observed)
    index_base <- get_decorrelate_order(rows_base, local$reorder_base, observed)
    index_base <- index_base[seq_len(local$size_base)]
  }

  if (!identical(local$method_new, "all")) {
    rows_new <- get_decorrelate_rows(newdata)
    index_new <- get_decorrelate_order(rows_new, local$reorder_new, newdata)
    groups <- ceiling(n_pred / local$size_new)
    index_new <- get_local_lowrank_blocks(index_new, newdata, groups, local$kmeans_new)
  }

  local$index <- list(base = index_base, new = index_new)

  if (local$parallel) {
    n_index <- if (identical(local$method_new, "all")) 1L else length(local$index$new)
    local$ncores <- if ("ncores" %in% names_local) {
      min(n_index, local$ncores, parallel::detectCores())
    } else {
      min(n_index, parallel::detectCores())
    }
  }

  local
}

get_local_list_conditional <- function(local, object, newdata = NULL) {
  if (is.null(newdata)) newdata <- object$ssn.object$obs[FALSE, , drop = FALSE]
  random_options <- c("size_random", "groups_random", "method_random")
  random_settings <- if (is.list(local)) Filter(Negate(is.null), local[intersect(names(local), random_options)]) else list()
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
  if (is.null(local)) {
    if (object$n > 5000 || NROW(newdata) > 5000) {
      local <- TRUE
      message("Because the data size or number of conditional simulations exceeds 5,000, we are setting local = TRUE to perform computationally efficient approximations. To override this behavior and compute the exact solution, rerun with local = FALSE. Be aware that setting local = FALSE may result in exceedingly long computational times.")
    } else local <- FALSE
  }
  exact <- identical(local, FALSE)
  if (identical(local, TRUE)) local <- list()
  if (exact) local <- list(approximation = "low-rank", method_base = "all", method_new = "all")
  if (!is.list(local)) stop("local must be FALSE, TRUE, or a list.", call. = FALSE)
  if (is.null(local[["approximation"]])) local$approximation <- "vecchia"
  approximation <- local[["approximation"]]
  if (length(approximation) != 1L || is.na(approximation) || !approximation %in% c("low-rank", "vecchia")) {
    stop('local$approximation must be "low-rank" or "vecchia".', call. = FALSE)
  }
  if (length(random_settings) && (approximation != "vecchia" || is.null(object$random) || identical(local[["method"]], "all"))) {
    stop("local$size_random, local$groups_random, and local$method_random are only used for Vecchia simulation with random effects.", call. = FALSE)
  }
  if (approximation == "vecchia") {
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
    if (local[["method"]] == "all") return(get_local_list_conditional(FALSE, object, newdata))
  } else {
    local <- get_local_list_conditional_lowrank(local, object, newdata, object$n, NROW(newdata))
  }
  if (is.null(local[["parallel"]])) local$parallel <- FALSE
  local$exact <- exact || (approximation == "low-rank" &&
    local[["method_base"]] == "all" && local[["method_new"]] == "all")
  if (!local$exact) local$chunk_size <- get_conditional_chunk_size(local)
  local
}

get_local_list_conditional_vecchia <- function(local, object, newdata) {
  lowrank_options <- c("method_base", "size_base", "reorder_base", "method_new", "size_new", "reorder_new", "kmeans_new")
  if (any(names(local) %in% lowrank_options)) {
    stop('Base and block settings require local$approximation = "low-rank".', call. = FALSE)
  }
  if (is.null(local[["size"]])) local$size <- if (inherits(object, "ssn_glm")) 60L else 30L
  if ("size_response" %in% names(local)) stop("local$size_response is not supported; local$size controls both observed neighborhoods.", call. = FALSE)
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
  if (isTRUE(local[["parallel"]])) {
    warning('local$parallel is not used when local$approximation = "vecchia". Ignoring.', call. = FALSE)
  }
  local$parallel <- FALSE
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
