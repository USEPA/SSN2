get_local_list_estimation <- function(local, data, n, partition_factor) {

  if (is.logical(local)) {
    if (local) {
      local <- list()
    }
  }

  names_local <- names(local)

  # errors
  if (!"index" %in% names_local && "method" %in% names_local) {
    if (!local$method %in% c("random", "kmeans")) {
      stop("Invalid local method. Local method must be \"random\" or \"kmeans\".", call. = FALSE)
    }
  }

  if (!"index" %in% names_local && "var_adjust" %in% names_local) {
    if (!local$var_adjust %in% c("none", "theoretical", "empirical", "pooled")) {
      stop("Invalid local var_adjust. Local var_adjust must be \"none\", \"theoretical\", \"empirical\", or \"pooled\".", call. = FALSE)
    }
  }

  if ("index" %in% names_local) {
    # if index is a factor and there are levels in the factor not in the observed
    # data, the code will fail. Storing as character prevents this (acts as droplevels)
    if (is.factor(local$index)) {
      local$index <- as.character(local$index)
    }
    local$size <- NULL
    local$groups <- NULL
    local$method <- NULL
  } else {
    if (!"size" %in% names_local) {
      if ("groups" %in% names_local) {
        local$size <- ceiling(n / local$groups)
      } else {
        local$size <- 200
        local$groups <- ceiling(n / local$size)
      }
    } else {
      local$groups <- ceiling(n / local$size)
    }
    if (!"method" %in% names_local) {
      local$method <- "kmeans"
    }
    local$index <- get_local_estimation_index(local, data, n)
  }

  # setting var adjust
  if (!"var_adjust" %in% names_local) {
    if (n <= 1e5) {
      local$var_adjust <- "theoretical"
    } else {
      message('var_adjust was not specified and the sample size exceeds 100,000, so the default var_adjust value is being changed from "theoretical" to "none". To override this behavior, rerun and set var_adjust in local. Be aware that setting var_adjust to "theoretical" may result in exceedingly long computational times.')
      local$var_adjust <- "none"
    }

  } # "none", "empirical", "theoretical", and "pooled"

  # setting partition factor
  local$partition_factor <- partition_factor

  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$parallel) {
    n_index <- length(unique(local$index))
    if ("ncores" %in% names_local) {
      cores_available <- parallel::detectCores()
      local$ncores <- min(n_index, local$ncores, cores_available)
    } else {
      local$ncores <- parallel::detectCores()
      local$ncores <- min(n_index, local$ncores)
    }
  }

  local
}

get_local_estimation_index <- function(local, data, n) {
  if (local$method == "random") {
    index <- sample(rep(seq_len(local$groups), times = local$size)[seq_len(n)])
  } else if (local$method == "kmeans") {
    # any extra elements in local (beyond the reserved names below) are
    # forwarded to kmeans() by value, e.g. to control nstart or algorithm
    kmeans_arg_names <- setdiff(names(local), c("size", "groups", "method", "index", "parallel", "ncores", "var_adjust"))
    kmeans_args <- local[kmeans_arg_names]
    x <- st_coordinates(data)
    # modifyList() so an explicit local$iter.max overrides the default
    # instead of being passed alongside it as a second named argument
    kmeans_call_args <- utils::modifyList(list(x = x, centers = local$groups, iter.max = 30), kmeans_args)
    index <- do.call("kmeans", kmeans_call_args)$cluster
  } else {
    stop("local$method must be random (the default) or kmeans")
  }
  index
}

get_local_list_prediction <- function(local) {
  # set local neighborhood size
  # method can be "all" (for all data),
  # or "covariance" (for local covariance neighborhoods)

  if (is.logical(local)) {
    if (local) {
      local <- list(method = "covariance", size = 200, chunk_size = 1000L, parallel = FALSE)
    } else {
      local <- list(method = "all", chunk_size = 1000L, parallel = FALSE)
    }
  }

  names_local <- names(local)

  # errors
  if ("method" %in% names_local) {
    if (!local$method %in% c("all", "covariance")) {
      stop("Invalid local method. Local method must be \"all\" or \"covariance\".", call. = FALSE)
    }
  }

  if (!"method" %in% names_local) {
    # local$method <- "all"
    local$method <- "covariance"
  }

  if (local$method %in% c("covariance") && !"size" %in% names_local) {
    local$size <- 200
  }

  if (!"chunk_size" %in% names_local) {
    local$chunk_size <- 1000L
  }
  if (!is.numeric(local$chunk_size) || length(local$chunk_size) != 1 || is.na(local$chunk_size) || local$chunk_size < 1) {
    stop("local$chunk_size must be a single positive number.", call. = FALSE)
  }
  local$chunk_size <- as.integer(local$chunk_size)

  if (!"parallel" %in% names_local) {
    local$parallel <- FALSE
    local$ncores <- NULL
  }

  if (local$parallel) {
    if (!"ncores" %in% names_local) {
      local$ncores <- parallel::detectCores()
    }
  }

  local
}

get_local_list_prediction_block <- function(local) {
  if (is.logical(local)) {
    if (length(local) != 1 || is.na(local)) {
      stop("local must be TRUE, FALSE, or a local-control list.", call. = FALSE)
    }
    local <- if (local) {
      list(method = "covariance", size = 4000L, method_new = "basis", size_new = 4000L,
           ordering = "pid", chunk_size = 1000L, parallel = FALSE)
    } else {
      list(method = "all", method_new = "basis", size_new = Inf,
           ordering = "pid", chunk_size = 1000L, parallel = FALSE)
    }
  }
  if (!is.list(local)) {
    stop("local must be TRUE, FALSE, or a local-control list.", call. = FALSE)
  }

  names_local <- names(local)
  if (!is.null(local$method) && !local$method %in% c("all", "covariance")) {
    stop("Invalid local method. Local method must be \"all\" or \"covariance\".", call. = FALSE)
  }
  if (!is.null(local$method_new) && !local$method_new %in% c("basis", "subset")) {
    stop("local$method_new must be \"basis\" or \"subset\".", call. = FALSE)
  }

  if (!"method" %in% names_local) local$method <- "covariance"
  if (identical(local$method, "covariance") && !"size" %in% names_local) local$size <- 4000L
  if (!"method_new" %in% names_local) local$method_new <- "basis"
  if (!"size_new" %in% names_local) local$size_new <- 4000L
  local$ordering <- get_decorrelate_ordering(local$ordering)
  if (!"chunk_size" %in% names_local) local$chunk_size <- 1000L
  if (!"parallel" %in% names_local) local$parallel <- FALSE

  numeric_names <- c("size_new", "chunk_size")
  if (identical(local$method, "covariance") || "size" %in% names(local)) {
    numeric_names <- c("size", numeric_names)
  }
  for (name in numeric_names) {
    value <- local[[name]]
    if (is.null(value) || length(value) != 1 || !is.numeric(value) || is.na(value) ||
        value < 1 || (name != "size_new" && !is.finite(value))) {
      stop("local$", name, " must be a single positive number.", call. = FALSE)
    }
    if (is.finite(value)) local[[name]] <- as.integer(value)
  }
  if (!is.logical(local$parallel) || length(local$parallel) != 1 || is.na(local$parallel)) {
    stop("local$parallel must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(local$parallel)) {
    if (!"ncores" %in% names_local) local$ncores <- parallel::detectCores()
    if (!is.numeric(local$ncores) || length(local$ncores) != 1 || is.na(local$ncores) || local$ncores < 1) {
      stop("local$ncores must be a single positive number.", call. = FALSE)
    }
    local$ncores <- as.integer(local$ncores)
  } else {
    local$ncores <- NULL
  }
  local
}
