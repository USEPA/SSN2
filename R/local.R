#' Parameters that control spatial indexing (not yet implemented for spatial indexing
#'   and function calls to this just act as placeholders)
#'
#' @param local Placeholder
#' @param data Placeholder
#' @param n  Placeholder
#' @param partition_factor Placeholder
#'
#' @noRd
get_local_list_estimation <- function(local, data, n, partition_factor) {
  # size can be an integer and sets group size
  # alternatively, set the number of groups
  # index overrides size and groups
  # set var_adjust as "none", "theoretical", "empirical", and "pooled"

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
        local$size <- 100
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
    if (n <= 100000) {
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

#' A helper to get spatial indexes (not yet implemented for spatial indexing
#'   and function calls to this just act as placeholders)
#'
#' @param local Placeholder
#' @param data Placeholder
#' @param n  Placeholder
#'
#' @noRd
get_local_estimation_index <- function(local, data, n) {
  if (local$method == "random") {
    index <- sample(rep(seq_len(local$groups), times = local$size)[seq_len(n)])
  } else if (local$method == "kmeans") {
    kmeans_args <- setdiff(names(local), c("size", "groups", "method", "index", "parallel", "ncores", "var_adjust"))
    x <- st_coordinates(data)
    index <- do.call("kmeans", c(list(x = x, centers = local$groups, iter.max = 30), kmeans_args))$cluster
  } else {
    stop("local$method must be random (the default) or kmeans")
  }
  index
}


#' Parameters that control local neighborhood prediction(not yet implemented
#'   and function calls to this just act as placeholders)
#'
#' @param local Placeholder
#'
#' @noRd
get_local_list_prediction <- function(local) {
  # set local neighborhood size
  # method can be "all" (for all data),
  # or "covariance" (for local covariance neighborhoods)

  if (is.logical(local)) {
    if (local) {
      local <- list(method = "covariance", size = 2000, parallel = FALSE)
    } else {
      local <- list(method = "all", parallel = FALSE)
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
    local$size <- 2000
  }

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
