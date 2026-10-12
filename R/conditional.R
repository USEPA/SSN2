#' Conditionally simulate from a model
#'
#' @description Conditionally simulate prediction data from a model object.
#'
#' @param object A fitted model object.
#' @param newdata A character string naming a single prediction data set
#'   (accessible via \code{object$ssn.object$preds}). May be omitted when
#'   \code{"newdata"} is not requested; otherwise, omission uses the
#'   \code{.missing} prediction set, if present.
#'   \code{"all"}/multiple prediction sets are not supported
#'   (a joint draw is only meaningful over one coherent set of locations).
#' @param output The output type, which can be any subset of
#'   \code{c("newdata", "beta", "object")} (\code{"all"} is shorthand for this full vector).
#'   For generalized linear models, \code{"latent"} returns simulated observed-site
#'   link values, including offsets..
#'   \code{"random"} returns realized random effects for models with random effects. For [ssn_lm()] fits only, when
#'   \code{simulate_covparams = TRUE}, \code{output} can also include
#'   \code{"cov"}, \code{"ssn"}, \code{"tailup"}, \code{"taildown"},
#'   \code{"euclid"}, \code{"nugget"}, or \code{"randcov"}, matching
#'   [vcov.SSN2()]'s \code{type} values -- see \code{simulate_covparams}
#'   below. The default is \code{"newdata"}. See Details for more.
#' @param samples The number of conditional simulations. The default is
#'   \code{1,000}.
#' @param local For \code{ssn_lm()} and \code{ssn_glm()} model objects, an optional logical or list
#'   controlling the big data approximation.
#'   If omitted, \code{local} is set to \code{TRUE} when the observed or
#'   requested prediction sample size exceeds 5,000, and \code{FALSE} otherwise.
#'   If \code{FALSE}, no big data approximation is implemented.
#'   Before using local approximations, it is useful to create \code{.bmat}
#'   distance matrices using \code{ssn_create_bigdist()}.
#'   If a list is provided, \code{local$approximation} selects which big data
#'   approximation is used and can take on the values
#'   \code{"low-rank"} or \code{"vecchia"} (the default):
#'   \itemize{
#'     \item \code{"low-rank"}: a base sample is drawn from the observed
#'       data, \code{newdata} is split into blocks, and each block is
#'       simulated conditional on the base sample (blocks are
#'       assumed conditionally independent of one another given the base
#'       sample, coefficients, and realized random effects). The base-sample settings (\code{method_base}/\code{size_base}/
#'       \code{reorder_base}) and the \code{newdata}-blocking settings
#'       (\code{method_new}/\code{size_new}/\code{reorder_new}/\code{kmeans_new})
#'       are separate from one another.
#'       \itemize{
#'         \item \code{method_base}: Whether the observed data conditioned on is
#'           restricted to a base sample. If \code{method_base = "all"}, no
#'           big data approximation is applied to the observed data (all of it is
#'           conditioned on). If \code{method_base = "base"}, the observed data
#'           is subset to \code{size_base} observations (ordered via
#'           \code{reorder_base}) to form the base sample. The default is
#'           \code{"base"}.
#'         \item \code{reorder_base}: The data reordering approach used to reorder
#'           the observed data prior to subsetting to obtain a base sample.
#'           If \code{reorder_base = "none"}, no reordering
#'           is applied to the observed data. If \code{reorder_base = "random"}, the observed data order is
#'           randomly reshuffled. If \code{reorder_base = "grts"}, the observed data order is
#'           randomly generated using the GRTS algorithm for spatially balanced
#'           sampling via \code{spsurvey::grts()}. The default is \code{"grts"}.
#'         \item \code{size_base}: The number of observed data observations used for the base sample.
#'           The default is 5,000. See Details for more.
#'         \item \code{method_new}: Whether \code{newdata} is split into blocks
#'           for simulation. If \code{method_new = "all"}, no big data
#'           approximation is applied to \code{newdata} (it is simulated all at
#'           once). If \code{method_new = "base"}, \code{newdata} is split into
#'           blocks of (approximately) \code{size_new} observations each (ordered
#'           via \code{reorder_new}, and optionally grouped via \code{kmeans_new}),
#'           with each block simulated conditional on the base sample, coefficients, and realized random effects
#'           independently of every other block. The default is \code{"base"}.
#'         \item \code{reorder_new}: The data reordering approach used to reorder
#'           \code{newdata} before splitting it into blocks. If \code{reorder_new = "none"}, no reordering
#'           is applied to \code{newdata}. If \code{reorder_new = "random"}, \code{newdata} is
#'           randomly reshuffled. The default is \code{"random"}.
#'         \item \code{kmeans_new}: Whether \code{newdata} observations
#'           should be assigned to blocks based on k-means clustering
#'           on the coordinates, with clusters of size approximately equal to
#'           \code{size_new}. The default is \code{FALSE} when \code{reorder_new = "none"}
#'           and \code{TRUE} otherwise.
#'         \item \code{size_new}: The (approximate) number of observations used
#'           for each block. The default is 1,000. See Details for more.
#'         \item \code{chunk_size}: For big data \code{ssn_lm()} and \code{ssn_glm()}
#'           simulation, the maximum number of \code{samples} columns processed simultaneously.
#'           The default is 1,000. \code{chunk_size} is purely computational
#'           and does not affect correctness.
#'         \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'           parallel package is automatically used. The default is \code{FALSE}.
#'           Only used when \code{local$approximation} is \code{"low-rank"}.
#'         \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'           parallelize over. The default is the number of available cores on your machine.
#'       }
#'       If \code{local$approximation} is \code{"low-rank"}, defaults for the remaining \code{"low-rank"}
#'       settings are chosen such that \code{local} is transformed into
#'       \code{list(approximation = "low-rank", method_base = "base", size_base = 5000,
#'       reorder_base = "grts", method_new = "base", size_new = 1000,
#'       reorder_new = "random", kmeans_new = TRUE, parallel = FALSE)}.
#'     \item \code{"vecchia"}: every \code{newdata} location is simulated one
#'       at a time (in some order over \code{newdata}), each conditional on
#'       neighbors selected from all observed data plus each already-simulated
#'       \code{newdata} location (not a single shared base sample).
#'       Unlike low-rank simulation, prediction blocks are not assumed conditionally
#'       independent given a shared base. This matches \code{local = FALSE} when
#'       \code{method = "all"}; \code{method = "covariance"} truncates the
#'       conditioning set to a fixed number of neighbors sorted by absolute covariance.
#'       No parallelization exists because the algorithm is inherently sequential,
#'       as each new observation can depend on previous ones.
#'       For generalized linear models, the observed latent field is drawn first using local normal (i.e., Gaussian) conditionals.
#'       \itemize{
#'         \item \code{method}: The neighbor-selection rule used to build each
#'           location's conditioning set once it exceeds \code{size} candidates
#'           (all observed data plus every already-simulated \code{newdata}
#'           location). Values are \code{"all"} or \code{"covariance"}
#'           (the \code{size} candidates with the highest absolute spatial covariance).
#'           True partitions restrict eligible candidates. The default is
#'           \code{"covariance"}. For both linear and generalized linear models,
#'           \code{method = "all"} uses the exact path, as with \code{local = FALSE}.
#'         \item \code{size}: The number of neighbors used when \code{method}
#'           is \code{"covariance"}. The default is 30 for
#'           \code{ssn_lm()} and 60 for \code{ssn_glm()}.
#'         \item \code{ordering}: The order \code{newdata} locations are
#'           simulated in -- \code{"pid"}, \code{"maxmin"}, \code{"middleout"},
#'           \code{"outsidein"}, \code{"coordinate"}, \code{"grts"},
#'           \code{"random"}, or \code{"none"} (same options as
#'           \code{ssn_decorrelate()}'s \code{ordering} argument). The default is
#'           \code{"pid"}, which orders by \code{pid} and then by network ID.
#'           The same rule is applied separately to observed generalized linear model sites.
#'           \code{"maxmin"}, \code{"middleout"}, \code{"outsidein"}, and
#'           \code{"coordinate"} require \pkg{GPvecchia}; \code{"grts"} requires \pkg{spsurvey}.
#'         \item \code{method_random}: The method used to assign observations to
#'           spatial covariance blocks when conditioning shared random effects.
#'           If \code{method_random = "random"}, observations are randomly assigned
#'           to blocks. If \code{method_random = "kmeans"}, observations are assigned
#'           based on k-means clustering on the coordinates with \code{groups_random}
#'           clusters. The default is \code{"kmeans"}.
#'           Both grouping methods have a random component; use \code{set.seed()}
#'           for reproducible assignments. These blocks are separate from model-fitting groups.
#'           Random-effect block options require an active Vecchia approximation
#'           and a model with random effects.
#'         \item \code{size_random}: The target number of observations in each
#'           random-effect block. The default is 1,000 when neither \code{size_random}
#'           nor \code{groups_random} is supplied.
#'         \item \code{groups_random}: The number of random-effect blocks. An explicit
#'           \code{size_random} takes precedence and sets \code{groups_random} to
#'           \code{ceiling(n / size_random)}, where \code{n} is the observed sample size.
#'           If only \code{groups_random} is supplied, \code{size_random} is
#'           \code{ceiling(n / groups_random)}.
#'         \item \code{chunk_size}: For local \code{ssn_lm()} and \code{ssn_glm()}
#'           simulation, with or without random effects, the
#'           maximum number of \code{samples} columns processed simultaneously.
#'           The default is 1,000. \code{chunk_size} is purely computational and does
#'           not affect correctness.
#'         \item \code{parallel}: Not used for Vecchia; \code{TRUE} warns and is ignored.
#'       }
#'   }
#'       With \code{local = TRUE} or \code{local = list()}, the default is
#'       \code{list(approximation = "vecchia", method = "covariance")},
#'       with \code{size = 30} for \code{ssn_lm()} and \code{size = 60} for \code{ssn_glm()}.
#'       Other lists without \code{approximation} also use Vecchia.
#'       For both approximations, \code{chunk_size} must be a positive integer
#'       and is separate from neighborhood and prediction-block sizes.
#'       Low-rank base and block settings require \code{approximation = "low-rank"};
#'       Vecchia settings (\code{method}, \code{size}, and \code{ordering}) are ignored
#'       when \code{approximation = "low-rank"}.
#' @param simulate_covparams For [ssn_lm()] model objects, whether to also
#'   simulate new covariance parameters for each sample. \code{simulate_covparams}
#'   requires \code{object$vcov$cov} to be specified during model fitting by
#'   selecting \code{ddf = "satterthwaite"}. \code{simulate_covparams = TRUE}
#'   should generally not be used for sample sizes greater than 500
#'   given its computational inefficiencies. The default is
#'   \code{FALSE}.
#' @param ... Other arguments. Not used (needed for generic consistency).
#'
#' @details
#'
#'   If \code{"newdata"} is in \code{output},
#'   conditional simulations are returned for each row of the named prediction set.
#'   If \code{"beta"} is in \code{output},
#'   conditional simulations are returned for each fixed effect
#'   (i.e., element of \code{coef(object)}). If \code{"object"} is in \code{output},
#'   observed responses (\code{ssn_lm()}) or fitted link values (including offsets)
#'   (\code{ssn_glm()}) are returned in observed-row order,
#'   repeated across simulation columns. If \code{"latent"} is requested, the
#'   actual observed link-scale draws are returned, with offsets included.
#'   \code{type} affects only \code{newdata}.
#'   For example, \code{c("newdata", "beta")} returns
#'   the conditional simulations both for \code{newdata} and for the
#'   fixed effects. If a covariance-parameter name
#'   (\code{"cov"}, \code{"ssn"}, \code{"tailup"}, \code{"taildown"},
#'   \code{"euclid"}, \code{"nugget"}, or \code{"randcov"}) is in
#'   \code{output} (only available when \code{simulate_covparams = TRUE}),
#'   the simulated covariance-parameter draws themselves are returned.
#'
#'   \code{local} Details: When \code{local$approximation} is \code{"low-rank"}, the big
#'   data approximation works by assigning \code{size_base}
#'   observations to a base sample and then simulating data for the base sample.
#'   The remaining observations are assigned to blocks. For each block, data
#'   are simulated from the conditional distribution given the base sample.
#'   Observations from the same block share conditional covariance while
#'   observations from distinct blocks are assumed conditionally independent
#'   (given the base sample). Parallelization generally further speeds up
#'   computations. When \code{local$approximation} is \code{"vecchia"}, no such
#'   independence assumption is made -- see the \code{local} argument above
#'   for details.
#'
#' @return If \code{output = "newdata"}, an a x b matrix of conditional simulations
#'   for each row in \code{newdata}, where a is the
#'   number of rows in \code{newdata} and b is the number of samples.
#'   If \code{output = "beta"}, an p x b matrix of conditional simulations for each
#'   element in \code{coef(object)}, where p is the
#'   number of fixed effects and b is the number of samples.
#'   If \code{output = "object"}, an n x b matrix of observed responses
#'   (\code{ssn_lm()}) or fitted link values (including offsets)
#'   (\code{ssn_glm()}), where n is the number of observed rows and b is the number of samples. If
#'   \code{output} is \code{"cov"}, \code{"ssn"}, \code{"tailup"},
#'   \code{"taildown"}, \code{"euclid"}, \code{"nugget"}, or
#'   \code{"randcov"}, a (covariance parameter) x \code{samples} matrix of
#'   the simulated covariance parameters is returned instead --
#'   \code{"cov"} returns every free covariance parameter, \code{"ssn"} the
#'   tailup, taildown, Euclidean, and nugget subset,
#'   \code{"tailup"}/\code{"taildown"}/\code{"euclid"}/\code{"nugget"} each
#'   their own individual subset, and \code{"randcov"} the random-effect
#'   subset. \code{output = "latent"} returns an n x b matrix of simulated observed
#'   link values (including offsets), regardless of \code{type}.
#'   \code{output = "random"} returns a q x b matrix of realized random-effect coefficients,
#'   separated by partition factor (when relevant).
#'
#'   If \code{output} has more than one element, a list is returned with the
#'   respective elements named according to the relevant \code{output}. For example,
#'   \code{output = c("newdata", "beta")} returns a list with elements
#'   \code{"newdata"} and \code{"beta"}, each containing the relevant output
#'   for \code{output = "newdata"} and \code{output = "beta"}, respectively.
#'
#' @name conditional.SSN2
#' @method conditional ssn_lm
#' @order 1
#' @export
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, predpts = "CapeHorn", overwrite = TRUE)
#' ssn_create_distmat(mf04p, predpts = "CapeHorn", overwrite = TRUE, among_predpts = TRUE)
#'
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' set.seed(1)
#' draws <- conditional(ssn_mod, "CapeHorn")
#' head(draws) # rows are prediction locations; columns are simulations
conditional.ssn_lm <- function(object, newdata, output = "newdata", samples = 1000, local, simulate_covparams = FALSE, ...) {
  if (!is.numeric(samples) || length(samples) != 1L || !is.finite(samples) || samples < 1 || samples != floor(samples)) {
    stop("samples must be a positive integer.", call. = FALSE)
  }
  if (missing(local)) local <- NULL
  if (!is.logical(simulate_covparams) || length(simulate_covparams) != 1 || is.na(simulate_covparams)) {
    stop("simulate_covparams must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(simulate_covparams) && !is.null(object$local_index)) {
    stop("simulate_covparams = TRUE is not supported for models fitted with 'local'; use simulate_covparams = FALSE to simulate with fitted covariance parameters held fixed.", call. = FALSE)
  }
  covparam_outputs <- c("cov", "ssn", "tailup", "taildown", "euclid", "nugget", "randcov")
  if ("all" %in% output) {
    output <- c("newdata", "beta", "object", if ("random" %in% output) "random")
  }
  output_valid <- is.character(output) && length(output) > 0 && !anyNA(output) &&
    all(output %in% c("newdata", "beta", "object", "random", covparam_outputs))
  if (!output_valid) {
    stop("output must be any subset of \"newdata\", \"beta\", \"object\", \"random\" (or \"cov\", \"ssn\", \"tailup\", \"taildown\", \"euclid\", \"nugget\", or \"randcov\" when simulate_covparams = TRUE), or \"all\".", call. = FALSE)
  }
  if ("random" %in% output && is.null(object$random)) {
    stop('output = "random" requires a model with random effects.', call. = FALSE)
  }
  newdata_name <- if ("newdata" %in% output) resolve_conditional_newdata_name(object, newdata) else NULL
  conditioning <- get_local_list_conditional(local, object,
    if (!is.null(newdata_name)) object$ssn.object$preds[[newdata_name]] else NULL)
  if (!conditioning$exact && isTRUE(simulate_covparams)) {
    simulate_covparams <- FALSE
    message("simulate_covparams = TRUE is not used when a big-data approximation (local) is specified; setting simulate_covparams = FALSE.")
  }
  if (any(output %in% covparam_outputs) && !isTRUE(simulate_covparams)) {
    stop("output can only be \"cov\", \"ssn\", \"tailup\", \"taildown\", \"euclid\", \"nugget\", or \"randcov\" when simulate_covparams = TRUE.", call. = FALSE)
  }

  if (!conditioning$exact) {
    if (identical(conditioning$approximation, "low-rank")) {
      val <- draw_conditional_lowrank(object, newdata_name, conditioning, samples, output)
      if ("random" %in% output) val$random <- check_conditional_random_output(val$random)
      return(if (length(output) == 1) val[[output]] else val[output])
    }
    val <- draw_conditional_vecchia_model(object, newdata_name, conditioning, samples, output)
    if ("random" %in% output) val$random <- check_conditional_random_output(val$random)
    return(if (length(output) == 1) val[[output]] else val[output])
  }

  if (!isTRUE(simulate_covparams)) {
    val <- draw_conditional_exact(object, newdata_name, samples, output)
    if ("random" %in% output) val$random <- check_conditional_random_output(val$random)
    return(if (length(output) == 1) val[[output]] else val[output])
  }
  context <- get_conditional_context(object, newdata_name, local = TRUE)

  if (object$n > 500) {
    warning("simulate_covparams = TRUE may result in exceedingly long computational times when observed data sample sizes are greater than 500.", call. = FALSE)
  }
  sw <- get_satterthwaite_cached(object, method = "numeric")
  if (is.null(sw$vcov_theta)) {
    stop("simulate_covparams = TRUE requires a numerically positive definite covariance-parameter Hessian, which was not available for this fit (see vcov(object, type = \"cov\")).", call. = FALSE)
  }
  val <- draw_conditional_covparams(object, context, sw, samples, output)
  if ("random" %in% output) val$random <- check_conditional_random_output(val$random)
  if (length(output) == 1) val[[output]] else val[output]
}

#' @rdname conditional.SSN2
#' @method conditional ssn_glm
#' @order 2
#' @export
#'
#' @param type For \code{ssn_glm()} model objects, the scale of the conditional
#'   simulations for \code{newdata}.
#'   When \code{type = "link"}, the predicted means on the
#'   link scale are returned. When \code{type = "response"}, the predicted means
#'   on the response scale are returned. When \code{type = "new"}, a new observation
#'   is simulated from the appropriate response distribution with mean equal to
#'   the mean on the response scale and dispersion equal to the dispersion parameter
#'   from \code{object}. The default is \code{"link"}.
#' @param newdata_size The \code{size} value for each observation in \code{newdata}
#'   used when predicting for the binomial family, with a default value of 1.
conditional.ssn_glm <- function(object, newdata, output = "newdata", type = c("link", "response", "new"), samples = 1000, local, newdata_size, ...) {
  if (!is.numeric(samples) || length(samples) != 1L || !is.finite(samples) || samples < 1 || samples != floor(samples)) {
    stop("samples must be a positive integer.", call. = FALSE)
  }
  if (missing(local)) local <- NULL
  if (missing(newdata_size)) newdata_size <- NULL
  if ("all" %in% output) {
    output <- c("newdata", "beta", "object", if ("latent" %in% output) "latent", if ("random" %in% output) "random")
  }
  if (!is.character(output) || !length(output) || anyNA(output) || any(!output %in% c("newdata", "beta", "object", "latent", "random"))) {
    stop('output must be "newdata", "beta", "object", "latent", "random", or "all".', call. = FALSE)
  }
  if ("random" %in% output && is.null(object$random)) {
    stop('output = "random" requires a model with random effects.', call. = FALSE)
  }
  type <- match.arg(type)
  newdata_name <- if ("newdata" %in% output) resolve_conditional_newdata_name(object, newdata) else NULL
  conditioning <- get_local_list_conditional(local, object,
    if (!is.null(newdata_name)) object$ssn.object$preds[[newdata_name]] else NULL)
  if (conditioning$exact) {
    val <- draw_conditional_exact(object, newdata_name, samples, output, type, newdata_size)
    if ("random" %in% output) val$random <- check_conditional_random_output(val$random)
    return(if (length(output) == 1) val[[output]] else val[output])
  }
  if (identical(conditioning$approximation, "low-rank")) {
    val <- draw_conditional_lowrank(object, newdata_name, conditioning, samples, output, type, newdata_size)
    if ("random" %in% output) val$random <- check_conditional_random_output(val$random)
    return(if (length(output) == 1) val[[output]] else val[output])
  }
  val <- draw_conditional_vecchia_model(object, newdata_name, conditioning, samples, output, type, newdata_size)
  if ("random" %in% output) val$random <- check_conditional_random_output(val$random)
  if (length(output) == 1) val[[output]] else val[output]
}

#' Reject unsupported model classes for \code{conditional()}
#'
#' @param object A fitted model object.
#'
#' @return \code{NULL}, invisibly, if \code{object} is an \code{ssn_lm} or
#'   \code{ssn_glm} object; otherwise an error.
#'
#' @noRd
validate_conditional_scope <- function(object) {
  if (!inherits(object, "ssn_lm") && !inherits(object, "ssn_glm")) {
    stop("conditional() is only implemented for ssn_lm() and ssn_glm() model objects.", call. = FALSE)
  }
  invisible(NULL)
}

#' Resolve and validate \code{conditional()}'s single required \code{newdata} name
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata The name of one prediction set in
#'   \code{object$ssn.object$preds}.
#'
#' @return The resolved, validated single prediction-set name.
#'
#' @noRd
resolve_conditional_newdata_name <- function(object, newdata) {
  if (missing(newdata) || is.null(newdata)) {
    if (!".missing" %in% names(object$ssn.object$preds)) {
      stop("No missing data to predict. newdata must be specified in the newdata argument or object$ssn.object$preds$.missing must be non-NULL.", call. = FALSE)
    }
    newdata <- ".missing"
  }
  resolved <- resolve_newdata_name(object, newdata)
  if (length(resolved) != 1) {
    stop("conditional() requires a single named prediction set (not \"all\" or multiple sets); resolved to ", length(resolved), " sets.", call. = FALSE)
  }
  if (!resolved %in% names(object$ssn.object$preds)) {
    stop("\"", resolved, "\" is not a valid prediction set name in object$ssn.object$preds.", call. = FALSE)
  }
  resolved
}

#' Build the model-type-agnostic pieces of a \code{conditional()} context
#'
#' Assembles \code{newdata}/its model matrix/offset and the fitted model's
#' own design matrix, shared by both Gaussian and GLM conditional simulation.
#' When \code{local} is not \code{TRUE}, also builds the dense exact-path
#' pieces (observed-data covariance factor, observed-to-prediction and
#' prediction-to-prediction covariance) used by
#' \code{\link{get_conditional_cov}()}; local (big-data) conditional
#' simulation builds these pieces itself, per-row, in
#' \code{\link{prepare_conditional_vecchia}()} instead.
#'
#' @param object A fitted \code{ssn_lm}/\code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set to simulate.
#' @param local Whether a big-data approximation is active (skips the dense
#'   exact-path pieces when \code{TRUE}).
#'
#' @return A list with \code{object}, \code{newdata_name}, \code{newdata},
#'   \code{x0}, \code{newdata_offset}, \code{Xmat}, \code{add_newdata_rows},
#'   and (when \code{local} is not \code{TRUE}) \code{cov_lowchol_base},
#'   \code{C0}, \code{SqrtSigInv_X}, \code{SqrtSigInv_C0}, \code{Sigma22}.
#'
#' @noRd
get_conditional_base_context <- function(object, newdata_name, local = FALSE) {
  validate_conditional_scope(object)

  Xmat <- model.matrix(object)
  if (is.null(newdata_name)) {
    return(list(object = object, newdata_name = NULL,
      newdata = object$ssn.object$obs[FALSE, , drop = FALSE],
      x0 = Xmat[FALSE, , drop = FALSE], Xmat = Xmat,
      newdata_offset = NULL, add_newdata_rows = FALSE))
  }
  pn <- get_prediction_newdata(object, newdata_name)
  newdata <- pn$newdata
  add_newdata_rows <- pn$add_newdata_rows

  if (NROW(newdata) == 0) {
    stop("The prediction set \"", newdata_name, "\" has no rows.", call. = FALSE)
  }

  nm <- get_newdata_model_matrix(object, newdata)
  newdata <- nm$newdata
  x0 <- nm$newdata_model
  newdata_offset <- nm$offset

  base <- list(
    object = object, newdata_name = newdata_name, newdata = newdata,
    x0 = x0, newdata_offset = newdata_offset,
    Xmat = Xmat,
    add_newdata_rows = add_newdata_rows
  )

  if (isTRUE(local)) {
    return(base)
  }

  cov_matrix_val <- covmatrix(object)
  cov_lowchol_base <- t(chol(cov_matrix_val))

  C0 <- covmatrix(object, newdata_name, cov_type = "obs.pred")
  Sigma22 <- covmatrix(object, newdata_name, cov_type = "pred.pred")

  SqrtSigInv_X <- forwardsolve(cov_lowchol_base, Xmat)
  SqrtSigInv_C0 <- forwardsolve(cov_lowchol_base, C0)

  c(base, list(
    cov_lowchol_base = cov_lowchol_base,
    C0 = C0,
    SqrtSigInv_X = SqrtSigInv_X, SqrtSigInv_C0 = SqrtSigInv_C0,
    Sigma22 = Sigma22
  ))
}

#' Build a Gaussian \code{conditional()} context
#'
#' Extends \code{\link{get_conditional_base_context}()} with the observed
#' response/offset, fitted coefficients, and their covariance, needed for
#' Gaussian conditional simulation.
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param newdata_name The name of the prediction set to simulate.
#' @param local Whether a big-data approximation is active.
#'
#' @return The \code{\link{get_conditional_base_context}()} list, extended
#'   with \code{y_object}, \code{y} (offset-subtracted), \code{betahat}, and
#'   \code{cov_betahat}.
#'
#' @noRd
get_conditional_context <- function(object, newdata_name, local = FALSE) {
  context <- get_conditional_base_context(object, newdata_name, local = local)

  y <- model.response(model.frame(object))
  context$y_object <- y
  offset <- model.offset(model.frame(object))
  if (!is.null(offset)) {
    y <- y - offset
  }

  context$y <- y
  context$betahat <- coefficients(object)
  context$cov_betahat <- vcov(object)
  context
}

#' Build a GLM \code{conditional()} context
#'
#' Extends \code{\link{get_conditional_base_context}()} with the observed
#' response/offset, fitted latent linear predictor, coefficients and their
#' (both corrected and uncorrected) covariance, family, dispersion, and
#' binomial size, needed for GLM conditional simulation.
#'
#' @param object A fitted \code{ssn_glm} model object.
#' @param newdata_name The name of the prediction set to simulate.
#' @param local Whether a big-data approximation is active.
#'
#' @return The \code{\link{get_conditional_base_context}()} list, extended
#'   with \code{y}, \code{w} (fitted link-scale values), \code{w_free}
#'   (offset-subtracted \code{w}), \code{offset}, \code{betahat},
#'   \code{cov_betahat}, \code{cov_betahat_uncorrected}, \code{family},
#'   \code{dispersion}, and \code{size}.
#'
#' @noRd
get_conditional_context_glm <- function(object, newdata_name, local = FALSE) {
  context <- get_conditional_base_context(object, newdata_name, local = local)

  y <- model.response(model.frame(object))
  offset <- model.offset(model.frame(object))
  w <- fitted(object, type = "link")
  w_free <- if (!is.null(offset)) w - offset else w

  context$y <- y
  context$w <- w
  context$w_free <- w_free
  context$offset <- offset
  context$betahat <- coefficients(object)
  context$cov_betahat <- vcov(object)
  context$cov_betahat_uncorrected <- vcov(object, var_correct = FALSE)
  context$family <- object$family
  context$dispersion <- as.vector(coef(object, type = "dispersion"))
  context$size <- object$size
  context
}

#' Compute the exact-path conditional covariance among prediction locations
#'
#' @param context A conditional-simulation context from
#'   \code{\link{get_conditional_context}()}/
#'   \code{\link{get_conditional_context_glm}()}.
#' @return A list with \code{Sigma_cond} (the \code{m x m} conditional
#'   covariance matrix) and \code{chol_cond_cov} (its lower Cholesky factor).
#'
#' @noRd
get_conditional_cov <- function(context) {
  Sigma_cond <- context$Sigma22 - crossprod(context$SqrtSigInv_C0)
  Sigma_cond <- as.matrix(Matrix::forceSymmetric(Sigma_cond))
  chol_cond_cov <- t(chol(Sigma_cond))
  list(Sigma_cond = Sigma_cond, chol_cond_cov = chol_cond_cov)
}

#' Check whether a simulated covariance-parameter value is within its valid range
#'
#' Used by \code{\link{simulate_theta_draw_ssn}()}'s reject-and-redraw loop
#' for \code{simulate_covparams = TRUE}.
#'
#' @param name The parameter's free-parameter name (e.g. \code{"tailup_de"},
#'   \code{"euclid_rotate"}, \code{"nugget"}, or a random-effect variance
#'   name).
#' @param value The simulated value to check.
#' @param euclid_type The Euclidean covariance type, used to bound
#'   \code{"euclid_extra"}.
#'
#' @return \code{TRUE} if \code{value} is finite and within the field's valid
#'   range, \code{FALSE} otherwise.
#'
#' @noRd
is_valid_covparam_field <- function(name, value, euclid_type) {
  if (!is.finite(value)) {
    return(FALSE)
  }
  if (grepl("_de$", name)) {
    return(value >= 0)
  }
  if (grepl("_range$", name)) {
    return(value > 0)
  }
  if (identical(name, "euclid_rotate")) {
    return(value >= 0 && value <= pi)
  }
  if (identical(name, "euclid_scale")) {
    return(value > 0 && value <= 1)
  }
  if (identical(name, "euclid_extra")) {
    if (identical(euclid_type, "matern")) {
      return(value >= 0.2 && value <= 5)
    }
    if (identical(euclid_type, "cauchy")) {
      return(value > 0)
    }
    if (identical(euclid_type, "pexponential")) {
      return(value > 0 && value <= 2)
    }
    return(TRUE)
  }
  if (identical(name, "nugget")) {
    return(value >= 0)
  }
  # anything else is a random-effect variance component (spmodel's own
  # randcov_params() does not validate non-negativity either; the same
  # explicit check spmodel's try_build_theta() adds is applied here)
  value >= 0
}

#' Clamp a simulated covariance-parameter value to its nearest valid boundary
#'
#' Used by \code{\link{simulate_theta_draw_ssn}()} once its reject-and-redraw
#' loop exhausts \code{max_attempts}; the field-by-field boundary logic
#' mirrors \code{\link{is_valid_covparam_field}()}'s own field-specific valid
#' ranges.
#'
#' @param name The parameter's free-parameter name.
#' @param value The simulated value to clamp.
#' @param euclid_type The Euclidean covariance type, used to bound
#'   \code{"euclid_extra"}.
#'
#' @return The clamped value.
#'
#' @noRd
clamp_covparam_field <- function(name, value, euclid_type) {
  if (grepl("_de$", name)) {
    return(max(0, value))
  }
  if (grepl("_range$", name)) {
    return(max(0, value))
  }
  if (identical(name, "euclid_rotate")) {
    return(min(max(value, 0), pi))
  }
  if (identical(name, "euclid_scale")) {
    return(min(max(value, .Machine$double.eps), 1))
  }
  if (identical(name, "euclid_extra")) {
    if (identical(euclid_type, "matern")) {
      return(min(max(value, 0.2), 5))
    }
    if (identical(euclid_type, "cauchy")) {
      return(max(value, .Machine$double.eps))
    }
    if (identical(euclid_type, "pexponential")) {
      return(min(max(value, .Machine$double.eps), 2))
    }
    return(value)
  }
  if (identical(name, "nugget")) {
    return(max(0, value))
  }
  max(0, value)
}

#' Simulate one covariance-parameter draw via reject-and-redraw
#'
#' Draws \code{theta ~ N(theta_hat_free, vcov_theta)} repeatedly (up to
#' \code{max_attempts} times) until every free covariance-parameter field is
#' valid (see \code{\link{is_valid_covparam_field}()}); if every attempt
#' fails, clamps the last draw to its nearest valid boundary via
#' \code{\link{clamp_covparam_field}()} instead.
#'
#' @param theta_hat_free The fitted free covariance-parameter values (on the
#'   original, not optimization, scale).
#' @param vcov_theta_lowchol The lower Cholesky factor of the free
#'   covariance-parameter estimates' variance-covariance matrix.
#' @param cov_names_free Names of the free covariance-parameter fields.
#' @param euclid_type The Euclidean covariance type.
#' @param max_attempts The maximum number of reject-and-redraw attempts.
#'
#' @return A list with \code{theta} (the accepted or clamped draw),
#'   \code{exhausted} (whether every attempt was rejected), and
#'   \code{invalid_names} (which fields were invalid on the final attempt).
#'
#' @noRd
simulate_theta_draw_ssn <- function(theta_hat_free, vcov_theta_lowchol, cov_names_free,
                                     euclid_type, max_attempts = 50) {
  draw <- NULL
  for (attempt in seq_len(max_attempts)) {
    draw <- theta_hat_free + as.numeric(vcov_theta_lowchol %*% rnorm(length(theta_hat_free)))
    names(draw) <- cov_names_free
    valid <- vapply(cov_names_free, function(nm) is_valid_covparam_field(nm, draw[[nm]], euclid_type), logical(1))
    if (all(valid)) {
      return(list(theta = draw, exhausted = FALSE, invalid_names = character(0)))
    }
  }

  invalid_names <- cov_names_free[!valid]
  clamped <- vapply(cov_names_free, function(nm) clamp_covparam_field(nm, draw[[nm]], euclid_type), numeric(1))
  names(clamped) <- cov_names_free

  for (comp in c("tailup", "taildown", "euclid")) {
    range_name <- paste0(comp, "_range")
    de_name <- paste0(comp, "_de")
    if (range_name %in% cov_names_free && clamped[[range_name]] == 0) {
      clamped[[range_name]] <- .Machine$double.eps
      if (de_name %in% cov_names_free) {
        clamped[[de_name]] <- 0
      }
    }
  }

  list(theta = clamped, exhausted = TRUE, invalid_names = invalid_names)
}

#' Draw exact conditional samples with simulated covariance-parameter uncertainty
#'
#' Implements \code{conditional()}'s \code{simulate_covparams = TRUE} path:
#' for each replicate, simulates a covariance-parameter draw via
#' \code{\link{simulate_theta_draw_ssn}()}, refits \code{beta} and its
#' covariance under that draw, and (when requested) draws the conditional
#' \code{newdata} distribution given that replicate's \code{beta}/covariance
#' parameters.
#'
#' @param object A fitted \code{ssn_lm} model object.
#' @param context A Gaussian conditional-simulation context from
#'   \code{\link{get_conditional_context}()}.
#' @param sw A cached Satterthwaite context from
#'   \code{\link{get_satterthwaite_cached}(object, method = "numeric")}
#'   supplying \code{context} (free covariance-parameter names/values) and
#'   \code{vcov_theta} (their numeric variance-covariance matrix).
#' @param samples The number of simulated columns to draw.
#' @param output A character vector of requested outputs; any of
#'   \code{"object"}, \code{"beta"}, \code{"newdata"}, \code{"cov"},
#'   \code{"ssn"}, \code{"tailup"}, \code{"taildown"}, \code{"euclid"},
#'   \code{"nugget"}, \code{"randcov"}.
#' @param max_attempts The maximum reject-and-redraw attempts per replicate;
#'   see \code{\link{simulate_theta_draw_ssn}()}.
#'
#' @return A list with elements named by \code{output}.
#'
#' @noRd
draw_conditional_covparams <- function(object, context, sw, samples, output, max_attempts = 50) {
  sw_context <- sw$context
  vcov_theta <- sw$vcov_theta
  cov_names_free <- sw_context$cov_names_free_orig
  theta_hat_free <- sw_context$cov_val_free
  euclid_type <- remove_covtype(class(sw_context$params_object$euclid))

  vcov_theta_lowchol <- t(chol(vcov_theta))

  p <- NCOL(context$Xmat)
  m <- NROW(context$x0)
  betahat <- as.vector(context$betahat)
  need_newdata <- "newdata" %in% output
  need_random <- "random" %in% output
  need_beta <- need_newdata || need_random || "beta" %in% output
  design <- if (!is.null(object$random) && (need_newdata || need_random))
    get_conditional_random_design(object, context$newdata) else NULL
  new_random <- if (need_random) matrix(NA_real_, length(design$variance), samples) else NULL

  new_betahat <- if (need_beta) matrix(NA_real_, p, samples, dimnames = list(names(context$betahat), NULL)) else NULL
  new_val <- if (need_newdata) matrix(NA_real_, m, samples) else NULL
  new_cov <- matrix(NA_real_, length(cov_names_free), samples, dimnames = list(cov_names_free, NULL))

  exhausted_count <- 0L
  exhausted_params <- character(0)

  for (b in seq_len(samples)) {
    draw <- simulate_theta_draw_ssn(theta_hat_free, vcov_theta_lowchol, cov_names_free, euclid_type, max_attempts = max_attempts)
    if (draw$exhausted) {
      exhausted_count <- exhausted_count + 1L
      exhausted_params <- union(exhausted_params, draw$invalid_names)
    }
    new_cov[, b] <- draw$theta[cov_names_free]

    if (need_beta) {
      params_object_b <- fill_perturbed_params_object(draw$theta, sw_context)
      object_b <- object
      object_b$coefficients$params_object <- params_object_b

      de_total_b <- sum(params_object_b$tailup[["de"]], params_object_b$taildown[["de"]], params_object_b$euclid[["de"]])
      has_randcov_b <- !is.null(params_object_b$randcov)
      if (de_total_b == 0 && !has_randcov_b) {
        nugget_b <- params_object_b$nugget[["nugget"]]
        if (!(nugget_b > 0)) {
          stop(
            "A simulated covariance-parameter draw (replicate ", b, ") is degenerate ",
            "(zero total spatial variance and zero nugget); conditional simulation cannot ",
            "proceed for this replicate. This can happen when simulate_covparams = TRUE's ",
            "reject-and-redraw exhausts its attempts and clamps to a boundary far from the ",
            "fitted values -- inspect vcov(object, type = \"cov\") for ",
            "near-boundary covariance-parameter estimates.",
            call. = FALSE
          )
        }
        cov_lowchol_base_b <- Matrix::Diagonal(n = object$n, x = sqrt(nugget_b))
      } else {
        cov_lowchol_base_b <- t(chol(covmatrix(object_b)))
      }
      SqrtSigInv_X_b <- forwardsolve(cov_lowchol_base_b, context$Xmat)
      cov_betahat_b <- chol2inv(chol(as.matrix(Matrix::forceSymmetric(crossprod(SqrtSigInv_X_b)))))
      beta_b <- betahat + as.numeric(t(chol(cov_betahat_b)) %*% rnorm(p))
      new_betahat[, b] <- beta_b

      if (!is.null(design)) {
        design_b <- get_conditional_random_design(object_b, context$newdata)
        effect <- prepare_conditional_random(object_b, design_b)
        residual <- as.vector(context$y) - context$Xmat %*% beta_b
        a <- draw_conditional_random(effect, residual)
        if (need_random) new_random[, b] <- a
        if (need_newdata) {
          prediction <- prepare_conditional_exact_prediction(effect$spatial, context$newdata_name,
            t(effect$blocks[[1L]]$factor))
          new_val[, b] <- context$x0 %*% beta_b +
            draw_conditional_exact_prediction(prediction, residual - as.matrix(design_b$observed %*% a)) +
            as.matrix(design_b$newdata %*% a)
        }
      } else if (need_newdata) {
        C0_b <- covmatrix(object_b, context$newdata_name, cov_type = "obs.pred")
        Sigma22_b <- covmatrix(object_b, context$newdata_name, cov_type = "pred.pred")
        SqrtSigInv_C0_b <- forwardsolve(cov_lowchol_base_b, C0_b)
        cond_cov_b <- as.matrix(Matrix::forceSymmetric(Sigma22_b - crossprod(SqrtSigInv_C0_b)))
        chol_cond_cov_b <- t(chol(cond_cov_b))

        resid_b <- as.vector(context$y) - context$Xmat %*% beta_b
        cond_mu_b <- context$x0 %*% beta_b + crossprod(SqrtSigInv_C0_b, forwardsolve(cov_lowchol_base_b, resid_b))
        new_val[, b] <- as.numeric(cond_mu_b) + as.numeric(chol_cond_cov_b %*% rnorm(m))
      }
    }
  }

  if (exhausted_count > 0) {
    warning(
      exhausted_count, " of ", samples, " simulated covariance-parameter draws exhausted ", max_attempts, " ",
      "reject-and-redraw attempts and were clamped to their nearest valid boundary; affected ",
      "parameter(s): ", paste(exhausted_params, collapse = ", "), ".",
      call. = FALSE
    )
  }

  val <- list()
  if (need_random) val$random <- conditional_random_output(new_random, design)
  if ("object" %in% output) {
    val$object <- matrix(rep(context$y_object, times = samples), ncol = samples)
  }
  if (need_newdata) {
    if (!is.null(context$newdata_offset)) {
      new_val <- new_val + context$newdata_offset
    }
    if (context$add_newdata_rows) {
      rownames(new_val) <- context$object$missing_index
    }
    val$newdata <- new_val
  }
  if ("beta" %in% output) {
    val$beta <- new_betahat
  }
  if ("cov" %in% output) {
    val$cov <- new_cov
  }

  spcov_names_free <- intersect(get_spcov_field_names(), cov_names_free)
  if ("ssn" %in% output) {
    val$ssn <- new_cov[spcov_names_free, , drop = FALSE]
  }
  for (comp in intersect(c("tailup", "taildown", "euclid", "nugget"), output)) {
    component_names_free <- cov_names_free[startsWith(cov_names_free, comp)]
    val[[comp]] <- if (length(component_names_free) == 0) NULL else new_cov[component_names_free, , drop = FALSE]
  }
  if ("randcov" %in% output) {
    randcov_names_free <- setdiff(cov_names_free, spcov_names_free)
    val$randcov <- if (length(randcov_names_free) == 0) NULL else new_cov[randcov_names_free, , drop = FALSE]
  }

  val
}

#' Draw GLM response-family samples given the mean and dispersion
#'
#' Column-by-column family sampling shared by \code{conditional()}'s
#' \code{type = "new"} output and, indirectly, \code{ssn_rpois()}/
#' \code{ssn_rbinom()}/etc.'s own per-family sampling pattern.
#'
#' @param family The GLM family name (\code{"poisson"}, \code{"nbinomial"},
#'   \code{"binomial"}, \code{"beta"}, \code{"Gamma"}, or
#'   \code{"inverse.gaussian"}).
#' @param mu An \code{n x samples} matrix of natural-scale means
#'   (probabilities for \code{"binomial"}/\code{"beta"}, rates/means
#'   otherwise).
#' @param dispersion The family's dispersion parameter.
#' @param size The binomial trial size (only used when \code{family} is
#'   \code{"binomial"}).
#'
#' @return An \code{n x samples} matrix of response-family draws.
#'
#' @noRd
draw_glm_response <- function(family, mu, dispersion, size) {
  n <- NROW(mu)
  samples <- NCOL(mu)
  out <- matrix(NA_real_, n, samples)

  for (j in seq_len(samples)) {
    mu_j <- mu[, j]
    out[, j] <- if (family == "poisson") {
      rpois(n, mu_j)
    } else if (family == "nbinomial") {
      rnbinom(n, mu = mu_j, size = dispersion)
    } else if (family == "binomial") {
      rbinom(n, size, mu_j)
    } else if (family == "beta") {
      a <- mu_j * dispersion
      b <- (1 - mu_j) * dispersion
      val <- rbeta(n, shape1 = a, shape2 = b)
      pmin(pmax(val, 1e-4), 1 - 1e-4)
    } else if (family == "Gamma") {
      rgamma(n, shape = dispersion, scale = mu_j / dispersion)
    } else if (family == "inverse.gaussian") {
      if (!requireNamespace("statmod", quietly = TRUE)) {
        stop("Install the statmod package before using conditional() with family = \"inverse.gaussian\".", call. = FALSE)
      }
      dispersion_true <- 1 / (mu_j * dispersion)
      statmod::rinvgauss(n, mean = mu_j, dispersion = dispersion_true)
    }
  }
  out
}
