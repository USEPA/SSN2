
#' Apply the Spatial Decorrelation Transformation for Machine Learning Models
#'
#' @description Apply the spatial decorrelation transformation for SSN data,
#'   allowing for random effects, anisotropy, partition factors, and big data methods.
#'
#' @param formula A two-sided linear formula describing the fixed effect structure
#'   of the model, with the response to the left of the \code{~} operator and
#'   the terms on the right, separated by \code{+} operators. \code{.} on the
#'   right-hand side represents every variable in \code{data} except the geometry column.
#' @param ssn.object A spatial stream network object with class `SSN`.
#' @param tailup_type,taildown_type,euclid_type,nugget_type Covariance types
#'   specified as in [ssn_lm()].
#' @param tailup_params,taildown_params,euclid_params,nugget_params Known
#'   covariance parameter objects to be used in the spatial decorrelation
#'   transformation. See [tailup_params()], [taildown_params()],
#'   [euclid_params()], and [nugget_params()].
#' @param additive Additive-function variable for the tail-up covariance, as in [ssn_lm()].
#' @param anisotropy A logical indicating whether Euclidean (geometric) anisotropy should
#'   be modeled. Not required if \code{euclid_params} is provided with a nonzero
#'   \code{rotate} or a \code{scale} less than one. When \code{anisotropy} is
#'   \code{TRUE}, computational times can significantly increase. The default
#'   is \code{FALSE}.
#' @param random A one-sided linear formula describing the random effect structure
#'   of the model. Terms are specified to the right of the \code{~ operator}.
#'   Each term has the structure \code{x1 + ... + xn | g1/.../gm}, where \code{x1 + ... + xn}
#'   specifies the model for the random effects and \code{g1/.../gm} is the grouping
#'   structure. Separate terms are separated by \code{+} and must generally
#'   be wrapped in parentheses. Random intercepts are added to each model
#'   implicitly when at least one other variable is defined.
#'   If a random intercept is not desired, this must be explicitly
#'   defined (e.g., \code{x1 + ... + xn - 1 | g1/.../gm}). If only a random intercept
#'   is desired for a grouping structure, the random intercept must be specified
#'   as \code{1 | g1/.../gm}. Note that \code{g1/.../gm} is shorthand for \code{(1 | g1/.../gm)}.
#'   If only random intercepts are desired and the shorthand notation is used,
#'   parentheses can be omitted.
#' @param randcov_params An optional random effect covariance object
#'   to be used in the spatial decorrelation transformation. See [spmodel::randcov_params()].
#' @param partition_factor A one-sided linear formula with a single term
#'   specifying the partition factor. The partition factor assumes observations
#'   from different levels of the partition factor are uncorrelated.
#' @param algorithm The machine learning algorithm applied. Available options
#'   include \code{"ranger"}, \code{"randomForest"}, and \code{"xgboost"}.
#'   \code{"ranger"} specifies a random forest via [ranger::ranger()].
#'   \code{"randomForest"} specifies a random forest via [randomForest::randomForest()].
#'   \code{"xgboost"} specifies a boosted decision tree ensemble via [xgboost::xgboost()].
#' @param statistic The statistic used to evaluate fit in the test data. Available options
#'   include \code{"bias"} (mean bias), \code{"MSPE"} (mean-squared-prediction error),
#'    \code{"RMSPE"} (root-mean-squared-prediction error)
#'   and \code{"cor2"} (the predictive R-squared; i.e., the
#'   squared correlation between observations and predictions).
#' @param training A list controlling how the training and test data are assigned
#'   when evaluating test data performance.
#'   The following arguments detail this process:
#'   \itemize{
#'    \item \code{method}: The method used to evaluate test data performance.
#'      \code{"split"} will split \code{data} up into distinct training and test sets
#'      proportionally based on \code{p}. \code{"cv"} will split \code{data} up
#'      via k-fold cross validation based on \code{folds}, the number of folds.
#'    \item \code{p}: The proportion (a numeric vector between zero and one) of observations in \code{data} that should
#'      be assigned to the training data. The default is 0.8, which means that
#'      80% of the observations are assigned to the training data and 20% to the
#'      test data. Ignored if \code{training_index} or \code{test_index} are provided.
#'    \item \code{replicate}: The number of times to replicate \code{"split"} with different random training and test assignments.
#'    \item \code{folds}: The number of folds to use in cross-validation. Requires \code{method = "cv"}. Ignored if \code{folds_index} is specified. The default is 5, matching the 80/20 default split above.
#'    \item \code{training_index}: A numeric vector that specifies which rows (i.e., indices)
#'      of \code{data} should be assigned to the training data. If omitted, defaults
#'      to the rows which are not already included in \code{test_index}.
#'    \item \code{test_index}: A numeric vector that specifies which rows (i.e., indices)
#'      of \code{data} should be assigned to the test data. If omitted, defaults
#'      to the rows which are not already included in \code{training_index}.
#'    \item \code{folds_index}: A numeric vector that specifies which rows
#'      of \code{data} are associated with each cross-validation fold. Requires
#'      \code{method = "cv"}.
#'   }
#'   If omitted, \code{training} is transformed into
#'   \code{list(method = "split", p = 0.8, replicate = 1)}.
#' @param evaluate_test A logical indicating whether a grid should be constructed
#'   and evaluated when spatial decorrelation parameters are known (i.e.,
#'   \code{*_params} are all specified, and, if random effects are included, \code{randcov_params} is specified).
#'   If \code{TRUE}, constructs and evaluates the grid after assigning observations to
#'   test and training data sets. If any parameters (spatial or random effects) are estimated via a grid search,
#'   \code{evaluate_test} is set to \code{TRUE}.
#' @param ordering The data ordering applied: Available options
#'    include \code{"pid"}, \code{"grts"}, \code{"maxmin"}, \code{"middleout"},
#'   \code{"outsidein"}, \code{"coordinate"}, \code{"random"}, and \code{"none"}.
#'   \code{"pid"} applies ordering by point ID from the network geometry.
#'   \code{"grts"} applies ordering using a spatially balanced GRTS sample via \code{spsurvey::grts()}.
#'   \code{"maxmin"} applies maximum minimum distance ordering via \code{GPvecchia::order_maxmin_exact()}.
#'   \code{"middleout"} applies middle out ordering via \code{GPvecchia::order_middleout()}.
#'   \code{"outsidein"} applies middle out ordering via \code{GPvecchia::order_outsidein()}.
#'   \code{"coordinate"} applies middle out ordering via \code{GPvecchia::order_coordinate(..., coordinate = c(1, 2))},
#'   which orders from bottom-left to top-right of the spatial domain.
#'   \code{"random"} applies a completely random ordering.
#'   \code{"none"} applies no random ordering.
#'   The default is \code{"pid"}.
#' @param local An optional logical or list controlling the big data approximation.
#'   If omitted, \code{local} is set
#'   to \code{TRUE} or \code{FALSE} based on the sample size (the number of
#'   non-missing observations in \code{data}) -- if the sample size exceeds 5,000,
#'   \code{local} is set to \code{TRUE}. Otherwise it is set to \code{FALSE}.
#'   If \code{local} is \code{FALSE}, no big data approximation
#'   is implemented. If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{method}: The big data approximation method. If \code{method = "all"},
#'       all observations are used and \code{size} is ignored. 
#'       If \code{method = "covariance"}, the \code{size} data observations
#'       with the highest covariance with the observation requiring prediction are used.
#'       The default is \code{"covariance"}.
#'     \item \code{size}: The number of data observations to use when \code{method}
#'       is \code{"covariance"}. The default is 30.
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. This can significantly speed
#'       up computations even when \code{method = "all"} (i.e., no big data
#'       approximation is used), as predictions
#'       are spread out over multiple cores. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 30, method = "covariance", parallel = FALSE)}.
#' @param grid An explicit grid of parameter values by which to evaluate fit. The
#'   names of \code{grid} must contain all the names returned by \code{ssn_decorrelate_grid(formula, data, ...)}.
#' @param dense_grid A logical
#'   which controls the density of the constructed grid to be evaluated. If
#'   \code{dense_grid} is \code{TRUE}, a denser grid is used. If \code{dense_grid}
#'   is \code{FALSE}, a sparser grid is used. By default, \code{dense_grid} is \code{FALSE}.
#' @param x A fitted model from [ssn_decorrelate()].
#' @param ... Other arguments to the functions called by \code{algorithm}.
#'
#' @details
#'   The spatial decorrelation transformation is a preprocessing transformation
#'   that reduces the impacts of spatial dependence (i.e., covariance, correlation)
#'   on machine learning models. Predictions are made on the
#'   decorrelated scale and then recorrelated to account for spatial dependence.
#'   See Heaton et al., 2025 for details.
#' 
#'   \code{tailup_type} Details: Let \eqn{D} be a matrix of hydrologic distances,
#'   \eqn{W} be a diagonal matrix of weights from \code{additive}, \eqn{r = D / range},
#'   and \eqn{I} be
#'   an identity matrix. Then parametric forms for flow-connected
#'   elements of \eqn{R(zu)} are given below:
#'   \itemize{
#'     \item linear: \eqn{(1 - r) * (r <= 1) * W}
#'     \item spherical: \eqn{(1 - 1.5r + 0.5r^3) * (r <= 1) * W}
#'     \item exponential: \eqn{exp(-r) * W}
#'     \item mariah: \eqn{log(90r + 1) / 90r * (D > 0) + 1 * (D = 0) * W}
#'     \item epa: \eqn{(D - range)^2 * F * (r <= 1) * W / 16range^5}
#'     \item gaussian: \eqn{2 exp(-r^2) * (1 - pnorm(r * 2^{1/2})) * W}
#'     \item none: \eqn{I} * W
#'   }
#'
#'   Details describing the \code{F} matrix in the \code{epa} covariance are given in Garreta et al. (2010).
#'   Flow-unconnected elements of \eqn{R(zu)} are assumed uncorrelated.
#'   Observations on different networks are also assumed uncorrelated.
#'
#'   \code{taildown_type} Details: Let \eqn{D} be a matrix of hydrologic distances,
#'   \eqn{r = D / range},
#'   and \eqn{I} be an identity matrix. Then parametric forms for flow-connected
#'   elements of \eqn{R(zd)} are given below:
#'   \itemize{
#'     \item linear: \eqn{(1 - r) * (r <= 1)}
#'     \item spherical: \eqn{(1 - 1.5r + 0.5r^3) * (r <= 1)}
#'     \item exponential: \eqn{exp(-r)}
#'     \item mariah: \eqn{log(90r + 1) / 90r * (D > 0) + 1 * (D = 0)}
#'     \item epa: \eqn{(D - range)^2 * F1 * (r <= 1) / 16range^5}
#'     \item gaussian: \eqn{0}
#'     \item none: \eqn{I}
#'   }
#'
#'   Now let \eqn{A} be a matrix that contains the shorter of the two distances
#'   between two sites and the common downstream junction, \eqn{r1 = A / range},
#'   \eqn{B} be a matrix that contains the longer of the two distances between two sites and the
#'   common downstream junction, \eqn{r2 = B / range},  and \eqn{I} be an identity matrix.
#'   Then parametric forms for flow-unconnected elements of \eqn{R(zd)} are given below:
#'   \itemize{
#'     \item linear: \eqn{(1 - r2) * (r2 <= 1)}
#'     \item spherical: \eqn{(1 - 1.5r1 + 0.5r2) * (1 - r2)^2 * (r2 <= 1)}
#'     \item exponential: \eqn{0}
#'     \item mariah: \eqn{(log(90r1 + 1) - log(90r2 + 1)) / (90r1 - 90r2) * (A =/ B) + (1 / (90r1 + 1)) * (A = B)}
#'     \item epa: \eqn{(B - range)^2 * F2 * (r2 <= 1) / 16range^5}
#'     \item gaussian: \eqn{2 exp(-(B - A) / range) * (1 - pnorm(r * 2^{1/2})) * W}
#'     \item none: \eqn{I}
#'   }
#'
#'   Details describing the \code{F1} and \code{F2} matrices in the \code{epa}
#'   covariance are given in Garreta et al. (2010).
#'   Observations on different networks are assumed uncorrelated.
#'
#'  \code{euclid_type} Details: Let \eqn{D} be a matrix of Euclidean distances,
#'  \eqn{r = D / range}, and \eqn{I} be an identity matrix. Then parametric
#'  forms for elements of \eqn{R(ze)} are given below:
#'   \itemize{
#'     \item exponential: \eqn{exp(- r )}
#'     \item spherical: \eqn{(1 - 1.5r + 0.5r^3) * (r <= 1)}
#'     \item gaussian: \eqn{exp(- r^2 )}
#'     \item cubic: \eqn{(1 - 7r^2 + 8.75r^3 - 3.5r^5 + 0.75r^7) * (r <= 1)}
#'     \item pentaspherical: \eqn{(1 - 1.875r + 1.25r^3 - 0.375r^5) * (r <= 1)}
#'     \item circular: \eqn{1 - (2 / \pi) * (r * sqrt(1 - r^2) + \arcsin(r))} for \eqn{0 \le r \le 1}, and zero for \eqn{r > 1}
#'     \item wave: \eqn{sin(r) * (h > 0) / r + (h = 0)}
#'     \item jbessel: \eqn{Bj(h * range)}, Bj is Bessel-J function
#'     \item gravity: \eqn{(1 + r^2)^{-0.5}}
#'     \item rquad: \eqn{(1 + r^2)^{-1}}
#'     \item magnetic: \eqn{(1 + r^2)^{-1.5}}
#'     \item matern: \eqn{2^{1-extra} eta^{extra} K_{extra}(eta) / Gamma(extra)}, where \eqn{eta = \sqrt{2 extra} D / range}
#'     \item cauchy: \eqn{(1 + r^2)^{-extra}}
#'     \item pexponential: \eqn{exp(-D^{extra} / range)}
#'     \item none: \eqn{I}
#'   }
#'   The powered-exponential range has units of distance raised to \code{extra}.
#'
#'   \code{nugget_type} Details: Let \eqn{I} be an identity matrix and \eqn{0}
#'    be the zero matrix. Then parametric
#'    forms for elements the nugget variance are given below:
#'   \itemize{
#'     \item nugget: \eqn{I}
#'     \item none: \eqn{0}
#'   }
#'   In short, the nugget effect is modeled when \code{nugget_type} is \code{"nugget"}
#'   and omitted when \code{nugget_type} is \code{"none"}.
#'
#' \code{estmethod} Details: The various estimation methods are
#'   \itemize{
#'     \item \code{reml}: Maximize the restricted log-likelihood.
#'     \item \code{ml}: Maximize the log-likelihood.
#'   }
#'
#' \code{anisotropy} Details: By default, all Euclidean covariance parameters except \code{rotate}
#'   and \code{scale} as well as all random effect variance parameters
#'   are assumed unknown, requiring estimation. If either \code{rotate} or \code{scale}
#'   are given initial values other than 0 and 1 (respectively) or are assumed unknown
#'   in [euclid_initial()], \code{anisotropy} is implicitly set to \code{TRUE}.
#'   (Geometric) Anisotropy is modeled by transforming a Euclidean covariance function that
#'   decays differently in different directions to one that decays equally in all
#'   directions via rotation and scaling of the original Euclidean coordinates. The rotation is
#'   controlled by the \code{rotate} parameter in \eqn{[0, \pi]} radians. The scaling
#'   is controlled by the \code{scale} parameter in \eqn{[0, 1]}. The anisotropy
#'   correction involves first a rotation of the coordinates clockwise by \code{rotate} and then a
#'   scaling of the coordinates' minor axis by the reciprocal of \code{scale}. The Euclidean
#'   covariance is then computed using these transformed coordinates.
#'
#'  \code{random} Details: If random effects are used, the model
#'   can be written as \eqn{y = X \beta + Z1u1 + ... Zjuj + zu + zd + ze + n},
#'   where each Z is a random effects design matrix and each u is a random effect.
#'
#'  \code{partition_factor} Details: The partition factor can be represented in matrix form as \eqn{P}, where
#'   elements of \eqn{P} equal one for observations in the same level of the partition
#'   factor and zero otherwise. The covariance matrix involving only the
#'   spatial and random effects components is then multiplied element-wise
#'   (Hadamard product) by \eqn{P}, yielding the final covariance matrix.
#' 
#'   \code{training} Details: When \code{replicate} or the number of cross validation folds is at
#'     least two, there are separate grids evaluated for each replication (or fold). Statistics in each grid are
#'     averaged across replications (or folds) to determine a final grid ranked by \code{statistic}.
#'
#'   \code{local} Details: The big data approximation works by leveraging the
#'   conditional nature of the spatial decorrelation transformation via the
#'   Vecchia approximation. The Vecchia approximation enables efficient computation
#'   of the conditional distribution by considering only the \code{size} most relevant
#'   observations in the ordering (rather than using all the observations).
#'
#'   Observations with \code{NA} response values are removed for model
#'   fitting, but their values can be predicted afterwards by running
#'   \code{predict(object)}.
#'
#' @return A list with many elements that store information about the fitted model object:
#'   \itemize{
#'     \item \code{algorithm}: The machine learning algorithm used.
#'     \item \code{decorrelate_data}: The output of [ssn_decorrelate_data()] applied to \code{data}.
#'     \item \code{fit}: The fitted machine learning model object applied to the decorrelated data.
#'     \item \code{grid}: If used, the grid of spatial decorrelation parameters evaluated and their corresponding
#'       metrics when applied to the test data.
#'     \item \code{newdata}: The rows of \code{data} that have \code{NA} response values and are stored as prediction data.
#'     \item \code{training}: If used, the observations assigned to each training and test data set.
#'     \item \code{test}: If used, a list with the lowest (absolute) mean bias (bias), mean-squared-prediction error (MSPE),
#'       root-mean-squared-prediction error (RMSPE), and predictive R-squared (cor2).
#'   }
#'
#' @references Matthew J. Heaton, Andrew Millane, and Jake S. Rhodes. 2025. A Scalable
#'   Spatial Decorrelation Preprocessing Approach for Machine and Deep Learning.
#'   \emph{Journal of Data Science}. 1-15, DOI 10.6339/25-JDS1210
#'
#' @name ssn_decorrelate
#' @export
#' @examples
#' \donttest{
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path,
#'   predpts = "pred1km", overwrite = TRUE
#' )
#' ssn_create_distmat(mf04p, predpts = "pred1km", overwrite = TRUE)
#' fit <- ssn_decorrelate(
#'   Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential",
#'   taildown_type = "exponential", euclid_type = "exponential",
#'   additive = "afvArea", dense_grid = FALSE
#' )
#' tidy(fit$grid)
#' head(predict(fit, "pred1km"))
#' }
ssn_decorrelate <- function(formula, ssn.object,
                             tailup_type = "none", taildown_type = "none",
                             euclid_type = "none", nugget_type = "nugget",
                             tailup_params, taildown_params,
                             euclid_params, nugget_params,
                             additive, algorithm = "ranger", statistic = "RMSPE",
                             training, evaluate_test, anisotropy = FALSE,
                             random, randcov_params, partition_factor,
                             ordering,
                             local, grid, dense_grid, ...) {
  # a supplied *_params() object's own class is authoritative for type
  # (matching spmodel's decorrelate()) unless the caller names *_type
  # explicitly
  if (missing(tailup_type) && !missing(tailup_params)) tailup_type <- remove_covtype(class(tailup_params))
  if (missing(taildown_type) && !missing(taildown_params)) taildown_type <- remove_covtype(class(taildown_params))
  if (missing(euclid_type) && !missing(euclid_params)) euclid_type <- remove_covtype(class(euclid_params))
  if (missing(nugget_type) && !missing(nugget_params)) nugget_type <- remove_covtype(class(nugget_params))
  if (missing(tailup_params)) tailup_params <- NULL
  if (missing(taildown_params)) taildown_params <- NULL
  if (missing(euclid_params)) euclid_params <- NULL
  if (missing(nugget_params)) nugget_params <- NULL
  if (missing(additive)) additive <- NULL
  if (missing(random)) random <- NULL
  if (missing(randcov_params)) randcov_params <- NULL
  if (missing(partition_factor)) partition_factor <- NULL
  if (is.symbol(substitute(additive))) additive <- deparse1(substitute(additive))
  if (missing(training)) training <- NULL
  if (missing(local)) local <- NULL
  if (missing(grid)) grid <- NULL
  if (missing(dense_grid)) dense_grid <- FALSE
  if (missing(ordering)) ordering <- NULL
  check_decorrelate_dots(list(...))
  algorithm <- match.arg(algorithm, c("ranger", "randomForest", "xgboost"))
  statistic <- match.arg(statistic, c("RMSPE", "bias", "MSPE", "cor2"))
  ordering <- get_decorrelate_ordering(ordering)
  response_index <- get_decorrelate_response_index(formula, ssn.object)
  local <- get_decorrelate_local(local, length(response_index))

  pinned_tailup <- get_initial_from_params(tailup_type, tailup_params, tailup_initial)
  pinned_taildown <- get_initial_from_params(taildown_type, taildown_params, taildown_initial)
  pinned_euclid <- get_initial_from_params(euclid_type, euclid_params, euclid_initial)
  pinned_nugget <- get_initial_from_params(nugget_type, nugget_params, nugget_initial)
  randcov_initial_obj <- get_randcov_initial_from_params(randcov_params)

  initial_object <- get_initial_object(
    tailup_type, taildown_type, euclid_type, nugget_type,
    pinned_tailup, pinned_taildown, pinned_euclid, pinned_nugget
  )
  covariance_known <- get_decorrelate_covariance_known(
    initial_object, random, randcov_initial_obj, anisotropy
  )
  automatic_grid <- is.null(grid) && !covariance_known
  if (is.null(grid)) {
    if (covariance_known) {
      candidates <- list(specified = list(
        tailup_initial = initial_object$tailup_initial,
        taildown_initial = initial_object$taildown_initial,
        euclid_initial = initial_object$euclid_initial,
        nugget_initial = initial_object$nugget_initial,
        randcov_initial = randcov_initial_obj
      ))
    } else {
      candidates <- get_decorrelate_grid(
        formula, ssn.object, initial_object, additive, anisotropy, random,
        randcov_initial_obj, dense_grid, add_iid = FALSE
      )
    }
  } else {
    candidates <- get_decorrelate_grid_candidates(grid)
  }
  if (automatic_grid) candidates <- add_decorrelate_iid(candidates, random)

  if (missing(evaluate_test)) {
    evaluate_test <- !is.null(grid) || length(candidates) > 1L || automatic_grid
  }
  if (!is.logical(evaluate_test) || length(evaluate_test) != 1L || is.na(evaluate_test)) {
    stop("evaluate_test must be TRUE or FALSE.", call. = FALSE)
  }
  if (automatic_grid) evaluate_test <- TRUE

  if (evaluate_test) {
    response_index <- get_decorrelate_response_index(formula, ssn.object)
    training <- get_decorrelate_training(
      training, response_index, NROW(ssn.object$obs)
    )
    evaluation <- get_decorrelate_grid_evaluation(
      candidates, training, formula, ssn.object,
      tailup_type, taildown_type, euclid_type, nugget_type,
      additive, anisotropy, random, partition_factor, ordering,
      algorithm, statistic, local, list(...)
    )
    summary_grid <- get_decorrelate_grid_summary(evaluation, statistic)
    parameter_grid <- get_decorrelate_grid_parameters(candidates)
    summary_grid <- merge(parameter_grid, summary_grid, by = "candidate", sort = FALSE)
    if (identical(statistic, "cor2")) {
      summary_grid <- summary_grid[order(summary_grid[[statistic]], decreasing = TRUE), , drop = FALSE]
    } else if (identical(statistic, "bias")) {
      summary_grid <- summary_grid[order(abs(summary_grid[[statistic]])), , drop = FALSE]
    } else {
      summary_grid <- summary_grid[order(summary_grid[[statistic]]), , drop = FALSE]
    }
    rownames(summary_grid) <- NULL
    attr(summary_grid, "statistic") <- statistic
    best <- summary_grid$candidate[[1]]
    candidates <- candidates[best]
    test <- as.list(summary_grid[1, c("bias", "MSPE", "RMSPE", "cor2"), drop = FALSE])
    test$statistic <- statistic
    summary_grid$candidate <- NULL
    grid_result <- structure(summary_grid, class = c("ssn_decorrelate_grid", "data.frame"))
  } else {
    training <- NULL
    test <- NULL
    grid_result <- NULL
  }

  candidate <- candidates[[1]]
  data_args <- get_decorrelate_data_args(
    formula, ssn.object,
    candidate$tailup_initial, candidate$taildown_initial,
    candidate$euclid_initial, candidate$nugget_initial,
    additive, candidate$randcov_initial, partition_factor
  )
  data_args$ordering <- ordering
  data_args$local <- local
  decorrelate_data <- do.call(ssn_decorrelate_data, data_args)
  fit <- fit_decorrelate_algorithm(decorrelate_data$tX, decorrelate_data$ty, algorithm, list(...))

  structure(
    list(
      algorithm = algorithm, call = match.call(), decorrelate_data = decorrelate_data,
      fit = fit, grid = grid_result, newdata = if (".missing" %in% names(decorrelate_data$covariance_fit$ssn.object$preds)) ".missing" else NULL,
      test = test, training = training
    ),
    class = "ssn_decorrelate"
  )
}

#' Apply the Spatial Decorrelation Transformation to a Data Object
#'
#' @description Apply the spatial decorrelation transformation to a data object.
#'   This object contains the transformed explanatory and response variables
#'   which can be used to fit a machine learning model. This object also contains
#'   information needed to decorrelate prediction data.
#' @inheritParams ssn_decorrelate
#' @param tailup_params A [tailup_params()] object giving the known tailup
#'   covariance parameters, or omitted/`NULL` for no tailup covariance
#'   (equivalent to `tailup_type = "none"`).
#' @param taildown_params A [taildown_params()] object giving the known
#'   taildown covariance parameters, or omitted/`NULL` for no taildown
#'   covariance (equivalent to `taildown_type = "none"`).
#' @param euclid_params A [euclid_params()] object giving the known Euclidean
#'   covariance parameters, or omitted/`NULL` for no Euclidean covariance
#'   (equivalent to `euclid_type = "none"`). Anisotropy is modeled
#'   automatically whenever \code{euclid_params}'s \code{rotate} is nonzero
#'   or its \code{scale} is not one; there is no separate \code{anisotropy}
#'   argument.
#' @param nugget_params A [nugget_params()] object giving the known nugget
#'   covariance parameter. Required (there is no default): every active
#'   covariance parameter, including the nugget, must be supplied and known.
#' @param randcov_params A [spmodel::randcov_params()] object giving the
#'   known random-effect variances, or omitted/`NULL` for no random effects.
#'   The random effect grouping structure is derived automatically from its
#'   names; there is no separate \code{random} argument.
#' @details The spatial decorrelation transformation is a preprocessing
#'   transformation that reduces the impacts of spatial dependence (i.e.,
#'   covariance, correlation) on machine learning models. See
#'   [ssn_decorrelate()] and Heaton et al., 2025 for more details.
#' @return A list with many elements that store information about
#'   the fitted model object. Importantly, the list contains the following elements:
#'   \itemize{
#'     \item \code{X}: The original fixed effects design matrix (of explanatory variables)
#'     \item \code{y}: The original response variable
#'     \item \code{tX}: The spatially decorrelated transformed fixed effects design matrix
#'     \item \code{ty}: The spatially decorrelated transformed response variable
#'   }
#' @seealso [ssn_decorrelate()], [ssn_decorrelate_newdata()],
#'   [ssn_recorrelate_newdata()]
#' @references Matthew J. Heaton, Andrew Millane, and Jake S. Rhodes. 2025. A Scalable
#'   Spatial Decorrelation Preprocessing Approach for Machine and Deep Learning.
#'   \emph{Journal of Data Science}. 1-15, DOI 10.6339/25-JDS1210
#' @export
#' @examples
#' \donttest{
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#' transformed <- ssn_decorrelate_data(
#'   Summer_mn ~ ELEV_DEM, mf04p, additive = "afvArea",
#'   tailup_params = tailup_params("exponential", de = 1, range = 10000),
#'   nugget_params = nugget_params("nugget", nugget = 0.1)
#' )
#' head(cbind(response = transformed$y, transformed = transformed$ty))
#' }
ssn_decorrelate_data <- function(formula, ssn.object,
                                 tailup_params, taildown_params,
                                 euclid_params, nugget_params,
                                 additive,
                                 randcov_params, partition_factor,
                                 ordering, local, ...) {
  if (missing(tailup_params)) tailup_params <- NULL
  if (missing(taildown_params)) taildown_params <- NULL
  if (missing(euclid_params)) euclid_params <- NULL
  if (missing(nugget_params)) nugget_params <- NULL
  if (missing(additive)) additive <- NULL
  if (missing(randcov_params)) randcov_params <- NULL
  if (missing(partition_factor)) partition_factor <- NULL
  if (is.symbol(substitute(additive))) additive <- deparse1(substitute(additive))

  # random is derived from randcov_params's own names rather than accepted
  # as a separate argument, matching spmodel's decorrelate_data()
  random <- if (is.null(randcov_params)) NULL else reformulate(names(randcov_params))

  if (missing(local)) local <- NULL
  if (missing(ordering)) ordering <- NULL
  check_decorrelate_dots(list(...))
  ordering <- get_decorrelate_ordering(ordering)

  # decorrelation requires every parameter to be known, so type is read
  # directly from each supplied *_params() object's own class
  tailup_type <- if (is.null(tailup_params)) "none" else remove_covtype(class(tailup_params))
  taildown_type <- if (is.null(taildown_params)) "none" else remove_covtype(class(taildown_params))
  euclid_type <- if (is.null(euclid_params)) "none" else remove_covtype(class(euclid_params))
  nugget_type <- if (is.null(nugget_params)) "nugget" else remove_covtype(class(nugget_params))

  pinned_tailup <- get_initial_from_params(tailup_type, tailup_params, tailup_initial)
  pinned_taildown <- get_initial_from_params(taildown_type, taildown_params, taildown_initial)
  pinned_euclid <- get_initial_from_params(euclid_type, euclid_params, euclid_initial)
  pinned_nugget <- get_initial_from_params(nugget_type, nugget_params, nugget_initial)
  randcov_initial_obj <- get_randcov_initial_from_params(randcov_params)

  initial_object <- get_initial_object(
    tailup_type, taildown_type, euclid_type, nugget_type,
    pinned_tailup, pinned_taildown, pinned_euclid, pinned_nugget
  )
  # anisotropy is inferred purely from euclid_params's known rotate/scale
  # values rather than accepted as a separate argument, matching spmodel's
  # decorrelate_data()
  anisotropy <- get_anisotropy_corrected(FALSE, initial_object)
  initial_object$tailup_initial <- tailup_initial_NA(initial_object$tailup_initial)
  initial_object$taildown_initial <- taildown_initial_NA(initial_object$taildown_initial)
  initial_object$euclid_initial <- euclid_initial_NA(initial_object$euclid_initial, list(anisotropy = anisotropy))
  initial_object$nugget_initial <- nugget_initial_NA(initial_object$nugget_initial)
  check_decorrelate_known(initial_object, randcov_initial_obj)

  context <- get_decorrelate_context(
    formula, ssn.object, initial_object, additive, anisotropy, random,
    randcov_initial_obj, partition_factor, ...
  )
  covariance_fit <- context$covariance_fit
  total_var <- get_decorrelate_total_var(covariance_fit)
  X <- context$X
  y <- context$y
  local <- get_decorrelate_local(local, length(y))
  if (any(!is.finite(y)) || any(!is.finite(X))) {
    stop("Decorrelated observed data must have finite response and fixed-effect values.", call. = FALSE)
  }

  row_index <- if (NROW(covariance_fit$ssn.object$obs) == covariance_fit$n) NULL else covariance_fit$observed_index
  observed_rows_data <- if (is.null(row_index)) covariance_fit$ssn.object$obs else covariance_fit$ssn.object$obs[row_index, , drop = FALSE]
  rows <- get_decorrelate_rows(observed_rows_data)
  if (NROW(rows) != covariance_fit$n || NROW(X) != covariance_fit$n) {
    stop("The fitted SSN row map is not aligned with the covariance matrix.", call. = FALSE)
  }
  permutation <- get_decorrelate_order(rows, ordering, observed_rows_data)
  inverse_permutation <- order(permutation)

  transformed <- get_decorrelate_covariance_graph(
    covariance_fit, X, y, permutation, local
  )

  output <- list(
    terms = covariance_fit$terms,
    xlevels = covariance_fit$xlevels,
    contrasts = covariance_fit$contrasts,
    rows = rows,
    X = X,
    y = y,
    tX = transformed$tX[inverse_permutation, , drop = FALSE],
    ty = as.numeric(transformed$ty[inverse_permutation]),
    ordering = ordering,
    local = local,
    newdata = if (".missing" %in% names(covariance_fit$ssn.object$preds)) ".missing" else NULL,
    covariance_fit = covariance_fit,
    params = covariance_fit$coefficients$params_object,
    total_var = total_var,
    random = random,
    partition_factor = partition_factor
  )
  structure(output, class = "ssn_decorrelate_data")
}

#' Apply the Spatial Decorrelation Transformation to a Newdata Object for Prediction
#'
#' @description Apply the spatial decorrelation transformation to a newdata object.
#'   This object contains explanatory variables that are transformed for prediction
#'   according to some spatial decorrelation transformation.
#' @param object An [ssn_decorrelate_data()] object.
#' @param newdata A character vector that indicates the name of the prediction data set
#'   for which predictions are desired (accessible via \code{object$ssn.object$preds}).
#'   Note that the prediction data must be in the original SSN object used to fit the model.
#'   If \code{newdata} is omitted, predictions
#'   for all prediction data sets are returned. Note that the name \code{".missing"}
#'   indicates the prediction data set that contains the missing observations in the data used
#'   to fit the model.
#' @param local An optional logical or list controlling the big data approximation.
#'   If omitted, \code{local} is set
#'   to \code{TRUE} or \code{FALSE} based on the sample size (the number of
#'   non-missing observations in \code{data}) -- if the sample size exceeds 5,000,
#'   \code{local} is set to \code{TRUE}. Otherwise it is set to \code{FALSE}.
#'   If \code{local} is \code{FALSE}, no big data approximation
#'   is implemented. If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{method}: The big data approximation method. If \code{method = "all"},
#'       all observations are used and \code{size} is ignored. 
#'       If \code{method = "covariance"}, the \code{size} data observations
#'       with the highest covariance with the observation requiring prediction are used.
#'       The default is \code{"covariance"}.
#'     \item \code{size}: The number of data observations to use when \code{method}
#'       is \code{"distance"} or \code{"covariance"}. The default is 30.
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. This can significantly speed
#'       up computations even when \code{method = "all"} (i.e., no big data
#'       approximation is used), as predictions
#'       are spread out over multiple cores. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 30, method = "covariance", parallel = FALSE)}.
#' @param ... Other arguments.
#'
#' @return A list with many elements that store information about
#'   the fitted model object. Importantly, the list contains the following element:
#'   \itemize{
#'     \item \code{X_newdata}: The original fixed effects design matrix (of explanatory variables) for the prediction data.
#'     \item \code{tX_newdata}: The spatially decorrelated fixed effects design matrix for the prediction data.
#'   }
#' @export
#' @examples
#' \donttest{
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path,
#'   predpts = "CapeHorn", overwrite = TRUE
#' )
#' ssn_create_distmat(mf04p, predpts = "CapeHorn", overwrite = TRUE)
#' transformed <- ssn_decorrelate_data(
#'   Summer_mn ~ ELEV_DEM, mf04p, additive = "afvArea",
#'   tailup_params = tailup_params("exponential", de = 1, range = 10000),
#'   nugget_params = nugget_params("nugget", nugget = 0.1)
#' )
#' new_transformed <- ssn_decorrelate_newdata(transformed, "CapeHorn")
#' head(new_transformed$tX_newdata)
#' }
ssn_decorrelate_newdata <- function(object, newdata, local, ...) {
  if (!inherits(object, "ssn_decorrelate_data")) {
    stop("object must be returned by ssn_decorrelate_data().", call. = FALSE)
  }
  if (missing(newdata)) newdata <- object$newdata
  if (!is.character(newdata) || length(newdata) != 1L || is.na(newdata)) {
    stop("newdata must be one named SSN prediction set.", call. = FALSE)
  }
  if (!newdata %in% names(object$covariance_fit$ssn.object$preds)) {
    stop("newdata is not a prediction set in the decorrelation object.", call. = FALSE)
  }
  check_decorrelate_dots(list(...))
  if (missing(local) || is.null(local)) local <- object$local
  local <- get_decorrelate_local(local, length(object$y))

  prediction_data <- get_prediction_newdata(object$covariance_fit, newdata)$newdata
  matrix_object <- get_newdata_model_matrix(object$covariance_fit, prediction_data)
  X0 <- matrix_object$newdata_model
  offset0 <- matrix_object$offset
  if (is.null(offset0)) offset0 <- rep(0, NROW(X0))
  if (any(!is.finite(X0)) || any(!is.finite(offset0))) {
    stop("Decorrelated prediction data must have finite fixed-effect values and offsets.", call. = FALSE)
  }

  transformed <- get_decorrelate_covariance_prediction(
    object, newdata, prediction_data, X0, local
  )
  tX0 <- transformed$tX_newdata
  yscale <- transformed$yscale
  yoffset <- transformed$yoffset

  output <- list(
    rows = get_decorrelate_rows(prediction_data),
    X_newdata = X0,
    offset = offset0,
    tX_newdata = tX0,
    yscale = yscale,
    yoffset = yoffset,
    local = local
  )
  structure(output, class = "ssn_decorrelate_newdata")
}

#' Recorrelate Machine Learning Predictions
#'
#' @description Recorrelate machine learning predictions according to a
#'   spatial decorrelation transformation.
#'
#'
#' @param object An [ssn_decorrelate_newdata()] object.
#' @param ty_newdata Predictions from the machine learning model trained on the
#'   spatially decorrelated data and applied to spatially decorrelated newdata.
#'
#' @return A vector of predictions
#' @export
#' @examples
#' \donttest{
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path,
#'   predpts = "CapeHorn", overwrite = TRUE
#' )
#' ssn_create_distmat(mf04p, predpts = "CapeHorn", overwrite = TRUE)
#' transformed <- ssn_decorrelate_data(
#'   Summer_mn ~ ELEV_DEM, mf04p, additive = "afvArea",
#'   tailup_params = tailup_params("exponential", de = 1, range = 10000),
#'   nugget_params = nugget_params("nugget", nugget = 0.1)
#' )
#' learner <- lm.fit(transformed$tX, transformed$ty)
#' new_transformed <- ssn_decorrelate_newdata(transformed, "CapeHorn")
#' predictions <- ssn_recorrelate_newdata(
#'   new_transformed, as.numeric(new_transformed$tX_newdata %*% learner$coefficients)
#' )
#' head(predictions)
#' }
ssn_recorrelate_newdata <- function(object, ty_newdata) {
  if (!inherits(object, "ssn_decorrelate_newdata")) {
    stop("object must be returned by ssn_decorrelate_newdata().", call. = FALSE)
  }
  if (missing(ty_newdata) || !is.numeric(ty_newdata) || length(ty_newdata) != NROW(object$tX_newdata) ||
      any(!is.finite(ty_newdata))) {
    stop("ty_newdata must contain one finite decorrelated value per prediction row.", call. = FALSE)
  }
  output <- object$yscale * as.numeric(ty_newdata) + object$yoffset
  if (!is.null(object$offset)) {
    output <- output + object$offset
  }
  output
}

#' Create a Spatial Decorrelation Transformation Grid
#'
#' @description
#'  Create a spatial decorrelation transformation grid of initial parameters to be
#'   evaluated via a grid search.
#'
#' @inheritParams ssn_decorrelate
#'
#' @return A grid of spatial decorrelation parameters stored as a \code{data.frame}.
#' @export
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#' grid <- ssn_decorrelate_grid(
#'   Summer_mn ~ ELEV_DEM, mf04p, tailup_type = "exponential",
#'   additive = "afvArea", dense_grid = FALSE
#' )
#' tidy(grid)
#' # Supply edited rows as grid = custom_grid to ssn_decorrelate().
#' custom_grid <- grid[1:2, ]
#' custom_grid$tailup_range <- c(5000, 10000)
ssn_decorrelate_grid <- function(formula, ssn.object,
                                  tailup_type = "none", taildown_type = "none",
                                  euclid_type = "none", nugget_type = "nugget",
                                  tailup_params, taildown_params,
                                  euclid_params, nugget_params,
                                  additive, anisotropy = FALSE, random,
                                  randcov_params, dense_grid) {
  # a supplied *_params() object's own class is authoritative for type
  # (matching spmodel's decorrelate()) unless the caller names *_type
  # explicitly -- resolved here, in the direct calling frame, since
  # missing() can no longer distinguish "defaulted" from "explicitly
  # supplied" once these are forwarded to ssn_decorrelate_grid_internal()
  if (missing(tailup_type) && !missing(tailup_params)) tailup_type <- remove_covtype(class(tailup_params))
  if (missing(taildown_type) && !missing(taildown_params)) taildown_type <- remove_covtype(class(taildown_params))
  if (missing(euclid_type) && !missing(euclid_params)) euclid_type <- remove_covtype(class(euclid_params))
  if (missing(nugget_type) && !missing(nugget_params)) nugget_type <- remove_covtype(class(nugget_params))
  # additive's NSE resolution (bare column-name symbol vs. already-quoted
  # string) is only meaningful at this, the direct calling frame --
  # ssn_decorrelate_grid_internal() receives the already-resolved value and
  # must not re-substitute
  if (missing(additive)) {
    additive <- NULL
  } else if (is.symbol(substitute(additive))) {
    additive <- deparse1(substitute(additive))
  }
  if (missing(dense_grid)) dense_grid <- FALSE
  ssn_decorrelate_grid_internal(
    formula, ssn.object, tailup_type, taildown_type, euclid_type, nugget_type,
    tailup_params, taildown_params, euclid_params, nugget_params,
    additive, anisotropy, random, randcov_params, dense_grid,
    add_iid = TRUE, candidates = NULL
  )
}

#' Shared worker behind \code{ssn_decorrelate_grid()} and the automatic grid search
#'
#' Called by \code{\link{ssn_decorrelate_grid}()} (with \code{add_iid = TRUE},
#' \code{candidates = NULL}) and by \code{\link{ssn_decorrelate}()}'s
#' automatic-grid-construction path. Also accepts an already-known
#' \code{candidates} list or an existing grid \code{data.frame} (validated and
#' re-tidied via \code{\link{get_decorrelate_grid_table}()}) instead of
#' constructing one from \code{formula}/\code{ssn.object}, mirroring how
#' \code{\link{ssn_decorrelate}()} accepts a supplied \code{grid}.
#'
#' @param formula,ssn.object,tailup_type,taildown_type,euclid_type,nugget_type,tailup_params,taildown_params,euclid_params,nugget_params,additive,anisotropy,random,randcov_params,dense_grid
#'   See \code{\link{ssn_decorrelate_grid}()}. \code{additive} must already be
#'   resolved to a string or \code{NULL} (not a bare symbol).
#' @param add_iid Whether to append an untransformed baseline candidate.
#' @param candidates An already-known named list of covariance candidates, or
#'   an existing grid \code{data.frame}, used instead of constructing a grid
#'   from \code{formula}/\code{ssn.object} when non-\code{NULL}.
#'
#' @return A data frame of class \code{ssn_decorrelate_grid}.
#'
#' @noRd
ssn_decorrelate_grid_internal <- function(formula, ssn.object,
                                           tailup_type = "none", taildown_type = "none",
                                           euclid_type = "none", nugget_type = "nugget",
                                           tailup_params, taildown_params,
                                           euclid_params, nugget_params,
                                           additive, anisotropy = FALSE, random,
                                           randcov_params, dense_grid = FALSE,
                                           add_iid, candidates) {
  # a supplied *_params() object's own class is authoritative for type
  # (matching spmodel's decorrelate()) unless the caller names *_type
  # explicitly
  if (missing(tailup_type) && !missing(tailup_params)) tailup_type <- remove_covtype(class(tailup_params))
  if (missing(taildown_type) && !missing(taildown_params)) taildown_type <- remove_covtype(class(taildown_params))
  if (missing(euclid_type) && !missing(euclid_params)) euclid_type <- remove_covtype(class(euclid_params))
  if (missing(nugget_type) && !missing(nugget_params)) nugget_type <- remove_covtype(class(nugget_params))
  if (missing(tailup_params)) tailup_params <- NULL
  if (missing(taildown_params)) taildown_params <- NULL
  if (missing(euclid_params)) euclid_params <- NULL
  if (missing(nugget_params)) nugget_params <- NULL
  if (missing(additive)) additive <- NULL
  if (missing(random)) random <- NULL
  if (missing(randcov_params)) randcov_params <- NULL
  if (missing(candidates)) candidates <- NULL
  if (!missing(formula) && is.list(formula) && missing(ssn.object)) candidates <- formula
  if (is.null(candidates)) {
    pinned_tailup <- get_initial_from_params(tailup_type, tailup_params, tailup_initial)
    pinned_taildown <- get_initial_from_params(taildown_type, taildown_params, taildown_initial)
    pinned_euclid <- get_initial_from_params(euclid_type, euclid_params, euclid_initial)
    pinned_nugget <- get_initial_from_params(nugget_type, nugget_params, nugget_initial)
    randcov_initial_obj <- get_randcov_initial_from_params(randcov_params)
    initial_object <- get_initial_object(
      tailup_type, taildown_type, euclid_type, nugget_type,
      pinned_tailup, pinned_taildown, pinned_euclid, pinned_nugget
    )
    candidates <- get_decorrelate_grid(
      formula, ssn.object, initial_object, additive, anisotropy, random,
      randcov_initial_obj, dense_grid, add_iid
    )
  }
  get_decorrelate_grid_table(candidates)
}

#' Validate/tidy a candidates list or grid \code{data.frame} into a printable grid
#'
#' Shared tail of \code{\link{ssn_decorrelate_grid_internal}()}: converts an
#' already-known \code{candidates} list (or re-validates an existing grid
#' \code{data.frame}, via \code{\link{get_decorrelate_grid_candidates}()})
#' into the tidy, printable \code{ssn_decorrelate_grid} table returned by
#' \code{\link{ssn_decorrelate_grid}()}.
#'
#' @param candidates A named list of covariance candidates, or a grid
#'   \code{data.frame} as returned by \code{\link{ssn_decorrelate_grid}()}.
#'
#' @return A data frame of class \code{ssn_decorrelate_grid}.
#'
#' @noRd
get_decorrelate_grid_table <- function(candidates) {
  candidates <- get_decorrelate_grid_candidates(candidates)
  grid <- get_decorrelate_grid_parameters(candidates)
  grid$candidate <- NULL
  structure(grid, class = c("ssn_decorrelate_grid", "data.frame"))
}

#' @rdname predict.SSN2
#' @export
#' @method predict ssn_decorrelate
predict.ssn_decorrelate <- function(object, newdata, local, ...) {
  if (!inherits(object, "ssn_decorrelate")) {
    stop("object must be returned by ssn_decorrelate().", call. = FALSE)
  }
  if (missing(newdata)) newdata <- object$newdata
  if (missing(local)) local <- NULL
  check_decorrelate_dots(list(...))
  if (is.null(newdata) || !length(newdata)) {
    stop("No prediction set was supplied and the fitted object has no missing responses.", call. = FALSE)
  }
  if (identical(newdata, "all")) newdata <- names(object$decorrelate_data$covariance_fit$ssn.object$preds)
  if (!is.character(newdata) || any(!newdata %in% names(object$decorrelate_data$covariance_fit$ssn.object$preds))) {
    stop("newdata must name one or more prediction sets in the fitted SSN object.", call. = FALSE)
  }
  predict_one <- function(name) {
    transformed <- ssn_decorrelate_newdata(
      object$decorrelate_data, name, local = local
    )
    predicted <- predict_decorrelate_algorithm(object$fit, transformed$tX_newdata, object$algorithm)
    ssn_recorrelate_newdata(transformed, predicted)
  }
  result <- lapply(newdata, predict_one)
  names(result) <- newdata
  if (length(result) == 1L) result[[1]] else result
}

#' @rdname print.SSN2
#' @export
#' @method print ssn_decorrelate
print.ssn_decorrelate <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("\n")
  if (!is.null(x$test)) {
    statistics <- unlist(x$test[c("bias", "MSPE", "RMSPE", "cor2")])
    cat("stats:\n")
    print.default(format(statistics, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  params <- x$decorrelate_data$params
  for (component in c("tailup", "taildown", "euclid")) {
    coefficients <- params[[component]]
    type <- remove_covtype(class(coefficients))
    if (identical(type, "none")) next
    if (component == "euclid" && !x$decorrelate_data$covariance_fit$anisotropy) {
      coefficients <- coefficients[!names(coefficients) %in% c("rotate", "scale")]
    }
    label <- if (component == "euclid") "Euclidean" else component
    cat("\nCoefficients (", type, " ", label, " covariance):\n", sep = "")
    print.default(format(coefficients, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  if (!inherits(params$nugget, "nugget_none")) {
    cat("\nCoefficients (nugget):\n")
    print.default(format(params$nugget, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  if (length(params$randcov)) {
    cat("Coefficients (random effects):\n")
    print.default(format(params$randcov, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  invisible(x)
}

#' @rdname ssn_decorrelate_grid
#' @export
#' @method print ssn_decorrelate_grid
print.ssn_decorrelate_grid <- function(x, ...) {
  values <- if (is.data.frame(x)) x else get_decorrelate_grid_parameters(x)
  values$candidate <- NULL
  print.data.frame(values, ...)
  invisible(x)
}

#' @rdname ssn_decorrelate
#' @export
#' @method tidy ssn_decorrelate
tidy.ssn_decorrelate <- function(x, ...) {
  if (is.null(x$grid)) return(tibble::tibble())
  tidy.ssn_decorrelate_grid(x$grid, ...)
}

#' @rdname ssn_decorrelate_grid
#' @param x An object from \code{object$grid}.
#' @param sort_by Sort by a specific row in \code{x}. Fit statistics are
#'   \code{"bias"}, \code{"MSPE"}, \code{"RMSPE"}, and \code{"cor2"}.
#'   The default is \code{"MSPE"}.
#' @param decreasing Whether \code{sort_by} should sort in decreasing order?
#'   If \code{sort_by = "cor2"}, the default is \code{TRUE}; otherwise it is
#'   \code{FALSE}.
#' @export
#' @method tidy ssn_decorrelate_grid
tidy.ssn_decorrelate_grid <- function(x, sort_by, decreasing, ...) {
  if (missing(sort_by)) sort_by <- attr(x, "statistic")
  if (!is.data.frame(x)) x <- get_decorrelate_grid_parameters(x)
  x$candidate <- NULL
  if (!is.null(sort_by)) {
    if (length(sort_by) != 1L || is.na(sort_by) || !sort_by %in% names(x)) {
      stop("sort_by must be a variable in x.", call. = FALSE)
    }
    if (missing(decreasing)) decreasing <- identical(sort_by, "cor2")
    if (!is.logical(decreasing) || length(decreasing) != 1L || is.na(decreasing)) {
      stop("decreasing must be TRUE or FALSE.", call. = FALSE)
    }
    x <- x[order(x[[sort_by]], decreasing = decreasing), , drop = FALSE]
  }
  types <- c("tailup_type", "taildown_type", "euclid_type", "nugget_type")
  untransformed <- x$tailup_type == "none" & x$taildown_type == "none" &
    x$euclid_type == "none" & x$nugget_nugget == 1
  random_names <- names(x)[startsWith(names(x), "randcov_")]
  for (name in random_names) untransformed <- untransformed & x[[name]] == 0
  untransformed <- which(untransformed)
  if (all(c("euclid_rotate", "euclid_scale") %in% names(x)) &&
      all((is.na(x$euclid_rotate) | x$euclid_rotate == 0) &
          (is.na(x$euclid_scale) | x$euclid_scale == 1))) {
    x$euclid_rotate <- NULL
    x$euclid_scale <- NULL
  }
  if (length(untransformed)) {
    x[untransformed, types] <- "no transformation"
    parameters <- setdiff(names(x), c(types, "bias", "MSPE", "RMSPE", "cor2"))
    x[untransformed, parameters] <- NA
  }
  tibble::as_tibble(x)
}

#' Validate a candidates list (or convert and validate a grid \code{data.frame})
#'
#' @param grid A named list of covariance candidates, or a grid
#'   \code{data.frame} (converted via
#'   \code{\link{get_decorrelate_candidates_from_table}()}).
#'
#' @return The validated named list of candidates.
#'
#' @noRd
get_decorrelate_grid_candidates <- function(grid) {
  candidates <- if (is.data.frame(grid)) get_decorrelate_candidates_from_table(grid) else grid
  if (!is.list(candidates) || !length(candidates) || is.null(names(candidates)) ||
      anyNA(names(candidates)) || any(!nzchar(names(candidates))) || anyDuplicated(names(candidates))) {
    stop("candidates must be a non-empty named list of joint covariance candidates.", call. = FALSE)
  }
  required <- c("tailup_initial", "taildown_initial", "euclid_initial", "nugget_initial")
  for (candidate_name in names(candidates)) {
    candidate <- candidates[[candidate_name]]
    if (!is.list(candidate) || !all(required %in% names(candidate))) {
      stop("Each candidate must contain tailup_initial, taildown_initial, euclid_initial, and nugget_initial.", call. = FALSE)
    }
    check_decorrelate_known(candidate[required], candidate$randcov_initial)
  }
  candidates
}

#' Validate and normalize the \code{local} argument for decorrelation
#'
#' @param local A logical or named list; see the \code{local} argument to
#'   \code{\link{ssn_decorrelate}()}/\code{\link{ssn_decorrelate_data}()}.
#' @param n The observed sample size, used to resolve \code{local = NULL}.
#'
#' @return A resolved list with \code{method} (\code{"all"} or
#'   \code{"covariance"}) and, when relevant, \code{size}.
#'
#' @noRd
get_decorrelate_local <- function(local, n) {
  if (is.null(local)) local <- n > 5000L
  if (identical(local, FALSE)) return(list(method = "all"))
  if (identical(local, TRUE)) return(list(method = "covariance", size = 30L))
  if (!is.list(local) || !length(local) || is.null(names(local)) ||
      any(!names(local) %in% c("method", "size")) || anyDuplicated(names(local))) {
    stop("local must be FALSE, TRUE, or a named list with method and/or size.", call. = FALSE)
  }
  if (is.null(local$method)) local$method <- "covariance"
  if (length(local$method) != 1L || is.na(local$method) ||
      !local$method %in% c("all", "covariance")) {
    stop("local$method must be \"all\" or \"covariance\".", call. = FALSE)
  }
  if (identical(local$method, "all")) return(list(method = "all"))
  if (identical(local$method, "covariance")) {
    if (is.null(local$size)) local$size <- 30L
    if (length(local$size) != 1L || !is.numeric(local$size) || !is.finite(local$size) ||
        local$size < 1 || local$size != floor(local$size) || local$size > .Machine$integer.max) {
      stop("local$size must be one positive integer.", call. = FALSE)
    }
    local$size <- as.integer(local$size)
  }
  local[c("method", "size")]
}

#' Validate and normalize the \code{ordering} argument for decorrelation
#'
#' @param ordering A single string or \code{NULL}; see the \code{ordering}
#'   argument to \code{\link{ssn_decorrelate}()}/\code{\link{ssn_decorrelate_data}()}.
#'
#' @return The resolved ordering string, one of \code{"pid"}, \code{"none"},
#'   \code{"random"}, \code{"maxmin"}, \code{"middleout"}, \code{"outsidein"},
#'   \code{"coordinate"}, or \code{"grts"}.
#'
#' @noRd
get_decorrelate_ordering <- function(ordering) {
  if (is.null(ordering)) ordering <- "pid"
  if (length(ordering) != 1L || is.na(ordering) ||
      !ordering %in% c(
        "pid", "none", "random", "maxmin", "middleout", "outsidein",
        "coordinate", "grts"
      )) {
    stop(
      "ordering must be \"pid\", \"none\", \"random\", \"maxmin\", \"middleout\", ",
      "\"outsidein\", \"coordinate\", or \"grts\".",
      call. = FALSE
    )
  }
  ordering
}

#' Reject renamed/removed decorrelation arguments passed through \code{...}
#'
#' @param dots The result of \code{list(...)} from a decorrelation function.
#'
#' @return \code{NULL}, invisibly, if \code{dots} contains none of the
#'   rejected names; otherwise an error naming the current equivalent.
#'
#' @noRd
check_decorrelate_dots <- function(dots) {
  if ("conditioning" %in% names(dots)) {
    stop("conditioning is not an argument; use local = list(method = ..., size = ...).", call. = FALSE)
  }
  if ("covariance_args" %in% names(dots)) {
    stop("covariance_args is not an argument; supply grid or initial values.", call. = FALSE)
  }
}

#' Check that every covariance parameter is supplied and known
#'
#' Decorrelation requires every active covariance parameter (and, if random
#' effects are used, every random-effect variance) to be fully known; this
#' rejects any \code{initial_object}/\code{randcov_initial} that still has an
#' unknown or missing field.
#'
#' @param initial_object A joint covariance initial-value object with
#'   \code{tailup_initial}, \code{taildown_initial}, \code{euclid_initial},
#'   and \code{nugget_initial}.
#' @param randcov_initial A random-effect variance initial-value object, or
#'   \code{NULL}.
#'
#' @return \code{NULL}, invisibly, if every parameter is known; otherwise an
#'   error.
#'
#' @noRd
check_decorrelate_known <- function(initial_object, randcov_initial = NULL) {
  component_names <- c("tailup_initial", "taildown_initial", "euclid_initial", "nugget_initial")
  for (component_name in component_names) {
    component <- initial_object[[component_name]]
    if (is.null(component)) {
      stop("All covariance components must have fully known initial values for decorrelation.", call. = FALSE)
    }
    is_none <- any(grepl("_none$", class(component)))
    if (!is_none && (!length(component$initial) || !length(component$is_known) ||
        !all(component$is_known) || any(!is.finite(component$initial)))) {
      stop("All covariance parameters must be supplied and known for decorrelation.", call. = FALSE)
    }
  }
  if (is.null(randcov_initial) && !is.null(initial_object$randcov_initial)) {
    stop("All random-effect covariance parameters must be supplied and known for decorrelation.", call. = FALSE)
  }
  if (!is.null(randcov_initial) &&
      (!length(randcov_initial$initial) || !length(randcov_initial$is_known) ||
       !all(randcov_initial$is_known) || any(!is.finite(randcov_initial$initial)))) {
    stop("All random-effect covariance parameters must be supplied and known for decorrelation.", call. = FALSE)
  }
  invisible(NULL)
}

#' Build the row map used by decorrelation/simulation ordering and diagnostics
#'
#' @param data An SSN observed or prediction data frame (with \code{netgeom}
#'   information).
#' @param index An optional integer index selecting/reordering \code{data}'s
#'   rows first, or \code{NULL} to use every row in place.
#'
#' @return A data frame with \code{row} (position in \code{data}, after
#'   \code{index}), \code{NetworkID}, and \code{pid}.
#'
#' @noRd
get_decorrelate_rows <- function(data, index = NULL) {
  if (!is.null(index)) data <- data[index, , drop = FALSE]
  netgeom <- ssn_get_netgeom(data, reformat = TRUE)
  network_id <- as.character(netgeom$NetworkID)
  pid <- as.integer(netgeom$pid)
  data.frame(
    row = seq_len(NROW(data)), NetworkID = network_id, pid = pid,
    stringsAsFactors = FALSE
  )
}

#' Order rows for sequential decorrelation/simulation conditioning
#'
#' \code{"none"} preserves row order, \code{"random"} shuffles rows,
#' \code{"pid"} orders deterministically by point ID then \code{NetworkID}
#' (breaking remaining ties by row position). \code{"maxmin"},
#' \code{"middleout"}, \code{"outsidein"}, and \code{"coordinate"} (via
#' \code{GPvecchia}) and \code{"grts"} (via \code{spsurvey}) order by the
#' observations' own coordinates, matching spmodel's
#' \code{get_decorrelate_order()} (\code{R/decorrelate_data.R}) -- these are
#' numerical conditioning orders only; they say nothing about stream
#' topology or hydrologic flow direction.
#'
#' @param rows A row map data frame from \code{\link{get_decorrelate_rows}()}.
#' @param ordering One of \code{"none"}, \code{"random"}, \code{"pid"},
#'   \code{"maxmin"}, \code{"middleout"}, \code{"outsidein"},
#'   \code{"coordinate"}, or \code{"grts"}.
#' @param data The \code{sf} rows \code{rows} was built from. Only used (and
#'   only required) for the coordinate-based ordering values.
#'
#' @return An integer permutation of \code{seq_len(NROW(rows))}.
#'
#' @noRd
get_decorrelate_order <- function(rows, ordering, data = NULL) {
  if (identical(ordering, "none")) return(seq_len(NROW(rows)))
  if (identical(ordering, "random")) return(sample.int(NROW(rows)))
  if (identical(ordering, "pid")) return(order(rows$pid, rows$NetworkID, rows$row))

  n <- NROW(rows)
  coords <- sf::st_coordinates(data)
  x <- coords[, 1]
  y <- coords[, 2]

  if (identical(ordering, "maxmin")) {
    if (!requireNamespace("GPvecchia", quietly = TRUE)) {
      stop("Install the GPvecchia package before using \"maxmin\" ordering.", call. = FALSE)
    }
    return(GPvecchia::order_maxmin_exact(cbind(x, y)))
  }
  if (identical(ordering, "middleout")) {
    if (!requireNamespace("GPvecchia", quietly = TRUE)) {
      stop("Install the GPvecchia package before using \"middleout\" ordering.", call. = FALSE)
    }
    return(GPvecchia::order_middleout(cbind(x, y)))
  }
  if (identical(ordering, "outsidein")) {
    if (!requireNamespace("GPvecchia", quietly = TRUE)) {
      stop("Install the GPvecchia package before using \"outsidein\" ordering.", call. = FALSE)
    }
    return(GPvecchia::order_outsidein(cbind(x, y)))
  }
  if (identical(ordering, "coordinate")) {
    if (!requireNamespace("GPvecchia", quietly = TRUE)) {
      stop("Install the GPvecchia package before using \"coordinate\" ordering.", call. = FALSE)
    }
    return(GPvecchia::order_coordinate(cbind(x, y), coordinate = c(1, 2)))
  }

  # ordering == "grts"
  if (!requireNamespace("spsurvey", quietly = TRUE)) {
    stop("Install the spsurvey package before using \"grts\" ordering.", call. = FALSE)
  }
  frame <- sf::st_as_sf(
    data.frame(x = x, y = y, .index. = seq_len(n)),
    coords = c("x", "y"), crs = NA
  )
  samp <- spsurvey::grts(frame, n_base = n, projcrs_check = FALSE)
  samp$sites_base$.index.
}

#' Build the observed-data covariance context shared by decorrelation and simulation
#'
#' Assembles the observed design matrix/response (offset-subtracted),
#' resolves random effects/partition factor/anisotropy, and packages the
#' known covariance parameters into an \code{ssn_lm}-classed
#' \code{covariance_fit} object -- the shape the shared
#' \code{get_decorrelate_*} covariance-assembly helpers expect (matching
#' \code{covmatrix()}'s own fitted-model input shape).
#'
#' @param formula A model formula.
#' @param ssn.object A fitted-model-ready SSN object.
#' @param initial_object A joint covariance initial-value object with fully
#'   known parameters.
#' @param additive The additive function value column name, or \code{NULL}.
#' @param anisotropy Whether Euclidean anisotropy is active.
#' @param random A one- or two-sided random effect formula, or \code{NULL}.
#' @param randcov_initial A random-effect variance initial-value object, or
#'   \code{NULL}.
#' @param partition_factor A one-sided partition factor formula, or
#'   \code{NULL}.
#' @param ... Additional arguments passed to
#'   \code{\link{get_model_matrix_object}()}, such as \code{contrasts}.
#'
#' @return A list with \code{covariance_fit} (the \code{ssn_lm}-classed
#'   context), \code{X}, and \code{y} (offset-subtracted).
#'
#' @noRd
get_decorrelate_context <- function(formula, ssn.object, initial_object,
                                          additive, anisotropy, random,
                                          randcov_initial, partition_factor, ...) {
  response_name <- all.vars(formula)[[1]]
  observed_index <- which(!is.na(ssn.object$obs[[response_name]]))
  missing_index <- which(is.na(ssn.object$obs[[response_name]]))
  observed <- ssn.object$obs[observed_index, , drop = FALSE]
  sf_column_name <- attr(observed, "sf_column")
  matrix_object <- get_model_matrix_object(
    formula, observed, sf_column_name = sf_column_name, ...
  )
  observed <- matrix_object$obdata
  X <- matrix_object$X
  response <- as.numeric(model.response(matrix_object$obdata_model_frame))
  offset <- model.offset(matrix_object$obdata_model_frame)
  if (is.null(offset)) offset <- rep(0, length(response))
  y <- response - offset
  if (length(observed_index) != NROW(observed)) {
    stop("Decorrelated data cannot have missing fixed-effect values.", call. = FALSE)
  }

  partition_factor <- coerce_partition_factor(partition_factor, observed)
  partition_xlev <- get_partition_xlev(partition_factor, observed)
  random_context <- build_randcov_list(
    random, randcov_initial, observed, rep.int(1L, NROW(observed))
  )
  initial_object$randcov_initial <- random_context$randcov_initial
  params <- get_params_object_known(initial_object)
  anisotropy <- get_anisotropy_corrected(anisotropy, initial_object)
  ssn_restructured <- restruct_ssn_missing(ssn.object, observed_index, missing_index)

  covariance_fit <- structure(
    list(
      additive = additive,
      anisotropy = anisotropy,
      coefficients = list(params_object = params),
      contrasts = matrix_object$dots$contrasts,
      diagtol = 0,
      formula = matrix_object$formula,
      n = NROW(observed),
      observed_index = observed_index,
      missing_index = missing_index,
      partition_factor = partition_factor,
      partition_xlev = partition_xlev,
      random = random,
      random_xlev = random_context$randcov_xlev,
      ssn.object = ssn_restructured,
      terms = matrix_object$terms_val,
      xlevels = matrix_object$xlevels
    ),
    class = "ssn_lm"
  )
  list(covariance_fit = covariance_fit, X = X, y = y)
}

#' Apply the sequential covariance-neighbor decorrelation transform to observed data
#'
#' The training-time Vecchia-style engine behind \code{\link{ssn_decorrelate_data}()}'s
#' \code{local} (big-data) path: each ordered observation is standardized by
#' its conditional distribution given (up to) \code{local$size} of the
#' most-correlated earlier-ordered observations, via
#' \code{\link{get_decorrelate_conditional_values}()}.
#'
#' @param covariance_fit A \code{covariance_fit}-shaped list from
#'   \code{\link{get_decorrelate_context}()}.
#' @param X The observed design matrix.
#' @param y The observed (offset-subtracted) response.
#' @param permutation The row ordering to condition sequentially along, from
#'   \code{\link{get_decorrelate_order}()}.
#' @param local A resolved \code{local} list from
#'   \code{\link{get_decorrelate_local}()}.
#'
#' @return A list with \code{tX} and \code{ty} (in \code{permutation} order).
#'
#' @noRd
get_decorrelate_covariance_graph <- function(covariance_fit, X, y, permutation,
                                             local) {
  total_var <- get_decorrelate_total_var(covariance_fit)
  n <- length(permutation)
  tX <- matrix(NA_real_, nrow = n, ncol = NCOL(X), dimnames = list(NULL, colnames(X)))
  ty <- numeric(n)
  observed <- covariance_fit$ssn.object$obs

  for (position in seq_len(n)) {
    current <- permutation[[position]]
    variance <- get_decorrelate_marginal_variance(covariance_fit, observed[current, , drop = FALSE]) / total_var
    if (position == 1L) {
      scale <- get_decorrelate_positive_scale(variance, "observed")
      tX[position, ] <- X[current, ] / scale
      ty[position] <- y[current] / scale
      next
    }

    candidates <- permutation[seq_len(position - 1L)]
    cross_covariance <- get_decorrelate_observed_cross_covariance(
      covariance_fit, observed[current, , drop = FALSE], observed[candidates, , drop = FALSE]
    ) / total_var
    # method = "all" conditions on every earlier-ordered observation (the
    # exact transform, sequentially applied -- see get_decorrelated_value()
    # in spmodel's decorrelate_data.R); method = "covariance" truncates to
    # the local$size most-correlated of them
    keep <- if (identical(local$method, "all")) {
      seq_along(candidates)
    } else {
      get_decorrelate_covariance_neighbors(cross_covariance, local$size)
    }
    neighbor_index <- candidates[keep]
    covariance_neighbors <- get_decorrelate_observed_covariance(
      covariance_fit, observed[neighbor_index, , drop = FALSE]
    ) / total_var
    solved <- get_decorrelate_conditional_values(
      covariance_neighbors, cross_covariance[keep], variance,
      X[neighbor_index, , drop = FALSE], y[neighbor_index]
    )
    tX[position, ] <- (X[current, ] - solved$mean_x) / solved$scale
    ty[position] <- (y[current] - solved$mean_y) / solved$scale
  }

  list(tX = tX, ty = ty)
}

#' Apply the sequential covariance-neighbor decorrelation transform to prediction data
#'
#' The \code{ssn_decorrelate_newdata()} analogue of
#' \code{\link{get_decorrelate_covariance_graph}()}: each \code{newdata} row
#' is standardized by its conditional distribution given (up to)
#' \code{local$size} of the most-correlated observed rows, via
#' \code{\link{get_decorrelate_conditional_values}()}, processed in chunks to
#' avoid materializing a full dense observed-by-newdata covariance matrix at
#' once.
#'
#' @param object An \code{\link{ssn_decorrelate_data}()} object.
#' @param newdata_name The name of the prediction set being transformed.
#' @param newdata The prediction data frame for \code{newdata_name}.
#' @param X0 The prediction design matrix.
#' @param local A resolved \code{local} list from
#'   \code{\link{get_decorrelate_local}()}.
#'
#' @return A list with \code{tX_newdata}, \code{yscale}, and \code{yoffset}
#'   (each length/row count \code{n}).
#'
#' @noRd
get_decorrelate_covariance_prediction <- function(object, newdata_name, newdata, X0,
                                                   local) {
  n <- NROW(newdata)
  tX <- matrix(NA_real_, nrow = n, ncol = NCOL(X0), dimnames = list(NULL, colnames(X0)))
  yscale <- numeric(n)
  yoffset <- numeric(n)
  observed <- object$covariance_fit$ssn.object$obs

  chunks <- split(seq_len(n), ceiling(seq_len(n) / min(100L, n)))
  for (chunk in chunks) {
    cross_chunk <- get_block_obs_covariance(object$covariance_fit, newdata_name, chunk) / object$total_var
    for (i in chunk) {
      cross_covariance <- cross_chunk[match(i, chunk), ]
      keep <- if (identical(local$method, "all")) {
        seq_along(cross_covariance)
      } else {
        get_decorrelate_covariance_neighbors(cross_covariance, local$size)
      }
      covariance_neighbors <- get_decorrelate_observed_covariance(
        object$covariance_fit, observed[keep, , drop = FALSE]
      ) / object$total_var
      variance <- get_decorrelate_marginal_variance(
        object$covariance_fit, newdata[i, , drop = FALSE]
      ) / object$total_var
      solved <- get_decorrelate_conditional_values(
        covariance_neighbors, cross_covariance[keep], variance,
        object$X[keep, , drop = FALSE], object$y[keep]
      )
      tX[i, ] <- (X0[i, ] - solved$mean_x) / solved$scale
      yscale[[i]] <- solved$scale
      yoffset[[i]] <- solved$mean_y
    }
  }

  list(tX_newdata = tX, yscale = yscale, yoffset = yoffset)
}

#' Select the most-correlated neighbor indices for local decorrelation/simulation
#'
#' Ranks candidates by absolute covariance (rather than raw covariance,
#' since some covariance families can be negative) and keeps the top
#' \code{size}; ties are broken in favor of later candidate positions.
#'
#' @param cross_covariance A vector of covariances between the target and
#'   each candidate.
#' @param size The maximum number of neighbors to keep.
#'
#' @return An integer index into \code{cross_covariance} selecting the kept
#'   neighbors.
#'
#' @noRd
get_decorrelate_covariance_neighbors <- function(cross_covariance, size) {
  n <- length(cross_covariance)
  if (n == 0L) return(integer())
  keep_n <- min(as.integer(size), n)
  order(abs(as.numeric(cross_covariance)))[seq.int(n, length.out = keep_n, by = -1L)]
}

#' Solve the conditional decorrelation transform given a neighbor pool
#'
#' Shared low-level solver behind \code{\link{get_decorrelate_covariance_graph}()}
#' (training) and \code{\link{get_decorrelate_covariance_prediction}()}
#' (newdata): factors the neighbor pool's covariance matrix and uses it to
#' compute the conditional mean contribution and conditional standard
#' deviation for one target row.
#'
#' @param covariance_neighbors The neighbor pool's covariance matrix.
#' @param cross_covariance The target's covariance with each neighbor.
#' @param variance The target's marginal variance.
#' @param X The neighbor pool's design matrix rows.
#' @param y The neighbor pool's response values.
#'
#' @return A list with \code{mean_x}/\code{mean_y} (the part of the target's
#'   \code{X}/\code{y} predictable from the neighbor pool) and \code{scale}
#'   (the conditional standard deviation).
#'
#' @noRd
get_decorrelate_conditional_values <- function(covariance_neighbors, cross_covariance,
                                                variance, X, y) {
  lower <- tryCatch(t(chol(covariance_neighbors)), error = function(error) NULL)
  if (is.null(lower)) {
    stop("The covariance matrix for decorrelation neighbors is not positive definite.", call. = FALSE)
  }
  q <- forwardsolve(lower, as.numeric(cross_covariance))
  conditional_variance <- variance - sum(q^2)
  scale <- get_decorrelate_positive_scale(
    conditional_variance, "conditional", reference = variance
  )
  whitened_x <- forwardsolve(lower, X)
  whitened_y <- as.numeric(forwardsolve(lower, y))
  list(
    mean_x = as.numeric(crossprod(q, whitened_x)),
    mean_y = as.numeric(crossprod(q, whitened_y)),
    scale = scale
  )
}

#' Take the square root of a conditional variance, erroring if it is not positive
#'
#' @param variance The conditional variance to take the square root of.
#' @param label A short label for the error message (e.g. \code{"observed"},
#'   \code{"conditional"}, \code{"prediction"}).
#' @param reference A reference magnitude (usually the marginal variance)
#'   used to scale the numerical tolerance.
#'
#' @return \code{sqrt(variance)}.
#'
#' @noRd
get_decorrelate_positive_scale <- function(variance, label, reference = variance) {
  tolerance <- .Machine$double.eps * max(abs(reference), abs(variance), .Machine$double.xmin)
  if (!is.finite(variance) || variance <= tolerance) {
    stop(
      "The ", label, " conditional variance is not positive; check covariance parameters and duplicate locations.",
      call. = FALSE
    )
  }
  sqrt(variance)
}

#' Compute the total (spatial + nugget + random-effect) variance
#'
#' @param covariance_fit A \code{covariance_fit}-shaped list with
#'   \code{coefficients$params_object} and \code{diagtol}.
#'
#' @return The positive total variance, used to normalize covariances into
#'   correlations before decorrelation/simulation.
#'
#' @noRd
get_decorrelate_total_var <- function(covariance_fit) {
  params <- covariance_fit$coefficients$params_object
  total_var <- sum(get_spatial_nugget_var(params, covariance_fit$diagtol), params$randcov)
  if (!is.finite(total_var) || total_var <= 0) {
    stop("Decorrelation requires a positive total variance.", call. = FALSE)
  }
  total_var
}

#' Compute one row's marginal (spatial + nugget + random-effect) variance
#'
#' @param covariance_fit A \code{covariance_fit}-shaped list with
#'   \code{coefficients$params_object} and \code{diagtol}.
#' @param data A one-row data frame for the location whose marginal variance
#'   is needed.
#'
#' @return The marginal variance at \code{data}.
#'
#' @noRd
get_decorrelate_marginal_variance <- function(covariance_fit, data) {
  params <- covariance_fit$coefficients$params_object
  get_spatial_nugget_var(params, covariance_fit$diagtol) +
    randcov_newvar(params$randcov, data)
}

#' Build the cross-covariance matrix between two arbitrary row sets
#'
#' Assembles stream/Euclidean covariance plus, if active, random-effect and
#' partition-factor contributions, between two (possibly different) row sets
#' -- the SSN-correct analogue of \code{covmatrix(object, newdata, cov_type =
#' "pred.obs")}, used by the local decorrelation/simulation neighbor-pool
#' machinery for row sets that need not both be the full observed data.
#'
#' @param covariance_fit A \code{covariance_fit}-shaped list.
#' @param d1,d2 Data frames of rows to cross.
#'
#' @return A numeric vector of length \code{NROW(d1) * NROW(d2)}, in
#'   column-major (\code{d1}-varies-fastest) order.
#'
#' @noRd
get_decorrelate_observed_cross_covariance <- function(covariance_fit, d1, d2) {
  params <- covariance_fit$coefficients$params_object
  context <- list(
    ssn.object = covariance_fit$ssn.object, additive = covariance_fit$additive,
    anisotropy = covariance_fit$anisotropy
  )
  distance <- get_dist_object_bigdata_cross(
    d1, d2, params, context, backend = get_decorrelate_observed_backend(covariance_fit)
  )
  covariance <- get_cov_matrix_cross(params, distance, data_object = context)
  target_length <- NROW(d1) * NROW(d2)
  if (length(covariance) == 1L && target_length > 1L) {
    covariance <- matrix(covariance, nrow = NROW(d1), ncol = NROW(d2))
  }
  if (!is.null(params$randcov)) {
    covariance <- covariance + randcov_vector(params$randcov, d2, d1, xlev_list = covariance_fit$random_xlev)
  }
  partition <- partition_vector(covariance_fit$partition_factor, d2, d1, xlev = covariance_fit$partition_xlev)
  if (!is.null(partition)) covariance <- covariance * partition
  as.numeric(covariance)
}

#' Build the pool covariance matrix for one row set
#'
#' The single-row-set analogue of
#' \code{\link{get_decorrelate_observed_cross_covariance}()}: assembles the
#' full covariance matrix (including the nugget on the diagonal) among
#' \code{data}'s own rows, the SSN-correct analogue of \code{covmatrix(object,
#' newdata = data, cov_type = "pred.pred")} for an arbitrary row set.
#'
#' @param covariance_fit A \code{covariance_fit}-shaped list.
#' @param data A data frame of rows to build the covariance matrix for.
#'
#' @return A \code{NROW(data) x NROW(data)} covariance matrix.
#'
#' @noRd
get_decorrelate_observed_covariance <- function(covariance_fit, data) {
  params <- covariance_fit$coefficients$params_object
  context <- list(
    ssn.object = covariance_fit$ssn.object, additive = covariance_fit$additive,
    anisotropy = covariance_fit$anisotropy
  )
  distance <- get_dist_object_bigdata_cross(
    data, data, params, context, backend = get_decorrelate_observed_backend(covariance_fit)
  )
  covariance <- as.matrix(get_cov_matrix_cross(params, distance, data_object = context))
  if (!identical(dim(covariance), rep(NROW(data), 2L)) && length(covariance) == 1L) {
    covariance <- matrix(covariance, nrow = NROW(data), ncol = NROW(data))
  }
  if (!is.null(params$randcov)) {
    # see the identical note in get_decorrelate_observed_cross_covariance():
    # data is always a subset of covariance_fit$ssn.object$obs.
    covariance <- covariance + randcov_vector(params$randcov, data, data, xlev_list = covariance_fit$random_xlev)
  }
  partition <- partition_vector(covariance_fit$partition_factor, data, data, xlev = covariance_fit$partition_xlev)
  if (!is.null(partition)) covariance <- covariance * partition
  dependent_variance <- sum(
    params$tailup[["de"]], params$taildown[["de"]], params$euclid[["de"]]
  )
  nugget_variance <- get_spatial_nugget_var(params, covariance_fit$diagtol) - dependent_variance
  diag(covariance) <- diag(covariance) + nugget_variance
  as.matrix(covariance)
}

#' Select the appropriate distance backend for decorrelation/simulation covariance
#'
#' Called by the local/Vecchia/low-rank neighbor-pool machinery (never the
#' exact/default path), which reads many small, per-pool/per-target
#' submatrices per call -- the same repeated-small-read pattern that keeps
#' local model fitting's own backend resolution \code{.bmat}-preferred (see
#' \code{get_data_object_bigdata.R}). Prefers \code{.bmat} for the same
#' reason: the dense (\code{.RData}) reader has no partial-read capability
#' and \code{unserialize()}s an entire per-network matrix on every call, with
#' no caching, so repeated small reads pay that full cost each time;
#' \code{.bmat} supports indexed reads of the requested submatrices.
#'
#' @param covariance_fit A \code{covariance_fit}-shaped list.
#'
#' @return The distance backend selected by
#'   \code{\link{select_square_dist_backend}()} for the observed data, given
#'   which stream components (if any) are active.
#'
#' @noRd
get_decorrelate_observed_backend <- function(covariance_fit) {
  params <- covariance_fit$coefficients$params_object
  select_square_dist_backend(
    covariance_fit$ssn.object, "obs",
    inherits(params$tailup, "tailup_none"), inherits(params$taildown, "taildown_none"),
    prefer = "bigdata"
  )
}

#' Convert a fully-known \code{*_initial()} object into its \code{*_params()} equivalent
#'
#' Used to pass a candidate's fixed initial values through to
#' \code{\link{ssn_decorrelate_data}()}, which requires \code{*_params()}
#' objects. The reverse of \code{\link{get_initial_from_params}()}.
#'
#' @param initial_obj A \code{*_initial()} object, or \code{NULL}.
#' @param params_fn The corresponding \code{*_params()} constructor.
#'
#' @return The equivalent \code{*_params()} object, or \code{NULL} if
#'   \code{initial_obj} is \code{NULL}.
#'
#' @noRd
get_params_from_initial <- function(initial_obj, params_fn) {
  if (is.null(initial_obj)) return(NULL)
  type <- remove_covtype(class(initial_obj))
  do.call(params_fn, c(list(type), as.list(initial_obj$initial)))
}

#' Convert a \code{*_params()} object into its \code{*_initial()} equivalent
#'
#' Used at the \code{ssn_decorrelate()}/\code{ssn_decorrelate_grid()}/
#' \code{ssn_decorrelate_data()} public boundary, where known parameters are
#' supplied as \code{*_params()} objects but the internal fitting/grid
#' machinery expects \code{*_initial()} objects with \code{known = "given"}.
#' The reverse of \code{\link{get_params_from_initial}()}.
#'
#' @param type The covariance type (e.g. \code{"exponential"},
#'   \code{"none"}).
#' @param params_obj A \code{*_params()} object, or \code{NULL}.
#' @param initial_fn The corresponding \code{*_initial()} constructor.
#'
#' @return The equivalent \code{*_initial()} object (with every field
#'   \code{known}), or \code{NULL} if \code{params_obj} is \code{NULL}.
#'
#' @noRd
get_initial_from_params <- function(type, params_obj, initial_fn) {
  if (is.null(params_obj)) return(NULL)
  do.call(initial_fn, c(list(type), as.list(params_obj), list(known = "given")))
}

#' Convert a \code{randcov_params()} object into its \code{randcov_initial()} equivalent
#'
#' The random-effect analogue of \code{\link{get_initial_from_params}()}.
#'
#' @param randcov_params A \code{\link[spmodel:randcov_params]{randcov_params()}}
#'   object, or \code{NULL}.
#'
#' @return The equivalent \code{randcov_initial()} object (with every field
#'   known), or \code{NULL} if \code{randcov_params} is \code{NULL}.
#'
#' @noRd
get_randcov_initial_from_params <- function(randcov_params) {
  if (is.null(randcov_params)) return(NULL)
  do.call(randcov_initial, c(as.list(randcov_params), list(known = "given")))
}

#' Convert a \code{randcov_initial()} object into its \code{randcov_params()} equivalent
#'
#' The random-effect analogue of \code{\link{get_params_from_initial}()}; the
#' reverse of \code{\link{get_randcov_initial_from_params}()}.
#'
#' @param randcov_initial_obj A \code{randcov_initial()} object, or
#'   \code{NULL}.
#'
#' @return The equivalent \code{randcov_params()} object, or \code{NULL} if
#'   \code{randcov_initial_obj} is \code{NULL}.
#'
#' @noRd
get_randcov_params_from_initial <- function(randcov_initial_obj) {
  if (is.null(randcov_initial_obj)) return(NULL)
  do.call(randcov_params, as.list(randcov_initial_obj$initial))
}

#' Build the argument list for a \code{ssn_decorrelate_data()} call from \code{*_initial()} objects
#'
#' Converts a candidate's \code{*_initial()} objects to \code{*_params()} via
#' \code{\link{get_params_from_initial}()}/\code{\link{get_randcov_params_from_initial}()}
#' and packages every argument \code{ssn_decorrelate_data()} needs, ready for
#' \code{do.call()}.
#'
#' @param formula A model formula.
#' @param ssn.object A fitted-model-ready SSN object.
#' @param tailup_initial,taildown_initial,euclid_initial,nugget_initial A
#'   candidate's fully known covariance initial-value objects.
#' @param additive The additive function value column name, or \code{NULL}.
#' @param randcov_initial A random-effect variance initial-value object, or
#'   \code{NULL}.
#' @param partition_factor A one-sided partition factor formula, or
#'   \code{NULL}.
#'
#' @return A named list of arguments for \code{\link{ssn_decorrelate_data}()}.
#'   \code{ssn_decorrelate_data()} infers anisotropy and \code{random} itself
#'   (from \code{euclid_params}'s rotate/scale and \code{randcov_params}'s
#'   names, respectively), so neither is included here.
#'
#' @noRd
get_decorrelate_data_args <- function(formula, ssn.object,
                                      tailup_initial, taildown_initial, euclid_initial, nugget_initial,
                                      additive, randcov_initial, partition_factor) {
  list(
    formula = formula, ssn.object = ssn.object,
    tailup_params = get_params_from_initial(tailup_initial, tailup_params),
    taildown_params = get_params_from_initial(taildown_initial, taildown_params),
    euclid_params = get_params_from_initial(euclid_initial, euclid_params),
    nugget_params = get_params_from_initial(nugget_initial, nugget_params),
    additive = additive,
    randcov_params = get_randcov_params_from_initial(randcov_initial), partition_factor = partition_factor
  )
}

#' Check whether every covariance parameter (and random effect) is already known
#'
#' Used by \code{\link{ssn_decorrelate}()} to decide whether a grid search is
#' needed at all: \code{TRUE} only if \code{\link{check_decorrelate_known}()}
#' passes and, when \code{random} is given, \code{randcov_initial} supplies
#' every one of \code{random}'s terms by name.
#'
#' @param initial_object A joint covariance initial-value object.
#' @param random A one- or two-sided random effect formula, or \code{NULL}.
#' @param randcov_initial A random-effect variance initial-value object, or
#'   \code{NULL}.
#' @param anisotropy Whether Euclidean anisotropy is active.
#'
#' @return \code{TRUE} if every covariance parameter is known, \code{FALSE}
#'   otherwise (including if validation errors).
#'
#' @noRd
get_decorrelate_covariance_known <- function(initial_object, random, randcov_initial,
                                             anisotropy = FALSE) {
  is_known <- tryCatch({
    initial_object <- get_initial_NA_object(initial_object, list(anisotropy = anisotropy))
    check_decorrelate_known(initial_object, randcov_initial)
    randcov_names <- get_randcov_names(random)
    supplied_names <- unlist(lapply(names(randcov_initial$initial), function(name) {
      get_randcov_names(reformulate(name))
    }), use.names = FALSE)
    is.null(random) || (!is.null(randcov_initial) && setequal(randcov_names, supplied_names))
  }, error = function(error) FALSE)
  isTRUE(is_known)
}

#' Fit a machine learning algorithm to decorrelated training data
#'
#' @param X The decorrelated design matrix.
#' @param y The decorrelated response.
#' @param algorithm One of \code{"ranger"}, \code{"randomForest"}, or
#'   \code{"xgboost"}.
#' @param dots Additional algorithm-specific arguments (as a list) forwarded
#'   to the fitting function.
#'
#' @return The fitted model object from the requested algorithm's package.
#'
#' @noRd
fit_decorrelate_algorithm <- function(X, y, algorithm, dots) {
  if (identical(algorithm, "ranger")) {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Install the ranger package before using algorithm = \"ranger\".", call. = FALSE)
    }
    return(do.call(ranger::ranger, c(list(x = X, y = y), dots)))
  }
  if (identical(algorithm, "randomForest")) {
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("Install the randomForest package before using algorithm = \"randomForest\".", call. = FALSE)
    }
    return(do.call(randomForest::randomForest, c(list(x = X, y = y), dots)))
  }
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Install the xgboost package before using algorithm = \"xgboost\".", call. = FALSE)
  }
  do.call(xgboost::xgboost, c(list(x = X, y = y), dots))
}

#' Predict from a fitted decorrelation-learner algorithm
#'
#' @param fit A fitted model object from \code{\link{fit_decorrelate_algorithm}()}.
#' @param X The decorrelated design matrix to predict from.
#' @param algorithm One of \code{"ranger"}, \code{"randomForest"}, or
#'   \code{"xgboost"}.
#'
#' @return A numeric vector of predictions on the decorrelated scale.
#'
#' @noRd
predict_decorrelate_algorithm <- function(fit, X, algorithm) {
  if (identical(algorithm, "ranger")) return(as.numeric(predict(fit, data = X)$predictions))
  if (identical(algorithm, "randomForest")) return(as.numeric(predict(fit, newdata = X)))
  as.numeric(predict(fit, newdata = X))
}

#' Index the non-missing (fittable) response rows for decorrelation
#'
#' @param formula A model formula.
#' @param ssn.object A fitted-model-ready SSN object.
#'
#' @return An integer index into \code{ssn.object$obs} selecting rows with a
#'   non-missing response.
#'
#' @noRd
get_decorrelate_response_index <- function(formula, ssn.object) {
  response_name <- all.vars(formula[[2]])[[1]]
  which(!is.na(ssn.object$obs[[response_name]]))
}

#' Extract a formula's response variable, keeping missing values
#'
#' Used for cross-validation error calculations, where held-out rows'
#' original (possibly \code{NA}, but here always observed) responses are
#' needed independent of the fitted model's own removal of missing rows.
#'
#' @param formula A model formula.
#' @param data A data frame containing the response variable.
#'
#' @return A numeric vector, the response evaluated over every row of
#'   \code{data} (missing values retained).
#'
#' @noRd
get_decorrelate_response <- function(formula, data) {
  response_formula <- formula
  response_formula[[3]] <- 1
  frame <- model.frame(response_formula, data = data, na.action = na.pass)
  as.numeric(model.response(frame))
}

#' Validate and normalize the \code{training} argument for grid evaluation
#'
#' Resolves \code{training} into one or more training/test splits over the
#' fitted (non-missing-response) rows: \code{method = "split"} builds one or
#' more random splits (or uses explicitly supplied
#' \code{training_index}/\code{test_index}), and \code{method = "cv"} builds
#' k-fold splits (or uses an explicitly supplied \code{folds_index}).
#'
#' @param training A list describing the training design; see the
#'   \code{training} argument to \code{\link{ssn_decorrelate}()}, or
#'   \code{NULL} for the default single 80/20 split.
#' @param response_index An integer index into the original data selecting
#'   fitted (non-missing-response) rows, from
#'   \code{\link{get_decorrelate_response_index}()}.
#' @param original_n The number of rows in the original (unsubset) data,
#'   used to validate user-supplied row numbers.
#'
#' @return A list with \code{method}, \code{splits} (a named list of
#'   \code{training_index}/\code{test_index} pairs, indices into
#'   \code{response_index}), \code{response_index}, and \code{folds_index}
#'   (non-\code{NULL} only for \code{method = "cv"}).
#'
#' @noRd
get_decorrelate_training <- function(training, response_index, original_n) {
  if (is.null(training)) training <- list()
  if (!is.list(training)) stop("training must be a list.", call. = FALSE)
  method <- training[["method", exact = TRUE]]
  if (is.null(method)) method <- "split"
  n <- length(response_index)
  index <- seq_len(n)
  map_original <- function(values, label) {
    if (!is.numeric(values) || anyNA(values) || any(values != as.integer(values)) ||
        any(values < 1L | values > original_n) || anyDuplicated(values)) {
      stop("training$", label, " must contain unique valid original observation row numbers.", call. = FALSE)
    }
    mapped <- match(as.integer(values), response_index)
    if (anyNA(mapped)) {
      stop("training$", label, " cannot include rows whose response is missing.", call. = FALSE)
    }
    mapped
  }
  if (identical(method, "split")) {
    training_index <- training[["training_index", exact = TRUE]]
    test_index <- training[["test_index", exact = TRUE]]
    if (!is.null(training_index) || !is.null(test_index)) {
      train <- if (is.null(training_index)) NULL else map_original(training_index, "training_index")
      test <- if (is.null(test_index)) NULL else map_original(test_index, "test_index")
      if (is.null(train)) train <- setdiff(index, test)
      if (is.null(test)) test <- setdiff(index, train)
      splits <- list(list(training_index = train, test_index = test))
    } else {
      p <- training[["p", exact = TRUE]]
      if (is.null(p)) p <- 0.8
      replicate <- training[["replicate", exact = TRUE]]
      if (is.null(replicate)) replicate <- 1L
      replicate <- as.integer(replicate)
      if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1 || replicate < 1L) {
        stop("training$p must be between zero and one and training$replicate must be positive.", call. = FALSE)
      }
      splits <- lapply(seq_len(replicate), function(i) {
        train <- sample(index, floor(p * n))
        list(training_index = sort(train), test_index = setdiff(index, train))
      })
    }
  } else if (identical(method, "cv")) {
    folds_index <- training[["folds_index", exact = TRUE]]
    fold_alias <- training[["fold", exact = TRUE]]
    if (!is.null(folds_index) && !is.null(fold_alias)) {
      stop("Supply only one of training$folds_index and its training$fold alias.", call. = FALSE)
    }
    fold <- if (is.null(folds_index)) fold_alias else folds_index
    if (is.null(fold)) {
      folds <- training[["folds", exact = TRUE]]
      if (is.null(folds)) folds <- 5L
      folds <- as.integer(folds)
      if (folds < 2L || folds > n) stop("training$folds must be between 2 and the sample size.", call. = FALSE)
      fold <- sample(rep(seq_len(folds), length.out = n))
    } else if (length(fold) == original_n) {
      fold <- fold[response_index]
    }
    if (length(fold) != n || anyNA(fold)) {
      stop("training$folds_index must have one non-missing value per fitted response or original row.", call. = FALSE)
    }
    if (length(unique(fold)) < 2L) {
      stop("training$folds_index must define at least two non-empty folds.", call. = FALSE)
    }
    splits <- lapply(unique(fold), function(value) {
      test <- which(fold == value)
      list(training_index = setdiff(index, test), test_index = test)
    })
  } else {
    stop("training$method must be \"split\" or \"cv\".", call. = FALSE)
  }
  for (split in splits) {
    if (!length(split$training_index) || !length(split$test_index) ||
        anyDuplicated(c(split$training_index, split$test_index)) ||
        !all(c(split$training_index, split$test_index) %in% index)) {
      stop("Every training split must contain non-empty, disjoint training and test sets with valid observed rows.", call. = FALSE)
    }
  }
  names(splits) <- if (identical(method, "cv")) as.character(unique(fold)) else as.character(seq_along(splits))
  list(method = method, splits = splits, response_index = response_index,
       folds_index = if (identical(method, "cv")) fold else NULL)
}

#' Evaluate every covariance candidate over every training/test split
#'
#' Fits and evaluates each candidate (via
#' \code{\link{get_decorrelate_candidate_statistics}()}) on each training
#' split, so results can later be averaged per candidate across replications
#' or folds by \code{\link{get_decorrelate_grid_summary}()}.
#'
#' @param candidates A named list of covariance candidates.
#' @param training A resolved training-design list from
#'   \code{\link{get_decorrelate_training}()}.
#' @param formula,ssn.object,additive,anisotropy,random,partition_factor,ordering,algorithm,local
#'   See \code{\link{ssn_decorrelate}()}.
#' @param tailup_type,taildown_type,euclid_type,nugget_type,statistic Unused;
#'   accepted only for a call signature matching the caller's own arguments.
#' @param dots Additional algorithm-specific arguments (as a list).
#'
#' @return A data frame with one row per candidate/split, containing
#'   \code{candidate}, \code{split}, and the evaluation statistics from
#'   \code{\link{get_decorrelate_candidate_statistics}()}.
#'
#' @noRd
get_decorrelate_grid_evaluation <- function(candidates, training, formula, ssn.object,
                                             tailup_type, taildown_type, euclid_type, nugget_type,
                                             additive, anisotropy, random, partition_factor, ordering,
                                             algorithm, statistic, local, dots) {
  result <- lapply(names(candidates), function(candidate_name) {
    candidate <- candidates[[candidate_name]]
    args <- get_decorrelate_data_args(
      formula, ssn.object,
      candidate$tailup_initial, candidate$taildown_initial,
      candidate$euclid_initial, candidate$nugget_initial,
      additive, candidate$randcov_initial, partition_factor
    )
    do.call(rbind, lapply(names(training$splits), function(split_name) {
      split <- training$splits[[split_name]]
      statistics <- tryCatch(
        get_decorrelate_candidate_statistics(
          args, ssn.object, formula, training$response_index, split,
          ordering, local, algorithm, dots
        ),
        error = function(error) data.frame(
          bias = NA_real_, MSPE = NA_real_, RMSPE = NA_real_, cor2 = NA_real_,
          failure = conditionMessage(error)
        )
      )
      data.frame(candidate = candidate_name, split = split_name, statistics, row.names = NULL)
    }))
  })
  do.call(rbind, result)
}

#' Fit one candidate on a training split and evaluate it on the held-out test rows
#'
#' Marks the split's test rows missing, fits \code{\link{ssn_decorrelate_data}()}
#' and a learner on the remaining (training) rows, predicts the held-out rows
#' via \code{\link{ssn_decorrelate_newdata}()}/\code{\link{ssn_recorrelate_newdata}()},
#' and computes bias/MSPE/RMSPE/cor2 against the observed held-out response.
#'
#' @param args A candidate's \code{\link{ssn_decorrelate_data}()} argument
#'   list from \code{\link{get_decorrelate_data_args}()}.
#' @param ssn.object The full (unmodified) fitted-model-ready SSN object.
#' @param formula A model formula.
#' @param response_index An integer index into \code{ssn.object$obs}
#'   selecting fitted (non-missing-response) rows.
#' @param split A single training/test split (indices into
#'   \code{response_index}), from \code{\link{get_decorrelate_training}()}.
#' @param ordering The resolved row-ordering method.
#' @param local A resolved \code{local} list.
#' @param algorithm One of \code{"ranger"}, \code{"randomForest"}, or
#'   \code{"xgboost"}.
#' @param dots Additional algorithm-specific arguments (as a list).
#'
#' @return A one-row data frame with \code{bias}, \code{MSPE}, \code{RMSPE},
#'   \code{cor2}, and \code{failure} (\code{NA} unless an error occurred).
#'
#' @noRd
get_decorrelate_candidate_statistics <- function(args, ssn.object, formula,
                                                  response_index, split, ordering,
                                                  local, algorithm, dots) {
  training_original <- response_index[as.integer(split$training_index)]
  held_original <- response_index[as.integer(split$test_index)]
  held_ssn <- ssn.object
  excluded_original <- setdiff(response_index, training_original)
  response_variables <- all.vars(formula[[2]])
  for (variable in response_variables) {
    held_ssn$obs[[variable]][excluded_original] <- NA
  }
  args$ssn.object <- held_ssn
  args$ordering <- ordering
  args$local <- local
  object <- do.call(ssn_decorrelate_data, args)
  learner <- fit_decorrelate_algorithm(object$tX, object$ty, algorithm, dots)
  # missing_index reproduces exactly what get_decorrelate_context() computed
  # internally (via restruct_ssn_missing()) to build .missing from held_ssn$obs,
  # so held_original's positions within it can be found without a row-key lookup
  missing_index <- which(is.na(held_ssn$obs[[all.vars(formula)[[1]]]]))
  held_match <- match(held_original, missing_index)
  if (anyNA(held_match)) {
    stop("Held-out rows are not aligned with the decorrelation prediction set.", call. = FALSE)
  }
  missing_data <- object$covariance_fit$ssn.object$preds$.missing
  object$covariance_fit$ssn.object$preds$.missing <-
    missing_data[held_match, , drop = FALSE]
  transformed_test <- ssn_decorrelate_newdata(
    object, ".missing", local = local
  )
  decorated_prediction <- predict_decorrelate_algorithm(
    learner, transformed_test$tX_newdata, algorithm
  )
  prediction <- ssn_recorrelate_newdata(transformed_test, decorated_prediction)
  observed <- get_decorrelate_response(formula, ssn.object$obs)[held_original]
  error <- observed - prediction
  data.frame(
    bias = mean(error), MSPE = mean(error^2), RMSPE = sqrt(mean(error^2)),
    cor2 = if (length(prediction) < 2L || sd(prediction) == 0 || sd(observed) == 0) NA_real_ else cor(prediction, observed)^2,
    failure = NA_character_
  )
}

#' Summarize per-split evaluation statistics into a per-candidate ranked table
#'
#' Averages each candidate's statistics across its splits (folds/
#' replications), treats a candidate as failed if any of its splits errored,
#' drops candidates without a finite \code{statistic}, and orders the result
#' (best first: decreasing for \code{"cor2"}, increasing absolute value for
#' \code{"bias"}, increasing otherwise).
#'
#' @param evaluation The per-candidate/split data frame from
#'   \code{\link{get_decorrelate_grid_evaluation}()}.
#' @param statistic The statistic used to rank candidates: one of
#'   \code{"bias"}, \code{"MSPE"}, \code{"RMSPE"}, or \code{"cor2"}.
#'
#' @return A ranked data frame with one row per (non-failed, finite)
#'   candidate, its averaged statistics, and a \code{"statistic"} attribute.
#'
#' @noRd
get_decorrelate_grid_summary <- function(evaluation, statistic) {
  statistics <- c("bias", "MSPE", "RMSPE", "cor2")
  summary <- aggregate(evaluation[statistics], list(candidate = evaluation$candidate), mean, na.rm = TRUE)
  for (name in statistics) summary[[name]][!is.finite(summary[[name]])] <- NA_real_
  failed <- unique(evaluation$candidate[
    !is.na(evaluation$failure) & nzchar(evaluation$failure)
  ])
  summary[[statistic]][summary$candidate %in% failed] <- NA_real_
  valid <- is.finite(summary[[statistic]])
  if (!any(valid)) {
    failures <- unique(evaluation$failure[!is.na(evaluation$failure) & nzchar(evaluation$failure)])
    detail <- if (length(failures)) paste0(" Evaluation failures: ", paste(failures, collapse = "; ")) else ""
    stop("No covariance candidate produced a finite ", statistic, ".", detail, call. = FALSE)
  }
  summary <- summary[valid, , drop = FALSE]
  if (identical(statistic, "cor2")) {
    summary <- summary[order(summary[[statistic]], decreasing = TRUE), , drop = FALSE]
  } else if (identical(statistic, "bias")) {
    summary <- summary[order(abs(summary[[statistic]])), , drop = FALSE]
  } else {
    summary <- summary[order(summary[[statistic]]), , drop = FALSE]
  }
  rownames(summary) <- NULL
  attr(summary, "statistic") <- statistic
  summary
}

#' Convert a candidates list into a printable grid \code{data.frame}
#'
#' The inverse of \code{\link{get_decorrelate_candidates_from_table}()}: for
#' each candidate, flattens its \code{*_initial}/\code{randcov_initial}
#' objects into one row of type and parameter columns (\code{tailup_type},
#' \code{tailup_de}, \code{tailup_range}, ..., \code{randcov_<name>}, ...),
#' filling any column absent from a given candidate with \code{NA} so every
#' row shares the same columns.
#'
#' @param candidates A named list of covariance candidates.
#'
#' @return A data frame with one row per candidate and a \code{candidate}
#'   identifier column (removed by callers before returning the grid to
#'   users).
#'
#' @noRd
get_decorrelate_grid_parameters <- function(candidates) {
  rows <- lapply(names(candidates), function(candidate_name) {
    candidate <- candidates[[candidate_name]]
    values <- list(
      candidate = candidate_name,
      tailup_type = remove_covtype(class(candidate$tailup_initial)),
      taildown_type = remove_covtype(class(candidate$taildown_initial)),
      euclid_type = remove_covtype(class(candidate$euclid_initial)),
      nugget_type = remove_covtype(class(candidate$nugget_initial))
    )
    for (component in c("tailup", "taildown", "euclid", "nugget")) {
      initial <- candidate[[paste0(component, "_initial")]]
      initial <- switch(component,
        tailup = tailup_initial_NA(initial),
        taildown = taildown_initial_NA(initial),
        euclid = euclid_initial_NA(initial, list(anisotropy = TRUE)),
        nugget = nugget_initial_NA(initial)
      )$initial
      for (name in names(initial)) values[[paste(component, name, sep = "_")]] <- initial[[name]]
    }
    if (!is.null(candidate$randcov_initial)) {
      for (name in names(candidate$randcov_initial$initial)) {
        values[[paste0("randcov_", name)]] <- candidate$randcov_initial$initial[[name]]
      }
    }
    as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  })
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (name in missing) row[[name]] <- NA
    row[columns]
  })
  do.call(rbind, rows)
}
