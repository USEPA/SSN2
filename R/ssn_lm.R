#' Fitting Linear Models for Spatial Stream Networks
#'
#' @description This function works on spatial stream network objects to fit
#'   linear models with spatially autocorrelated errors using likelihood methods, allowing for
#'   non-spatial random effects, anisotropy, partition factors, big data methods, and more.
#'   The spatial formulation is described in Ver Hoef and Peterson (2010)
#'   and Peterson and Ver Hoef (2010).
#'
#' @param formula A two-sided linear formula describing the fixed effect structure
#'   of the model, with the response to the left of the \code{~} operator and
#'   the terms on the right, separated by \code{+} operators.
#' @param ssn.object A spatial stream network object with class \code{SSN}.
#' @param tailup_type The tailup covariance function type. Available options
#'   include \code{"linear"}, \code{"spherical"}, \code{"exponential"},
#'   \code{"mariah"}, \code{"epa"}, \code{"gaussian"}, and \code{"none"}. Parameterizations are
#'   described in Details.
#' @param taildown_type The taildown covariance function type. Available options
#'   include \code{"linear"}, \code{"spherical"}, \code{"exponential"},
#'   \code{"mariah"}, \code{"epa"}, \code{"gaussian"}, and \code{"none"}. Parameterizations are
#'   described in Details.
#' @param euclid_type The euclidean covariance function type. Available options
#'   include \code{"spherical"}, \code{"exponential"}, \code{"gaussian"},
#'   \code{"cosine"}, \code{"cubic"}, \code{"pentaspherical"}, \code{"wave"},
#'    \code{"jbessel"}, \code{"gravity"}, \code{"rquad"}, \code{"magnetic"}, and
#'    \code{"none"}. Parameterizations are
#'   described in Details.
#' @param nugget_type The nugget covariance function type. Available options
#'   include \code{"nugget"} or \code{"none"}. Parameterizations are
#'   described in Details.
#' @param tailup_initial An object from [tailup_initial()] specifying initial and/or
#'   known values for the tailup covariance parameters.
#' @param taildown_initial An object from [taildown_initial()] specifying initial and/or
#'   known values for the taildown covariance parameters.
#' @param euclid_initial An object from [euclid_initial()] specifying initial and/or
#'   known values for the euclidean covariance parameters.
#' @param nugget_initial An object from [nugget_initial()] specifying initial and/or
#'   known values for the nugget covariance parameters.
#' @param additive The name of the variable in \code{ssn.object} that is used
#'   to define spatial weights. Can be quoted or unquoted. For the tailup covariance functions, these additive
#'   weights are used for branching. Technical details that describe the role
#'   of the additive variable in the tailup covariance function are available
#'   in Ver Hoef and Peterson (2010).
#' @param estmethod The estimation method. Available options include
#'   \code{"reml"} for restricted maximum likelihood and \code{"ml"} for maximum
#'   likelihood. The default is \code{"reml"}.
#' @param anisotropy A logical indicating whether (geometric) anisotropy should
#'   be modeled. Not required if \code{spcov_initial} is provided with 1) \code{rotate}
#'   assumed unknown or assumed known and non-zero or 2) \code{scale} assumed unknown
#'   or assumed known and less than one. When \code{anisotropy} is \code{TRUE},
#'   computational times can significantly increase. The default is \code{FALSE}.
#' @param random A one-sided linear formula describing the random effect structure
#'   of the model. Terms are specified to the right of the \code{~ operator}.
#'   Each term has the structure \code{x1 + ... + xn | g1/.../gm}, where \code{x1 + ... + xn}
#'   specifies the model for the random effects and \code{g1/.../gm} is the grouping
#'   structure. Separate terms are separated by \code{+} and must generally
#'   be wrapped in parentheses. Random intercepts are added to each model
#'   implicitly when at least  one other variable is defined.
#'   If a random intercept is not desired, this must be explicitly
#'   defined (e.g., \code{x1 + ... + xn - 1 | g1/.../gm}). If only a random intercept
#'   is desired for a grouping structure, the random intercept must be specified
#'   as \code{1 | g1/.../gm}. Note that \code{g1/.../gm} is shorthand for \code{(1 | g1/.../gm)}.
#'   If only random intercepts are desired and the shorthand notation is used,
#'   parentheses can be omitted.
#' @param randcov_initial An optional object specifying initial and/or
#'   known values for the random effect variances. See [spmodel::randcov_initial()].
#' @param partition_factor A one-sided linear formula with a single term
#'   specifying the partition factor.  The partition factor assumes observations
#'   from different levels of the partition factor are uncorrelated.
#' @param local An optional logical or list controlling the big data approximation.
#'   \code{local} can only be used when big data distance matrices have been created
#'   using [ssn_create_bigdist()] and is most beneficial when the sample size is
#'   at least 5,000 \code{ssn_lm()} or 3,000 \code{ssn_glm()}.
#'   If a list is provided, the following arguments detail the big
#'   data approximation:
#'   \itemize{
#'     \item \code{index: }The group indexes. Observations in different
#'       levels of \code{index} are assumed to be uncorrelated for the
#'       purposes of estimation. If \code{index} is not provided, it is
#'       determined by specifying \code{method} and either \code{size} or \code{groups}.
#'     \item \code{method}: The big data approximation method used to determine \code{index}. Ignored
#'       if \code{index} is provided. If \code{method = "random"},
#'       observations are randomly assigned to \code{index} based on \code{size}.
#'       If \code{method = "kmeans"}, observations assigned to \code{index}
#'       based on k-means clustering on the coordinates with \code{groups} clusters. The default
#'       is \code{"kmeans"}. Note that both methods have a random component, which
#'       means that you may get different results from separate model fitting calls.
#'       To ensure consistent results, specify \code{index} or set a seed via
#'       \code{base::set.seed()}.
#'     \item \code{size}: The number of observations in each \code{index} group
#'       when \code{method} is \code{"random"}. If the number of observations
#'       is not divisible by \code{size}, some levels get \code{size - 1} observations.
#'       The default is 200.
#'     \item \code{groups: }The number of \code{index} groups. If \code{method}
#'       is \code{"random"}, \code{size} is \eqn{ceiling(n / groups)}, where
#'       \eqn{n} is the sample size. Automatically determined if \code{size}
#'       is specified. If \code{method} is \code{"kmeans"}, \code{groups}
#'       is the number of clusters.
#'     \item \code{var_adjust: }The approach for adjusting the variance-covariance
#'       matrix of the fixed effects. \code{"none"} for no adjustment, \code{"theoretical"}
#'       for the theoretically-correct adjustment,
#'       \code{"pooled"} for the pooled adjustment, and \code{"empirical"} for the
#'       empirical adjustment. The default is \code{"theoretical"} for samples sizes
#'       up to 100,000 and \code{"none"} for samples sizes exceeding 100,000.
#'     \item \code{parallel}: If \code{TRUE}, parallel processing via the
#'       parallel package is automatically used. The default is \code{FALSE}.
#'     \item \code{ncores}: If \code{parallel = TRUE}, the number of cores to
#'       parallelize over. The default is the number of available cores on your machine.
#'   }
#'   When \code{local} is a list, at least one list element must be provided to
#'   initialize default arguments for the other list elements.
#'   If \code{local} is \code{TRUE}, defaults for \code{local} are chosen such
#'   that \code{local} is transformed into
#'   \code{list(size = 200, method = "kmeans", var_adjust = "theoretical", parallel = FALSE)}.
#' @param ... Other arguments to \code{stats::optim()}.
#'
#' @details The linear model for spatial stream networks can be written as
#'   \eqn{y = X \beta + zu + zd + ze + n}, where \eqn{X} is the fixed effects design
#'   matrix, \eqn{\beta} are the fixed effects, \eqn{zu} is tailup random error,
#'   \eqn{zd} is taildown random error, and \eqn{ze} is Euclidean random error,
#'   and \eqn{n} is nugget random error. The tailup random errors capture spatial
#'   covariance moving downstream (and depend on downstream distance), the taildown
#'   random errors capture spatial covariance moving upstream (and depend on upstream)
#'   distance, the Euclidean random errors capture spatial covariance that depends on
#'   Euclidean distance, and the nugget random errors captures variability
#'   independent of spatial locations. The response \eqn{y} is modeled using a
#'   spatial covariance function expressed as
#'   \eqn{de(zu) * R(zu) + de(zd) * R(zd) + de(ze) * R(ze) + nugget * I}.
#'   \eqn{de(zu)}, \eqn{de(zu)}, and \eqn{de(zd)} represent the tailup, taildown, and Euclidean
#'   variances, respectively. \eqn{R(zu)}, \eqn{R(zd)}, and \eqn{R(ze)} represent the tailup,
#'   taildown, and Euclidean correlation matrices, respectively. Each correlation
#'   matrix depends on a range parameter that controls the distance-decay behavior
#'   of the correlation. \eqn{nugget} represents the nugget variance and
#'   \eqn{I} represents an identity matrix.
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
#'     \item cosine: \eqn{cos(r)}
#'     \item wave: \eqn{sin(r) * (h > 0) / r + (h = 0)}
#'     \item jbessel: \eqn{Bj(h * range)}, Bj is Bessel-J function
#'     \item gravity: \eqn{(1 + r^2)^{-0.5}}
#'     \item rquad: \eqn{(1 + r^2)^{-1}}
#'     \item magnetic: \eqn{(1 + r^2)^{-1.5}}
#'     \item none: \eqn{I}
#'   }
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
#'   and \code{scale} are assumed unknown, requiring estimation. If either \code{rotate} or \code{scale}
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
#'   can be written as \eqn{y = X \beta + W1\gamma 1 + ... Wj\gamma j + zu + zd + ze + n},
#'   where each Z is a random effects design matrix and each u is a random effect.
#'
#'  \code{partition_factor} Details: The partition factor can be represented in matrix form as \eqn{P}, where
#'   elements of \eqn{P} equal one for observations in the same level of the partition
#'   factor and zero otherwise. The covariance matrix involving only the
#'   spatial and random effects components is then multiplied element-wise
#'   (Hadmard product) by \eqn{P}, yielding the final covariance matrix.
#'
#' \code{local} Details: The big data approximation works by sorting observations into different levels
#'   of an index variable. Observations in different levels of the index variable
#'   are assumed to be uncorrelated for the purposes of model fitting. Sparse matrix methods are then implemented
#'   for significant computational gains. Parallelization generally further speeds up
#'   computations when data sizes are larger than a few thousand. Both the \code{"random"} and \code{"kmeans"} values of \code{method}
#'   in \code{local} have random components. That means you may get slightly different
#'   results when using the big data approximation and rerunning \code{ssn_lm()} with the same code. For consistent results,
#'   either set a seed via \code{base::set.seed()} or specify \code{index} to \code{local}.
#'
#'   Other Details: Observations with \code{NA} response values are removed for model
#'   fitting, but their values can be predicted afterwards by running
#'   \code{predict(object)}.
#'
#'
#' @return A list with many elements that store information about
#'   the fitted model object and has class \code{ssn_lm}. Many generic functions that
#'   summarize model fit are available for \code{ssn_lm} objects, including
#'   \code{AIC}, \code{AICc}, \code{anova}, \code{augment}, \code{coef},
#'   \code{cooks.distance}, \code{covmatrix}, \code{deviance}, \code{fitted}, \code{formula},
#'   \code{glance}, \code{glances}, \code{hatvalues}, \code{influence},
#'   \code{labels}, \code{logLik}, \code{loocv}, \code{model.frame}, \code{model.matrix},
#'   \code{plot}, \code{predict}, \code{print}, \code{pseudoR2}, \code{summary},
#'   \code{terms}, \code{tidy}, \code{update}, \code{varcomp}, and \code{vcov}.
#'
#'   This fitted model list contains the following elements:
#'   \itemize{
#'     \item \code{additive}: The name of the additive function value column.
#'     \item \code{anisotropy}: Whether euclidean anisotropy was modeled.
#'     \item \code{call}: The function call.
#'     \item \code{coefficients}: Model coefficients.
#'     \item \code{contrasts}: Any user-supplied contrasts.
#'     \item \code{cooks_distance}: Cook's distance values.
#'     \item \code{crs}: The geographic coordinate reference system.
#'     \item \code{deviance}: The model deviance.
#'     \item \code{diagtol}: A tolerance value that may be added to the diagonal
#'       of covariance matrices to encourage decomposition stability.
#'     \item \code{estmethod}: The estimation method.
#'     \item \code{euclid_max}: The maximum euclidean distance.
#'     \item \code{fitted}: Fitted values.
#'     \item \code{formula}: The model formula.
#'     \item \code{hatvalues}: The hat (leverage) values.
#'     \item \code{is_known}: An object that identifies which parameters are known.
#'     \item \code{local_index}: An index identifier used internally for sorting.
#'     \item \code{missing_index}: Which rows in the "obs" object had missing responses.
#'     \item \code{n}: The sample size.
#'     \item \code{npar}: The number of estimated covariance parameters.
#'     \item \code{observed_index}: Which rows in the "obs" object had observed responses.
#'     \item \code{optim}: The optimization output.
#'     \item \code{p}: The number of fixed effects.
#'     \item \code{partition_factor}: The partition factor formula.
#'     \item \code{pseudoR2}: The pseudo R-squared.
#'     \item \code{random}: The random effect formula.
#'     \item \code{residuals}: The residuals.
#'     \item \code{sf_column_name}: The name of the geometry columns \code{ssn.object}
#'     \item \code{ssn.object}: An updated \code{ssn.object}.
#'     \item \code{tail_max}: The maximum stream distance.
#'     \item \code{terms}: The model terms.
#'     \item \code{vcov}: Variance-covariance matrices
#'     \item \code{xlevels}: The levels of factors in the model matrix.
#'   }
#'
#'   These list elements are meant to be used with various generic functions
#'   (\code{e.g., residuals()} that operate on the model object.
#'   While possible to access elements of the fitted model list directly, we strongly
#'   advise against doing so when there is a generic available to return the element
#'   of interest. For example, we strongly recommend using \code{residuals()} to
#'   obtain model residuals instead of accessing the fitted model list directly via
#'   \code{object$residuals}.
#'
#' @note This function does not perform any internal scaling. If optimization is not
#'   stable due to large extremely large variances, scale relevant variables
#'   so they have variance 1 before optimization.
#'
#' @export
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' summary(ssn_mod)
#'
#' @references
#' Garreta, V., Monestiez, P. and Ver Hoef, J.M. (2010) Spatial modelling and
#' prediction on river networks: up model, down model, or hybrid?
#' \emph{Environmetrics} \bold{21(5)}, 439--456.
#'
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach
#' to geostatistical modeling in stream networks. \emph{Ecology} \bold{91(3)},
#' 644--651.
#'
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial
#' statistical models of stream networks (with discussion).
#' \emph{Journal of the American Statistical Association} \bold{105}, 6--18.
#' DOI: 10.1198/jasa.2009.ap08248.  Rejoinder pgs. 22--24.
ssn_lm <- function(formula, ssn.object,
                   tailup_type = "none", taildown_type = "none",
                   euclid_type = "none", nugget_type = "nugget",
                   tailup_initial, taildown_initial, euclid_initial, nugget_initial,
                   additive, estmethod = "reml", anisotropy = FALSE,
                   random, randcov_initial, partition_factor, local, ...) {
  # set defaults
  if (missing(tailup_initial)) tailup_initial <- NULL
  if (missing(taildown_initial)) taildown_initial <- NULL
  if (missing(euclid_initial)) euclid_initial <- NULL
  if (missing(nugget_initial)) nugget_initial <- NULL
  if (missing(additive)) additive <- NULL
  if (missing(random)) random <- NULL
  if (missing(randcov_initial)) randcov_initial <- NULL
  if (missing(partition_factor)) partition_factor <- NULL
  if (missing(local)) local <- NULL

  # fix additive
  if (is.symbol(substitute(additive))) { # or is.language
    additive <- deparse1(substitute(additive))
  }

  # set initials
  initial_object <- get_initial_object(
    tailup_type, taildown_type, euclid_type, nugget_type,
    tailup_initial, taildown_initial, euclid_initial, nugget_initial
  )

  # perform checks to return errors
  check_ssn_lm(initial_object, ssn.object, additive, estmethod)

  # get data object
  if (is.null(local) || (is.logical(local) && !local)) {
    data_object <- get_data_object(
      formula, ssn.object, additive, anisotropy,
      initial_object, random, randcov_initial, partition_factor, local, ...
    )
    if (data_object$n > 5000) message("Because the sample size exceeds 5000, consider setting local = TRUE to perform computationally efficient approximations. Ensure big data distance matrices have been created using ssn_create_bigdist().")
  } else {
    data_object <- get_data_object_bigdata(
      formula, ssn.object, additive, anisotropy,
      initial_object, random, randcov_initial, partition_factor, local, ...
    )
  }

  # add random to initial object
  initial_object$randcov_initial <- data_object$randcov_initial

  # get optim dotlist
  optim_dotlist <- get_optim_dotlist(...)

  if (data_object$parallel) {
    data_object$cl <- parallel::makeCluster(data_object$ncores)
    # invisible(clusterEvalQ(data_object$cl, library(Matrix)))
  }

  # covariance parameter estimation
  cov_est_object <- cov_estimate_gloglik(data_object, ssn.object, initial_object, estmethod, optim_dotlist)

  # compute model statistics
  if (is.null(local) || (is.logical(local) && !local)) {
    model_stats <- get_model_stats(cov_est_object, data_object, estmethod)
  } else {
    model_stats <- get_model_stats_bigdata(cov_est_object, data_object, estmethod)
  }

  # parallel cluster if necessary (add back when local implemented)
  if (data_object$parallel) {
    data_object$cl <- parallel::stopCluster(data_object$cl) # makes it NULL
  }

  # store index if necessary (add back when local implemented)
  if (is.null(local) || (is.logical(local) && !local)) { # local was stored as NULL in previous function call
    local_index <- NULL
  } else {
    local_index <- data_object$local_index
    local_index <- local_index[order(data_object$order_bigdata)]
    data_object$observed_index <- which(data_object$observed_index)
    data_object$missing_index <- which(data_object$missing_index)
  }



  output <- list(
    coefficients = model_stats$coefficients,
    fitted = model_stats$fitted,
    hatvalues = model_stats$hatvalues,
    residuals = model_stats$residuals,
    cooks_distance = model_stats$cooks_distance,
    vcov = model_stats$vcov,
    deviance = model_stats$deviance,
    pseudoR2 = model_stats$pseudoR2,
    p = data_object$p,
    n = data_object$n,
    npar = model_stats$npar,
    formula = formula,
    terms = data_object$terms,
    call = match.call(),
    estmethod = estmethod,
    anisotropy = data_object$anisotropy,
    optim = cov_est_object$optim_output,
    random = random,
    is_known = cov_est_object$is_known,
    partition_factor = partition_factor,
    observed_index = data_object$observed_index,
    missing_index = data_object$missing_index,
    contrasts = data_object$contrasts,
    xlevels = data_object$xlevels,
    local_index = local_index,
    sf_column_name = data_object$sf_column_name,
    crs = data_object$crs,
    ssn.object = data_object$ssn.object,
    additive = additive,
    diagtol = data_object$diagtol,
    tail_max = data_object$tail_max,
    euclid_max = data_object$euclid_max
  )

  output <- output[sort(names(output))]
  new_output <- structure(output, class = "ssn_lm")
  new_output
}
