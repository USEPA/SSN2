#' Fitting Generalized Linear Models for Spatial Stream Networks
#'
#' @description This function works on spatial stream network objects to fit
#'   generalized linear models with spatially autocorrelated errors using likelihood methods, allowing for
#'   non-spatial random effects, anisotropy, partition factors, big data methods, and more.
#'   The spatial formulation is described in Ver Hoef and Peterson (2010)
#'   and Peterson and Ver Hoef (2010).
#'
#' @inheritParams ssn_lm
#' @param family The generalized linear model family for use with \code{ssn_glm()}.
#'   Available options include \code{"Gaussian"}, \code{"poisson"},
#'   \code{"nbinomial"} (negative binomial), \code{"binomial"}, \code{"beta"},
#'   \code{"Gamma"}, and \code{"invgauss"}. When \code{family}
#'   is \code{"Gaussian"}, arguments are passed to and evaluated by [ssn_lm()].
#'   Can be quoted or unquoted. Note that the \code{family} argument
#'   only takes a single value, rather than the list structure used by [stats::glm].
#'   See Details for more.
#' @param dispersion_initial An object from [dispersion_initial()] specifying initial and/or
#'   known values for the tailup covariance parameters.
#'
#' @details The generalized linear model for spatial stream networks can be written as
#'   \eqn{g(\mu) = \eta = X \beta + zu + zd + ze + n}, where \eqn{\mu} is the expectation
#'   of the response given the random errors, \eqn{y}, \eqn{g()} is a function that links the mean
#'   and \eqn{\eta} (and is called a link function), \code{X} is the fixed effects design
#'   matrix, \eqn{\beta} are the fixed effects, \eqn{zu} is tailup random error,
#'   \eqn{zd} is taildown random error, and \eqn{ze} is Euclidean random error,
#'   and \eqn{n} is nugget random error.
#'
#'   There are six generalized linear model
#'   families available: \code{poisson} assumes \eqn{y} is a Poisson random variable
#'   \code{nbinomial} assumes \eqn{y} is a negative binomial random
#'   variable, \code{binomial} assumes \eqn{y} is a binomial random variable,
#'   \code{beta} assumes \eqn{y} is a beta random variable,
#'   \code{Gamma} assumes \eqn{y} is a gamma random
#'   variable, and \code{inverse.gaussian} assumes \eqn{y} is an inverse Gaussian
#'   random variable.
#'
#'   The supports for \eqn{y} for each family are given below:
#'   \itemize{
#'     \item family: support of \eqn{y}
#'     \item Gaussian: \eqn{-\infty < y < \infty}
#'     \item poisson: \eqn{0 \le y}; \eqn{y} an integer
#'     \item nbinomial: \eqn{0 \le y}; \eqn{y} an integer
#'     \item binomial: \eqn{0 \le y}; \eqn{y} an integer
#'     \item beta: \eqn{0 < y < 1}
#'     \item Gamma: \eqn{0 < y}
#'     \item inverse.gaussian: \eqn{0 < y}
#'   }
#'
#'   The generalized linear model families
#'   and the parameterizations of their link functions are given
#'   below:
#'   \itemize{
#'     \item family: link function
#'     \item Gaussian: \eqn{g(\mu) = \eta} (identity link)
#'     \item poisson: \eqn{g(\mu) = log(\eta)} (log link)
#'     \item nbinomial: \eqn{g(\mu) = log(\eta)} (log link)
#'     \item binomial: \eqn{g(\mu) = log(\eta / (1 - \eta))} (logit link)
#'     \item beta: \eqn{g(\mu) = log(\eta / (1 - \eta))} (logit link)
#'     \item Gamma: \eqn{g(\mu) = log(\eta)} (log link)
#'     \item inverse.gaussian: \eqn{g(\mu) = log(\eta)} (log link)
#'   }
#'
#'   The variance function of an individual \eqn{y} (given \eqn{\mu})
#'   for each generalized linear model family is given below:
#'   \itemize{
#'     \item family: \eqn{Var(y)}
#'     \item Gaussian: \eqn{\sigma^2}
#'     \item poisson: \eqn{\mu \phi}
#'     \item nbinomial: \eqn{\mu + \mu^2 / \phi}
#'     \item binomial: \eqn{n \mu (1 - \mu) \phi}
#'     \item beta: \eqn{\mu (1 - \mu) / (1 + \phi)}
#'     \item Gamma: \eqn{\mu^2 / \phi}
#'     \item inverse.gaussian: \eqn{\mu^2 / \phi}
#'   }
#'   The parameter \eqn{\phi} is a dispersion parameter that influences \eqn{Var(y)}.
#'   For the \code{poisson} and \code{binomial} families, \eqn{\phi} is always
#'   one. Note that this inverse Gaussian parameterization is different than a
#'   standard inverse Gaussian parameterization, which has variance \eqn{\mu^3 / \lambda}.
#'   Setting \eqn{\phi = \lambda / \mu} yields our parameterization, which is
#'   preferred for computational stability. Also note that the dispersion parameter
#'   is often defined in the literature as \eqn{V(\mu) \phi}, where \eqn{V(\mu)} is the variance
#'   function of the mean. We do not use this parameterization, which is important
#'   to recognize while interpreting dispersion estimates.
#'   For more on generalized linear model constructions, see McCullagh and
#'   Nelder (1989).
#'
#'   In the generalized linear model context, the tailup, taildown, Euclidean, and
#'   nugget covariance affect the modeled mean of an observation (conditional on
#'   these effects). On the link scale, the tailup random errors capture spatial
#'   covariance moving downstream (and depend on downstream distance), the taildown
#'   random errors capture spatial covariance moving upstream (and depend on upstream)
#'   distance, the Euclidean random errors capture spatial covariance that depends on
#'   Euclidean distance, and the nugget random errors captures variability
#'   independent of spatial locations. \eqn{\eta} is modeled using a
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
#'     \item exponential: \eqn{exp(-(r1 + r2))}
#'     \item mariah: \eqn{(log(90r1 + 1) - log(90r2 + 1)) / (90r1 - 90r2) * (A =/ B) + (1 / (90r1 + 1)) * (A = B)}
#'     \item epa: \eqn{(B - range)^2 * F2 * (r2 <= 1) / 16range^5}
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
#'  \code{random} Details: If random effects are used (the estimation method must be \code{"reml"} or
#'   \code{"ml"}), the model
#'   can be written as \eqn{g(\mu) = \eta = X \beta + W1\gamma 1 + ... Wj\gamma j + zu + zd + ze + n},
#'   where each Z is a random effects design matrix and each u is a random effect.
#'
#'  \code{partition_factor} Details: The partition factor can be represented in matrix form as \eqn{P}, where
#'   elements of \eqn{P} equal one for observations in the same level of the partition
#'   factor and zero otherwise. The covariance matrix involving only the
#'   spatial and random effects components is then multiplied element-wise
#'   (Hadmard product) by \eqn{P}, yielding the final covariance matrix.
#'
#'   Other Details: Observations with \code{NA} response values are removed for model
#'   fitting, but their values can be predicted afterwards by running
#'   \code{predict(object)}.
#'
#' @return A list with many elements that store information about
#'   the fitted model object and has class \code{ssn_glm}. Many generic functions that
#'   summarize model fit are available for \code{ssn_glm} objects, including
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
#'       of  ovariance matrices to encourage decomposition stability.
#'     \item \code{estmethod}: The estimation method.
#'     \item \code{euclid_max}: The maximum euclidean distance.
#'     \item \code{family}: The generalized linear model family
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
#'     \item \code{size}: The size of the binomial trials if relevant.
#'     \item \code{ssn.object}: An updated \code{ssn.object}.
#'     \item \code{tail_max}: The maximum stream distance.
#'     \item \code{terms}: The model terms.
#'     \item \code{vcov}: Variance-covariance matrices
#'     \item \code{xlevels}: The levels of factors in the model matrix.
#'     \item \code{y}: The response.
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
#' ssn_gmod <- ssn_glm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   family = "Gamma",
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' summary(ssn_gmod)
#'
#' @references
#' Garreta, V., Monestiez, P. and Ver Hoef, J.M. (2010) Spatial modelling and
#' prediction on river networks: up model, down model, or hybrid?
#' \emph{Environmetrics} \bold{21(5)}, 439--456.
#'
#' McCullagh P. and Nelder, J. A. (1989) \emph{Generalized Linear Models}. London: Chapman and Hall.
#'
#' Peterson, E.E. and Ver Hoef, J.M. (2010) A mixed-model moving-average approach
#' to geostatistical modeling in stream networks. \emph{Ecology} \bold{91(3)},
#' 644--651.
#'
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial
#' statistical models of stream networks (with discussion).
#' \emph{Journal of the American Statistical Association} \bold{105}, 6--18.
#' DOI: 10.1198/jasa.2009.ap08248.  Rejoinder pgs. 22--24.
ssn_glm <- function(formula, ssn.object, family,
                    tailup_type = "none", taildown_type = "none",
                    euclid_type = "none", nugget_type = "nugget",
                    tailup_initial, taildown_initial, euclid_initial, nugget_initial,
                    dispersion_initial, additive, estmethod = "reml", anisotropy = FALSE,
                    random, randcov_initial, partition_factor, ...) {
  # set defaults
  if (missing(tailup_initial)) tailup_initial <- NULL
  if (missing(taildown_initial)) taildown_initial <- NULL
  if (missing(euclid_initial)) euclid_initial <- NULL
  if (missing(nugget_initial)) nugget_initial <- NULL
  if (missing(additive)) additive <- NULL
  if (missing(dispersion_initial)) dispersion_initial <- NULL else family <- class(dispersion_initial)
  if (missing(random)) random <- NULL
  if (missing(randcov_initial)) randcov_initial <- NULL
  if (missing(partition_factor)) partition_factor <- NULL
  local <- NULL

  # fix family
  if (missing(family)) {
    stop("The family argument must be specified.", call. = FALSE)
  }
  if (is.symbol(substitute(family))) { # or is.language
    family <- deparse1(substitute(family))
  }

  # Call ssn_lm if necessary
  if (family == "Gaussian" || family == "gaussian") {
    call_val <- match.call()
    call_val[[1]] <- as.symbol("ssn_lm")
    call_list <- as.list(call_val)
    call_list <- call_list[-which(names(call_list) %in% c("family", "dispersion_initial"))]
    call_val <- as.call(call_list)
    object <- eval(call_val, envir = parent.frame())
    return(object)
  }

  # fix additive
  if (is.symbol(substitute(additive))) { # or is.language
    additive <- deparse1(substitute(additive))
  }

  # set initials
  initial_object <- get_initial_object_glm(
    tailup_type, taildown_type, euclid_type, nugget_type,
    tailup_initial, taildown_initial, euclid_initial, nugget_initial,
    family, dispersion_initial
  )

  # perform checks to return errors
  check_ssn_glm(initial_object, ssn.object, additive, estmethod)

  # get data object
  data_object <- get_data_object_glm(
    formula, ssn.object, family, additive, anisotropy,
    initial_object, random, randcov_initial, partition_factor, ...
  )

  # add random to initial object
  initial_object$randcov_initial <- data_object$randcov_initial

  # get optim dotlist
  optim_dotlist <- get_optim_dotlist(...)

  # # parallel cluster if necessary
  # if (data_object$parallel) {
  #   data_object$cl <- parallel::makeCluster(data_object$ncores)
  #   # invisible(clusterEvalQ(data_object$cl, library(Matrix)))
  # }


  # covariance parameter estimation
  cov_est_object <- cov_estimate_laploglik(data_object, ssn.object, initial_object, estmethod, optim_dotlist)

  model_stats_glm <- get_model_stats_glm(cov_est_object, data_object, estmethod)



  output <- list(
    coefficients = model_stats_glm$coefficients,
    fitted = model_stats_glm$fitted,
    hatvalues = model_stats_glm$hatvalues,
    residuals = model_stats_glm$residuals,
    cooks_distance = model_stats_glm$cooks_distance,
    vcov = model_stats_glm$vcov,
    deviance = model_stats_glm$deviance,
    pseudoR2 = model_stats_glm$pseudoR2,
    p = data_object$p,
    n = data_object$n,
    npar = model_stats_glm$npar,
    formula = formula,
    terms = data_object$terms,
    call = match.call(),
    estmethod = estmethod,
    # obdata = data_object$obdata,
    anisotropy = data_object$anisotropy,
    optim = cov_est_object$optim_output,
    random = random,
    is_known = cov_est_object$is_known,
    partition_factor = partition_factor,
    observed_index = data_object$observed_index,
    missing_index = data_object$missing_index,
    contrasts = data_object$contrasts,
    xlevels = data_object$xlevels,
    local_index = NULL,
    sf_column_name = data_object$sf_column_name,
    crs = data_object$crs,
    ssn.object = data_object$ssn.object,
    additive = additive,
    family = family,
    y = model_stats_glm$y,
    size = model_stats_glm$size,
    diagtol = data_object$diagtol,
    tail_max = data_object$tail_max,
    euclid_max = data_object$euclid_max
  )

  output <- output[sort(names(output))]
  new_output <- structure(output, class = "ssn_glm")
  new_output
}
