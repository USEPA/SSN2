#' Simulate random variables on a stream network
#'
#' @description Simulate random variables on a stream
#'   network with a specific mean and covariance structure. Designed to use
#'   \code{ssn_simulate()}, but individual simulation functions for each
#'   response distribution also exist.
#'
#' @param family The response distribution family. The default is \code{"Gaussian"}.
#' @param ssn.object A spatial stream network object with class \code{SSN}. Random
#'   variables are simulated for each row of \code{ssn.object$obs}.
#' @param network The spatial stream network to simulate on. Currently only
#'   allowed to be \code{"obs"} for the \code{ssn.object$obs} object.
#' @param tailup_params An object from [tailup_params()] specifying
#'   the tailup covariance parameters.
#' @param taildown_params An object from [taildown_params()] specifying
#'   the taildown covariance parameters.
#' @param euclid_params An object from [euclid_params()] specifying
#'   the Euclidean covariance parameters.
#' @param nugget_params An object from [nugget_params()] specifying
#'   the nugget covariance parameters.
#' @param additive The name of the variable in \code{ssn.object} that is used
#'   to define spatial weights. Can be quoted or unquoted. For the tailup covariance functions, these additive
#'   weights are used for branching. Technical details that describe the role
#'   of the additive variable in the tailup covariance function are available
#'   in Ver Hoef and Peterson (2010).
#' @param mean A numeric vector representing the mean. \code{mean} must have length 1
#'   (in which case it is recycled) or length equal
#'   to the number of rows in \code{ssn.object$obs}. The default is \code{0}.
#'   For non-Gaussian families, this is the mean of the latent normal random variable
#'   on the link scale. See Details for more.
#' @param samples The number of independent samples to generate. The default
#'   is \code{1}.
#' @param dispersion The dispersion value for the beta, gamma, negative binomial,
#'   and inverse Gaussian distributions. The default is \code{1}.
#' @param size A numeric vector representing the sample size for each binomial trial.
#'   The default is \code{1}, which corresponds to a Bernoulli trial for each observation.
#' @param randcov_params A [spmodel::randcov_params()] object.
#' @param partition_factor A formula indicating the partition factor.
#' @param local An optional
#'   logical or list controlling the big data approximation. If omitted, \code{local}
#'   is \code{TRUE} when the desired sample size exceeds
#'   5,000, and \code{FALSE} otherwise. \code{FALSE} uses exact simulation.
#'   \code{TRUE} or a list uses a local neighborhood (Vecchia) approach,
#'   simulating each observation
#'   conditional on neighbors selected from previously simulated
#'   observations. A list can contain:
#'   \itemize{
#'     \item \code{method}: The method used to build each
#'       location's local neighborhood set once it exceeds \code{size}
#'       already-simulated candidates. Values are \code{"all"} (no truncation,
#'       exact)
#'       or \code{"covariance"} (the \code{size} candidates with the
#'       highest covariance, in absolute value, with the location being
#'       simulated). This is the same convention as \code{predict()}'s own
#'       \code{local$method}. The default is \code{"covariance"}.
#'       \code{method = "all"} is an alias for \code{local = FALSE}.
#'     \item \code{size}: The number of neighbors used when \code{method}
#'       is \code{"covariance"}. The default is 30.
#'     \item \code{ordering}: The order locations are simulated in --
#'       \code{"pid"}, \code{"maxmin"}, \code{"middleout"},
#'       \code{"outsidein"}, \code{"coordinate"}, \code{"grts"},
#'       \code{"random"}, or \code{"none"} (these are the same options as
#'       \code{ssn_decorrelate()}'s \code{ordering} argument). The default is
#'       \code{"pid"}.
#'   }
#' @param ... Additional arguments passed to [ssn_rnorm()], including \code{local},
#'   when using \code{ssn_simulate()} or another simulation function.
#'
#' @details Random variables are simulated via the product of the covariance matrix's
#'   square (Cholesky) root and independent standard normal random variables
#'   with mean 0 and variance 1, with big data approximations available.
#'
#'   If not using \code{ssn_simulate()}, individual simulation functions for
#'   each response distribution do exist:
#'
#'   \itemize{
#'     \item \code{ssn_rnorm()}: Simulate from a Gaussian distribution
#'     \item \code{ssn_rpois()}: Simulate from a Poisson distribution
#'     \item \code{ssn_rnbinom()}: Simulate from a negative binomial distribution
#'     \item \code{ssn_rbinom()}: Simulate from a binomial distribution
#'     \item \code{ssn_rbeta()}: Simulate from a beta distribution
#'     \item \code{ssn_rgamma()}: Simulate from a gamma distribution
#'     \item \code{ssn_rinvgauss()}: Simulate from an inverse Gaussian distribution
#'
#'   }
#'
#' @return If \code{samples} is 1, a vector of random variables for each row of \code{ssn.object$obs}
#'   is returned. If \code{samples} is greater than one, a matrix of random variables
#'   is returned, where the rows correspond to each row of \code{ssn.object$obs} and the columns
#'   correspond to independent samples.
#'
#' @export
#'
#' @order 1
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' tailup <- tailup_params("exponential", de = 0.1, range = 200)
#' taildown <- taildown_params("exponential", de = 0.4, range = 300)
#' euclid <- euclid_params("spherical", de = 0.2, range = 1000, rotate = 0, scale = 1)
#' nugget <- nugget_params("nugget", nugget = 0.1)
#' ssn_simulate("gaussian", mf04p, "obs", tailup, taildown, euclid, nugget, additive = "afvArea")
#' ssn_simulate("poisson", mf04p, "obs", tailup, taildown, euclid, nugget, additive = "afvArea")
#'
#' @references
#' Ver Hoef, J.M. and Peterson, E.E. (2010) A moving average approach for spatial
#' statistical models of stream networks (with discussion).
#' \emph{Journal of the American Statistical Association} \bold{105}, 6--18.
#' DOI: 10.1198/jasa.2009.ap08248.  Rejoinder pgs. 22--24.
ssn_simulate <- function(family = "Gaussian", ssn.object, network = "obs", tailup_params,
                         taildown_params, euclid_params, nugget_params,
                         additive, mean = 0, samples = 1, dispersion = 1, size = 1,
                         randcov_params, partition_factor, ...) {
  # fix family
  if (is.symbol(substitute(family))) { # or is.language
    family <- deparse1(substitute(family))
  }

  fn_calls <- c(
    "Gaussian" = "ssn_rnorm",
    "gaussian" = "ssn_rnorm",
    "binomial" = "ssn_rbinom",
    "beta" = "ssn_rbeta",
    "poisson" = "ssn_rpois",
    "nbinomial" = "ssn_rnbinom",
    "Gamma" = "ssn_rgamma",
    "inverse.gaussian" = "ssn_rinvgauss"
  )

  call_val <- match.call()
  call_val[[1]] <- as.symbol(fn_calls[family])
  call_list <- as.list(call_val)
  call_list$family <- NULL
  call_val <- as.call(call_list)
  ssn_simulate_val <- eval(call_val, envir = parent.frame())
  ssn_simulate_val
}
