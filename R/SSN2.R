#' @keywords internal
"_PACKAGE"

# this is where src was used
# @useDynLib SSN2
#' @importFrom Matrix bdiag chol chol2inv colSums crossprod determinant diag Diagonal
#'   fac2sparse forceSymmetric Matrix print qr rankMatrix solve sparseMatrix t tcrossprod which
#' @importFrom RSQLite dbConnect dbDisconnect dbExistsTable dbListTables dbReadTable dbRemoveTable dbWriteTable SQLite
#' @importFrom generics tidy glance augment
#' @importFrom graphics abline legend par points title
#' @importFrom sf st_as_sf st_bbox st_centroid st_coordinates st_crs st_delete st_drop_geometry
#'   st_geometry_type st_intersects st_read st_set_crs st_set_geometry st_write
#' @importFrom spmodel AICc AUROC conditional covmatrix dispersion_initial dispersion_params glances kcv loocv pseudoR2 randcov_initial randcov_params satterthwaite varcomp
#' @importFrom stats aggregate AIC anova BIC coef coefficients complete.cases confint cooks.distance cor dbeta
#'   delete.response dbinom deviance dist dgamma dnbinom dpois fitted fitted.values formula .getXlevels hatvalues
#'   influence kmeans lm logLik model.frame model.matrix model.offset model.response na.fail na.omit na.pass nobs pchisq
#'   pf pnorm predict printCoefmat pt qnorm qqnorm qqline qt quantile rbinom resid residuals
#'   reformulate rbeta rgamma rnbinom rnorm rpois rstandard sd terms uniroot var vcov
#' @importFrom tibble tibble as_tibble
#' @importFrom utils combn read.table tail
#' @importFrom withr local_dir
#' @importFrom parallel clusterApply clusterEvalQ detectCores makeCluster parLapply parLapply stopCluster
#' @import doParallel
#' @importFrom doParallel registerDoParallel
#' @import foreach
#' @import itertools
#' @import iterators
#' @importFrom filematrix fm.create fm.load fm.open close
NULL

## @import doParallel registerDoParallel
## #' @importFrom parallel detectCores makeCluster parLapply stopCluster
