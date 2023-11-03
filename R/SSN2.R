#' @keywords internal
"_PACKAGE"

#' @useDynLib SSN2
#' @importFrom Matrix bdiag chol chol2inv colSums crossprod determinant diag Diagonal
#'   forceSymmetric Matrix print qr rankMatrix solve sparseMatrix t tcrossprod which
#' @importFrom RSQLite dbConnect dbDisconnect dbExistsTable dbListTables dbReadTable dbRemoveTable dbWriteTable SQLite
#' @importFrom generics tidy glance augment
#' @importFrom graphics abline legend par points title
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @importFrom sf st_as_sf st_bbox st_centroid st_coordinates st_crs st_delete st_drop_geometry
#'   st_geometry_type st_intersects st_read st_set_crs st_set_geometry st_write
#' @importFrom spmodel AICc covmatrix dispersion_initial dispersion_params glances loocv pseudoR2 randcov_initial randcov_params varcomp
#' @importFrom stats AIC anova coef coefficients complete.cases confint cooks.distance cor dbeta
#'   delete.response dbinom deviance dist dgamma dnbinom dpois fitted fitted.values formula .getXlevels hatvalues
#'   influence kmeans lm logLik model.frame model.matrix model.offset model.response na.omit na.pass pchisq
#'   pnorm predict printCoefmat pt qnorm qqnorm qqline qt quantile rbinom resid residuals
#'   reformulate rbeta rgamma rnbinom rnorm rpois rstandard terms var vcov
#' @importFrom tibble tibble as_tibble
#' @importFrom utils read.table tail
NULL
