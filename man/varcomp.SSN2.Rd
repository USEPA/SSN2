% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varcomp.R
\name{varcomp.SSN2}
\alias{varcomp.SSN2}
\alias{varcomp.ssn_lm}
\alias{varcomp.ssn_glm}
\title{Variability component comparison}
\usage{
\method{varcomp}{ssn_lm}(object, ...)

\method{varcomp}{ssn_glm}(object, ...)
}
\arguments{
\item{object}{A fitted model object from \code{\link[=ssn_lm]{ssn_lm()}} or \code{\link[=ssn_glm]{ssn_glm()}}.}

\item{...}{Other arguments. Not used (needed for generic consistency).}
}
\value{
A tibble that partitions the the total variability by the fixed effects
and each variance parameter. The proportion of variability explained by the
fixed effects is the pseudo R-squared obtained by \code{psuedoR2()}. The
remaining proportion is spread accordingly among each variance parameter:
\code{"tailup_de"}, \code{"taildown_de"}, \code{"euclid_de"}, \code{"nugget"},
and if random effects are used, each named random effect. For \code{ssn_glm()},
models, only the variances on the link scale are considered (i.e., the variance
function of the response is omitted).
}
\description{
Compare the proportion of total variability explained by the fixed effects
and each variance parameter.
}
\examples{
# Copy the mf04p .ssn data to a local directory and read it into R
# When modeling with your .ssn object, you will load it using the relevant
# path to the .ssn data on your machine
copy_lsn_to_temp()
temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(temp_path, overwrite = TRUE)

ssn_mod <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM,
  ssn.object = mf04p,
  tailup_type = "exponential",
  additive = "afvArea"
)
varcomp(ssn_mod)
}
