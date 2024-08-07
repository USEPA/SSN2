% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssn_create_distmat.R
\name{ssn_create_distmat}
\alias{ssn_create_distmat}
\title{Calculate Hydrologic Distances for an \code{SSN} object}
\usage{
ssn_create_distmat(
  ssn.object,
  predpts = NULL,
  overwrite = FALSE,
  among_predpts = FALSE,
  only_predpts = FALSE
)
}
\arguments{
\item{ssn.object}{An \code{SSN} object}

\item{predpts}{name of prediction points in an \code{SSN} object.
When a vector with length greater than one, each name is iterated upon.
Default is NULL.}

\item{overwrite}{Logical. If \code{TRUE}, overwrite existing distance
matrices. Defaults to \code{FALSE}.}

\item{among_predpts}{Logical. If \code{TRUE}, compute the pairwise distances
between the prediction sites. Defaults to \code{FALSE}.}

\item{only_predpts}{Logical. If \code{TRUE}, only compute distances for
prediction sites. Defaults to \code{FALSE}.}
}
\value{
The \command{ssn_create_distmat} function creates a collection
of hierarchical directories in the \code{ssn$path} directory,
which store the pairwise distances between sites associated with
the \code{SSN} object. See details section for additional information.
}
\description{
Creates a collection of (non-symmetric) matrices
containing pairwise downstream hydrologic distances between sites
in an \code{SSN} object
}
\details{
A distance matrix that contains the hydrologic distance
between any two sites in \code{SSN} object is needed to fit a spatial
statistical model using the tail-up and tail-down autocovariance
functions described in Ver Hoef and Peterson (2010). These models
are implemented in R via \command{ssn_lm} and \command{ssn_glm} in
the\code{SSN2} package. The hydrologic distance information needed to
model the covariance between flow-connected (i.e. water flows
from one location to the other) and flow-unconnected (i.e. water
does not flow from one location to the other, but they reside on
the same network) locations differs. The total hydrologic
distance is a directionless measure; it represents the hydrologic
distance between two sites, ignoring flow direction. The
hydrologic distance from each site to a common downstream stream
junction is used when creating models for flow-unconnected pairs,
which we term downstream hydrologic distance. In contrast, the
total hydrologic distance is used for modeling flow-connected
pairs, which we term total hydrologic distance.

A downstream hydrologic distance matrix provides enough
information to meet the data requirements for both the tail-up and
tail-down models. When two locations are flow-connected, the
downstream hydrologic distance from the upstream location to the
downstream location is greater than zero, but it is zero in the
other direction. When two locations are flow-unconnected the
downstream hydrologic distance will be greater than zero in both
directions. A site's downstream hydrologic distance to itself is
equal to zero. The format of the downstream hydrologic distance
matrix is efficient because distance information needed to fit
both the tail-up and tail-down models is only stored once. As an
example, a matrix containing the total hydrologic distance between
sites is easily calculated by adding the downstream distance
matrix to its transpose.

The downstream hydrologic distances are calculated based on the
binaryIDs and stored as matrices. The matrices are stored in a
directory named \sQuote{distance}, which is created by the
\command{ssn_create_distmat} function within the .ssn directory. The distance
directory will always contain at least one directory named
\sQuote{obs}, which contains a number of .RData files, one for each
network that has observed sites residing on it. The naming
convention for the files is based on the netID number
(e.g. dist.net1.RData). Each matrix in the \sQuote{obs} folder
contains the information to form a square matrix, which contains
the downstream hydrologic distance between each pair of observed
sites on the network. Direction is preserved, with columns
representing the FROM site and rows representing the TO site. Row
and column names correspond to the pid attribute for each site.

If the argument \code{predpts} is specified in the call to the
function, the downstream hydrologic distances between the observed
and prediction sites will also be computed. A new directory is
created within the distance directory, with the name corresponding
to the names attribute for the preds
(e.g. \code{attributes(ssn.object$preds)$names}). A sequence of
.RData files is created within this directory, similar to the
structure for the observed sites, except that two objects are
stored for each network that contains \emph{both} observed and
prediction sites. The letters \code{a} and \code{b} are used in
the naming convention to distinguish between the two objects
(e.g. dist.net1.a and dist.net1.b). The matrices that these
objects represent are not necessarily square. In matrices of type
\code{a}, rows correspond to observed locations and columns to
prediction locations. In contrast, rows correspond to prediction
locations and columns to observed locations in matrices of type
\code{b}. Direction is also preserved, with columns representing
the FROM site and rows representing the TO site in both object
types. Again, row and column names correspond to the pid attribute
for each site.

If \code{among_predpts = TRUE}, the downstream
hydrologic distances will also be computed between prediction
sites, for each network. Again these are stored within the distance
directory with the name corresponding to the prediction points
dataset. The naming convention for these prediction to prediction
site distance matrices is the same as the distance matrices stored
in the \sQuote{obs} directory (e.g. dist.net1.RData). These extra
distance matrices are needed to perform block Kriging using
\code{\link[SSN2]{predict.ssn_lm}}.

If \code{only_predpts = TRUE}, the downstream
hydrologic distances will not be calculated between observed sites
themselves. Pairwise distances will only be calculated for observed
and prediction locations and. Pairwise distances between prediction
locations will also be calculated if \code{among_predpts = TRUE}.
}
\examples{
## Copy the MiddleForke04.ssn data to a local temporary directory.
## Only needed for this example.
copy_lsn_to_temp()
## Import SSN data
mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
  predpts = c("pred1km.gpkg", "CapeHorn"),
  overwrite = TRUE
)

## Create distance matrices for observations and one set of prediction sites
## Include hydrologic distance matrices among prediction sites.
ssn_create_distmat(mf04p,
  predpts = "pred1km", overwrite = TRUE,
  among_predpts = TRUE
)

## Create distance matrices for an additional set of prediction points.
## Distance matrices for observations and pred1km prediction sites are
## not recalculated.
ssn_create_distmat(mf04p,
  predpts = "CapeHorn", overwrite = TRUE,
  among_predpts = TRUE, only_predpts = TRUE
)
}
