#' MiddleFork04.ssn: Middle Fork 2004 stream temperature dataset
#'
#' The \code{MiddleFork04.ssn} data folder contains the spatial, attribute,
#' and topological information needed to construct an SSN object using
#' the \code{SSN2} package.
#'
#' The \code{MiddleFork04.ssn} folder contains five shapefiles:
#'
#' \itemize{
#'    \item edges: polyline shapefile representing the stream network
#'    \item sites: point shapefile representing the observed site locations
#'    \item pred1km: point shapefile representing prediction site locations at
#'      approximately 1km intervals throughout the stream network
#'    \item Knapp: point shapefile representing prediction site locations on the Knapp River
#'    \item CapeHorn: point shapefile representing prediction site locations on the Cape Horn River
#' }
#'
#' The \code{MiddleFork04.ssn} includes one text file, \code{netID1.txt}, which contains the
#' topological information for the stream network in the Middle Fork 2004
#' dataset.
#'
#' The distance folder contains four folders that store the hydrologic
#' distance matrices for each of the point shapefiles (\code{obs}, \code{CapeHorn},
#' \code{Knapp}, and \code{pred1km}). See [ssn_create_distmat()] for a
#' detailed description of the distance matrix file structure.
#'
#' Attribute data is also stored within each of the spatial
#' datasets. The column names are defined as follows:
#'
#'  \code{edges}:
#'   \itemize{
#'     \item COMID: Common identifier of an NHD feature or relationship
#'     \item GNIS_Name: Feature name as found in the Geographic Names Information System
#'     \item REACHCODE: Unique identifier for a reach. The first 8 digits contain the identfier for the HUC8 and the last 6 digits are a unique within-HUC8 identifier for the reach
#'     \item FTYPE: three-digit integer used to classify hydrography features in the NHD and define subtypes
#'     \item FCODE: Numeric code that contains the feature type and its attributes as found in the NHDFCode lookup table
#'     \item CDRAINAG: Cumulative drainage area (km2) for the lowermost location on the edge
#'     \item AREAWTMAP: Area weighted mean annual precipitation (mm) at the lowermost location on the edge
#'     \item SLOPE: Slope of the edge (cm/cm)
#'     \item h2oAreaKm2: Watershed area (km2) for the lowermost location on the line segment
#'     \item rid: Reach identifier
#'     \item areaPI: Segment proportional influence value, calculated using watershed area (h2oAreaKm2)
#'     \item afvArea: Additive function value, calculated using areaPI
#'     \item upDist: Distance from the stream outlet (most downstream location in the the stream network) to the uppermost location on the line segment
#'     \item Length: Length of line segment (m)
#'     \item netID: Network identifier
#' }
#'   \code{sites}:
#'   \itemize{
#'     \item STREAMNAME: Stream name
#'     \item COMID: Common identifier of an NHD feature or relationship
#'     \item CDRAINAG: Cumulative drainage area (km2)
#'     \item AREAWTMAP: Area weighted mean annual precipitation (mm) at lowermost
#'       location on the line segment where the site resides
#'     \item SLOPE: Slope of the line segment (cm/cm) where the site resides
#'     \item ELEV_DEM: Elevation at the site based on a 30m DEM
#'     \item Source: Source of the data - relates to the ID field of the source table
#'     \item Summer_mn: Overall summer mean termperature (C) of the deployment
#'     \item MaxOver20: Binary variable: 1 represents the maximum summer temperature
#'       was greater than 20C and 0 indicates that it was less than 20C
#'     \item C16: Number of times daily stream temperature exceeded 16C
#'     \item C20: Number of times daily stream temperature exceeded 20C
#'     \item C24: Number of times daily stream temperature exceeded 24C
#'     \item FlowCMS: Average stream flow (cubic meters per sec) for August,
#'       by year, from 1950-2010 across 9 USGS gauges in the region
#'     \item AirMEANc: Average mean air temperature (C) from July 15 - August 31,
#'       from 1980-2009 across 10 COOP air stations within the domain
#'     \item AirMWMTc: Average maximum air temperature (C) from July 15 - August 31,
#'       from 1980-2009 across 10 COOP air stations within the domain.
#'       MWMT = maximum 7-day moving average of the maximum daily temperature
#'       (i.e. maximum of all the 7-day maximums)
#'     \item NEAR_X: x coordinate
#'     \item NEAR_Y: y coordinate
#'     \item rid: Reach identifier of the edge the site resides on
#'     \item ratio: Site ratio value; provides the proportional distance along the edge to the site location
#'     \item upDist: Distance upstream from the stream outlet (m)
#'     \item afvArea: Additive function value calculated using waterhsed area (h2oAreaKm2)
#'     \item locID: Location identifier
#'     \item netID: Stream network identifier
#'     \item pid: Point identifier
#' }
#'
#'   \code{pred1km}, \code{CapeHorn}, and \code{Knapp}:
#'   \itemize{
#'     \item COMID: Common identifier of an NHD feature or relationship
#'     \item GNIS_Name: Feature name of the edge the site resides on, as found in the Geographic Names Information System
#'     \item CDRAINAG: Cumulative drainage area (km2)
#'     \item AREAWTMAP: Area weighted mean annual precipitation (mm) at lowermost location on the line segment where the site resides
#'     \item SLOPE: Slope of the line segment (cm/cm) where the site resides
#'     \item ELEV_DEM: Elevation at the site based on a 30m DEM
#'     \item NEAR_X: x coordinate
#'     \item NEAR_Y: y coordinate
#'     \item rid: Reach identifier of the edge the site resides on
#'     \item ratio: Site ratio value; provides the proportional distance along the edge to the site location
#'     \item upDist: Distance upstream from the stream outlet (m)
#'     \item afvArea: Additive function value calculated using watershed area (h2oAreaKm2)
#'     \item locID: Location identifier
#'     \item netID: Stream network identifier
#'     \item pid: Point identifier
#'     \item FlowCMS: Average stream flow (cubic meters per sec) for August, by year, from 1950-2010 across 9 USGS gauges in the region
#'     \item AirMEANc: Average mean air temperature (C) from July 15 - August 31, from 1980-2009 across 10 COOP air stations within the domain
#'     \item AirMWMTc: Average maximum air temperature (C) from July 15 - August 31, from 1980-2009 across 10 COOP air stations within the domain. MWMT = maximum 7-day moving average of the maximum daily temperature(i.e. maximum of all the 7-day maximums)
#'  }
#'
#' @source \code{edges} are a modified version of the United States
#'   National Hydrography Dataset
#'   (http://nhd.usgs.gov/). \code{sites}, \code{pred1km}, \code{CapeHorn}
#'   and \code{Knapp} are unpublished United States Forest Service data.
#'
#' @docType data
#'
#' @name MiddleFork04.ssn
#'
#' @seealso [mf04p] for the Middle For 04 data as an \code{SSN} object.
NULL


#' Imported SSN object from the MiddleFork04.ssn data folder
#'
#'  The MiddleFork04.ssn data folder contains the spatial, attribute,
#'  and topological information needed to construct a spatial stream
#'  network object using the SSN2 package. \code{mf04p} was created
#'  using [ssn_import()].
#'
#' @seealso
#' [MiddleFork04.ssn] for details about the contents of \code{mf04p}.
#'   [ssn_import()] to convert a .ssn object to an \code{SSN} object in R.
#'   [ssn_create_distmat] for details about the distance matrix file structure.
"mf04p"
