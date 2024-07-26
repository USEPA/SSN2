#' Read relevant sf objects from user shalpefiles or geopackages
#'
#' @param fn Path to shapefile or geopackage
#'
#' @noRd
get_sf_obj <- function(fn) {
  f.ext.shp <- substr(fn, nchar(fn) - 3, nchar(fn)) == ".shp"
  f.ext.gpkg <- substr(fn, nchar(fn) - 4, nchar(fn)) == ".gpkg"

  ## Check if shapefile or geopackage exists----------
  ## If shapefile specified, check if it exists
  if (f.ext.shp) {
    shp.exists <- file.exists(fn)
    gpkg.exists <- FALSE

    fn <- substr(fn, 1, nchar(fn) - 4)
  }
  ## If geopackage specified, check if it exists
  if (f.ext.gpkg) {
    shp.exists <- FALSE
    gpkg.exists <- file.exists(fn)

    fn <- substr(fn, 1, nchar(fn) - 5)
  }

  ## If no file extension specified, check for both
  if (!f.ext.shp & !f.ext.gpkg) {
    shp.exists <- file.exists(paste0(fn, ".shp"))
    gpkg.exists <- file.exists(paste0(fn, ".gpkg"))
  }

  ## If both exist without file extension, return error
  if (gpkg.exists & shp.exists) {
    ## sfobj<- st_read(paste0(fn, ".gpkg"), quiet = TRUE)
    stop(paste0(
      dirname(fn), " contains both ", basename(fn), ".shp and ",
      basename(fn), ".gpkg."
    ))
  }
  ## Read in shapefile
  if (shp.exists & !gpkg.exists) {
    sfobj <- st_read(paste0(fn, ".shp"), quiet = TRUE)
  }
  ## Read in geopackage
  if (gpkg.exists & !shp.exists) {
    sfobj <- st_read(paste0(fn, ".gpkg"), quiet = TRUE)
  }
  ## Neither exists
  if (!gpkg.exists & !shp.exists) {
    stop(paste0(
      basename(fn), " not found in ", dirname(fn), ". ",
      basename(fn), " must reside in path and be in .shp or .gpkg format."
    ))
  }

  return(sfobj)
}
