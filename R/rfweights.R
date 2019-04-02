#' Derive weights from a fine-scale precipitation climatology
#' @description Weights for downscaling are computed interpolating a fine-scale
#' climatology to the target grid and dividing it by an averaged version of itself.
#' A suitable climatology could be represented for example by a fine-scale
#' precipitation climatology from a high-resolution regional climate model
#' (see e.g. Terzago et al. 2018), a local high-resolution gridded climatology
#' from observations, or a reconstruction such as those which can be downloaded
#' from the WORLDCLIM (http://www.worldclim.org) or CHELSA (http://chelsa-climate.org)
#' websites. The latter data will need to be converted to NetCDF format before
#' being used (see for example the GDAL tools (https://www.gdal.org).
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @references Terzago, S. et al. (2018). NHESS 18(11), 2825–2840.
#' @param z matrix with a spatial field of fine-scale precipitation climatology.
#' @param lon vector of longitudes of the high-resolution climatology
#' @param lat vector of latitudes of the high-resolution climatology
#' @param lonc vector of longitudes of the coarse field to downscale
#' @param latc vector of latitudes of the coarse field to downscale
#' @param nf refinement factor for downscaling (the coarse resolution is increased
#' by this factor).
#' The number of longitudes and latitudes is expected to be equal.
#' @param fsmooth logical to compute weights against a smooth average.
#' If false box averaging is used.
#' @return The matrix of weights with dimensions `c(lonc*nf, latc*nf)`
#' @export
rfweights <- function(z, lon, lat, lonc, latc, nf, fsmooth=TRUE) {

  if (length(lonc) != length(latc)) {
    stop(paste("The target grid must have the same",
               "number of longitudes and latitudes."))
  }
  nas <- length(lonc)
  grid <- lon_lat_fine(lonc, latc, nf)
  zr <- remapcon(lon, lat, z, grid$lon, grid$lat)
  if (fsmooth) {
    ww <- zr / smoothconv(zr, nas);
  } else {
    ww <- zr / interpola(agg(zr, nas), nas * nf)
  }
  return(ww)
}
