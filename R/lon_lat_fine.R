#' Interpolate longitude and latitude vectors to higher resolution by a factor `nf`.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param lon vector of longitudes.
#' @param lat vector of latitudes.
#' @param nf factor by which to increase resolution.
#' @return List with elements `lon` and `lat` with resolution increased by factor `nf`.
#' @export
#' @import utils
#' @examples
#' lon <- 5:9
#' lat <- 43:47
#' nf <- 4
#' grid <- lon_lat_fine(lon, lat, nf)
#' grid$lon
#' #  [1] 4.625 4.875 5.125 5.375 5.625 5.875 6.125 6.375 6.625 6.875 7.125 7.375
#' # [13] 7.625 7.875 8.125 8.375 8.625 8.875 9.125 9.375
#' grid$lat
#' #  [1] 42.625 42.875 43.125 43.375 43.625 43.875 44.125 44.375 44.625 44.875
#' # [11] 45.125 45.375 45.625 45.875 46.125 46.375 46.625 46.875 47.125 47.375
lon_lat_fine <- function(lon, lat, nf) {

  nlon <- length(lon)
  nlat <- length(lat)
  lone <- c(2 * lon[1] - lon[2], lon, 2 * tail(lon, 1) - tail(lon, 2)[1])
  late <- c(2 * lat[1] - lat[2], lat, 2 * tail(lat, 1) - tail(lat, 2)[1])
  xc <- 0:(nlon + 1)
  yc <- 0:(nlat + 1)
  xf <- seq(0.5 + 1 / (2 * nf), nlon + 0.5 - 1 / (2 * nf), 1 / nf)
  yf <- seq(0.5 + 1 / (2 * nf), nlat + 0.5 - 1 / (2 * nf), 1 / nf)
  lonf <- approx(xc, lone, xf, method = "linear")
  latf <- approx(yc, late, yf, method = "linear")

  return(list(lon = lonf$y, lat = latf$y))
}
