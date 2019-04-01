#' Interpolate longitude and latitude vectors to higher resolution by a factor `nf`.
#' @param lon Vector of longitudes
#' @param lon Vector of latitudes
#' @return List with elements `lon` and `lat` at finer resolution.
#' @export
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
