#' Compute logarithmic slope of a spatial spectrum
#' @param sx Vector containing input spectrum starting from `k=1`
#' @param kmin Wavenumber from which to start fit
#' @return The spatial spectral slope minus one
#' @export
#' f = initmetagauss(1.7, 64)
#' z = metagauss(f)
#' fx <- fft2d(z)
#' sx <- fitslope(fx)
#' print(sx)
#' # 1.640373
fitslope <- function(fx, kmin = 1) {
  k <- 1:length(fx)
  k <- k[kmin:length(k)]
  rel <- lm(log(fx[kmin:length(k)]) ~ log(k))
  sx <- rel$coefficients[2]
  return(abs(sx) - 1)
}
