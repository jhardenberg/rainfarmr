#' Compute logarithmic slope of a spatial power spectrum
#'
#' @description Uses a linear fit to derive the log-log slope of a
#' Fourier power spectrum. 
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param fx vector containing input power spectrum starting from `k=1`.
#' @param kmin minimum wavenumber for logarithmic fit range.
#' @param kmax maximum wavenumber for logarithmic fit range.
#' @return The spatial spectral slope minus one. The slope is returned as the
#' logarithmic slope of k*|A(k)|^2 where |A(k)|^2 are the squared spectral
#' amplitudes provided in input..
#' @export
#' @examples
#' # Make a synthetic rainfall field with prescribed logarithmic spectral slope
#' f = initmetagauss(1.7, 64)
#' r = metagauss(f)
#' # Check spectral slope of the resulting field
#' fx <- fft2d(r)
#' print(fitslope(fx))
#' # 1.640373
fitslope <- function(fx, kmin = 1, kmax = length(fx)) {
  k <- 1:length(fx)
  rel <- lm(log(fx[kmin:kmax]) ~ log(k[kmin:kmax]))
  sx <- rel$coefficients[2]
  return(abs(sx) - 1)
}
