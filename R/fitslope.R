#' Compute logarithmic slope of a spatial spectrum
#' @description Uses a linear fit to derive the log-log slope of a Fourier spectrum. 
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param fx vector containing input spectrum starting from `k=1`.
#' @param kmin minimum wavenumber for logarithmic fit range.
#' @param kmax maximum wavenumber for logarithmic fit range.
#' @return The spatial spectral slope minus one. The slope is interpreted as the logarithmic slope of k*|A(k)|^2 where A(k) are the spectral amplitudes of the input field.
#' @export
#' @examples
#' f = initmetagauss(1.7, 64)
#' z = metagauss(f)
#' fx <- fft2d(z)
#' sx <- fitslope(fx)
#' print(sx)
#' # 1.640373
fitslope <- function(fx, kmin = 1, kmax = 0) {
  if (kmax == 0) {
    kmax=length(fx)
  } 
  k <- 1:length(fx)
  k <- k[kmin:length(k)]
  rel <- lm(log(fx[kmin:length(k)]) ~ log(k))
  sx <- rel$coefficients[2]
  return(abs(sx) - 1)
}
