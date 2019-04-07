#' Generate the spectral amplitudes for a metagaussian field
#'
#' @description A Fourier spectrum with prescribed logarithimic slope and
#' zero phases is constructed.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param sx spectral slope for the output field.
#' The convention is that this is the slope of `k*|A(k)|^2`.
#' @param ns size of the output field.
#' @return Output complex field in Fourier space with specified spectral slope,
#'  with dimensions `c(ns, ns)`.
#' @export
#' @examples
#' # Make a synthetic rainfall field with prescribed logarithmic spectral slope
#' f = initmetagauss(1.7, 64)
#' r = metagauss(f)
#' # Check slope of the resulting field
#' fx <- fft2d(r)
#' fitslope(fx)
#' # 1.640373
initmetagauss <- function(sx, ns) {
    sx <- abs(sx)
    kx <- rep(c(0:(ns / 2), (- ns / 2 + 1):-1), ns)
    dim(kx) <- c(ns, ns)
    kx <- kx ^ 2 + t(kx) ^ 2
    kx[1, 1] <- 1
    kx <- (kx ^ (- (sx + 1) / 4))
    kx[1, 1] <- 0
    kx <- kx / sqrt(sum(abs(kx) ^ 2)) * ns * ns
    return(kx)
}
