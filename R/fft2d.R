#' Compute spatial Fourier power spectrum
#'
#' @description The Fourier power spectrum of the input field is computed and
#' averaged over shells (bins with width 1) of the modulus of the wavenumber.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param z matrix or array of input data with spatial dimensions `c(N, N)`.
#' If `z` is an array the power spectra are averaged over the third dimension.
#' The first two dimensions must be the same (the input fields must be square).
#' @return Spectral power (average square of absolute value of spectral
#' amplitudes) for wavenumbers `k=1:(N/2)`.
#' @export
#' @examples
#' # Make a synthetic rainfall field with prescribed logarithmic spectral slope
#' f = initmetagauss(1.7, 64)
#' r = metagauss(f)
#' # Check slope of the resulting field
#' fx <- fft2d(r)
#' fitslope(fx)
#' # 1.640373
fft2d <- function(z) {
  ns <- dim(z)[1]
  ny <- dim(z)[2]
  if (ns != ny) {
    stop("The input matrix must be square.")
  }
  ns2 <- floor(ns / 2)
  nt <- dim(z)[3]
  if (is.na(nt)) {
    dim(z) <- c(ns, ns, 1)
    nt <- 1
  }
  kx <- rep(c(0:(ns / 2), (-ns / 2 + 1):-1), ns)
  dim(kx) <- c(ns, ns)
  km <- sqrt(kx ^ 2 + t(kx) ^ 2)
  nn <- integer(ns)
  fx <- numeric(ns)
  for (it in 1:nt) {
    zf <- abs(fft(z[ , , it]) / (ns * ns)) ^ 2
    zf0 <- zf
    zf[ns2 + 1, ] <- zf0[ns2 + 1, ] / 2
    zf[ , ns2 + 1] <- zf0[ , ns2 + 1] / 2
    for (i in 1:(ns * ns)) {
      ik <- floor(km[i] + 1.5)
      fx[ik] <- fx[ik] + zf[i]
      nn[ik] <- nn[ik] + 1
    }
  }
  fx <- fx[2:(ns2 + 1)] / nn[2:(ns2 + 1)];
  return(fx)
}
