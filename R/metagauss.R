#' Generate a metagaussian field
#'
#' @description Random Fourier phases are added to the input spectrum `f` and 
#' an inverse FFT transform to real space is performed.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param f matrix with complex spectral amplitudes generated with the
#' [initmetagauss()] function.
#' @return A metagaussian field with random Fourier phases.
#' @export
#' @examples
#' f = initmetagauss(1.7, 64)
#' z = metagauss(f)
#' sd(z)
#' # [1] 1.000122
metagauss <- function(f) {
    ns <- dim(f)[1]
    ph <- matrix(rnorm(ns * ns), ns, ns)
    ph <- fft(ph)
    ph <- ph / abs(ph)
    ph[1, 1] <- 0.
    ph <- f * ph
    ph <- Re(fft(ph, inverse = TRUE)) / (ns * ns)
    return(ph)
}
