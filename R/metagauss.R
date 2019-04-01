#' Generate a metagaussian field multiplying the spectral field `f`
#' with random phases and performing an inverse FFT transform to real space.
#' @param f Input field generated with the `initmetagauss()`function.
#' @return The input field with random Fourier phases
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
