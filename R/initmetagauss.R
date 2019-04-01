#' Generate the spectral amplitudes for a metagaussian field
#' @param sx Spectral slope for the output field.
#' The convention is that this is the slope of `k*|A(k)|^2`
#' @param ns Size of the output field
#' @return Output field with specified spectral slope
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
