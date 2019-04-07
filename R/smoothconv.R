#' Smoothening using convolution with a circular kernel
#'
#' @description The input field is convolved with a circular kernel with equal
#' weights. Takes into account missing values.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param z matrix with the input field to smoothen, with dimensions `c(ns, ns)`
#' @param nas the smoothing kernel uses a radius `(ns/nas)/2`
#' @return The smoothened field.
#' @export
#' @examples
#' z <- rnorm(64 * 64)
#' dim(z) <- c(64, 64)
#' zs <- smoothconv(z, 8)
#' sd(zs)
#' # [1] 0.07910996
smoothconv <- function(z, nas) {
  imask <- !is.finite(z)
  z[imask] <- 0.
  ns <- dim(z)[1]

  sdim <- (ns / nas) / 2  # the smoothing sigma has width half a large pixel
  mask <- matrix(0., ns, ns)
  for (i in 1:ns) {
    for (j in 1:ns) {
      kx <- i - 1
      ky <- j - 1
      if (i > ns / 2 + 1) {
        kx <- i - ns - 1
      }
      if (j > ns / 2 + 1) {
        ky <- j - ns - 1
      }
      r2 <- kx * kx + ky * ky
      mask[i, j] <- exp( - (r2 / (sdim * sdim)) / 2)
    }
  }
  fm <- fft(mask)
  zf <- Re(fft(fm * fft(z), inverse = TRUE)) / sum(mask) / length(fm)
  if (sum(imask) > 0) {
           z[!imask] <- 1.0
           zf <- zf / (Re(fft(fm * fft(z), inverse = TRUE)) /
                       sum(mask) / length(fm))
  }
  zf[imask] <- NA
  return(zf)
}
