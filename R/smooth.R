#' Smoothening using convolution with a circular kernel
#'
#' @description The input field is convolved with a circular kernel with equal
#' weights. Takes into account missing values.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param z matrix with the input field to smoothen, with dimensions `c(ns, ns)`
#' @param sdim the smoothing kernel radius in pixel
#' @return The smoothened field.
#' @export
#' @examples
#' z <- rnorm(64 * 64)
#' dim(z) <- c(64, 64)
#' zs <- smooth(z, 8)
#' sd(zs)
#' # [1] 0.1334648

smooth <- function(z, sdim) {
  imask <- !is.finite(z)
  z[imask] <- 0.
  nsx <- dim(z)[1]
  nsy <- dim(z)[2]

  mask <- matrix(0., nsx, nsy)
  for (i in 1:nsx) {
    for (j in 1:nsy) {
      kx <- i - 1
      ky <- j - 1
      if (i > nsx / 2 + 1) {
        kx <- i - nsx - 1
      }
      if (j > nsy / 2 + 1) {
        ky <- j - nsy - 1
      }
      r2 <- kx * kx + ky * ky
      #mask[i, j] <- exp( - (r2 / (sdim * sdim)) / 2)
      mask[i, j] <- (r2 <= (sdim * sdim))
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
