#' Smoothen field with a circular kernel using convolution.
#' Takes into account missing values.
#' @param zi The input field to smoothen, with dimensions `c(ns, ns)`
#' @param nas The smoothing kernel uses a radius `(ns/nas)/2`
#' @return The smoothened field.
#' @export
smoothconv <- function(zi, nas) {
  imask <- !is.finite(zi)
  zi[imask] <- 0.
  ns <- dim(zi)[1]

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
  zf <- Re(fft(fm * fft(zi), inverse = TRUE)) / sum(mask) / length(fm)
  if (sum(imask) > 0) {
           zi[!imask] <- 1.0
           zf <- zf / (Re(fft(fm * fft(zi), inverse = TRUE)) /
                       sum(mask) / length(fm))
  }
  zf[imask] <- NA
  return(zf)
}
