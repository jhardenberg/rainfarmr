#' Spectral merging of a coarse field and of a fine field at a given wavenumber
#'
#' @description The Fourier spectra of the input fields are merged at wavenumber `kmax`.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param ra matrix containing a coarse field of size `c(nas, nas)`.
#' @param r matrix containing a field at higher resolution.
#' @param kmax wavenumber to use for merging (default `nas/2`).
#' @return The merged field, in physical space.
#' @export
#' @examples 
#' # Make a coarse field with power-law Fourier spectrum 
#' fa <- initmetagauss(1.7, 8)
#' ra <- metagauss(fa)
#' # Make a fine power-law Fourier spectrum
#' f <- initmetagauss(1.7, 32)
#' r <- metagauss(f)
#' # Merge the two fields in spectral space
#' rm <- mergespec(ra, r, kmax = 4)
#' # Check spectral slope of the resulting field
#' fx <- fft2d(rm)
#' fitslope(fx)
#' # 1.678554


mergespec <- function(ra, r, kmax = 0) {
  nx <- dim(r)[1]
  nax <- dim(ra)[1]
  nax2 <- nax / 2
  if (kmax == 0) {
    kmax <- nax2
  }

  pstr <- fft2d(r) * (nx * nx) ^ 2
  pstra <- fft2d(ra) * (nax * nax) ^ 2
  c <- pstra[kmax] / pstr[kmax]
  r <- r * sqrt(c)

  DFTf <- fft(r)
  DFTr <- fft(ra)

  # fftshift
  DFTr2 <- matrix(0., nx, nx)
  DFTr2[1:nax2, 1:nax2] <- DFTr[1:nax2, 1:nax2]
  DFTr2[nx - nax2 + (1:nax2), 1:nax2] <- DFTr[nax2 + (1:nax2), 1:nax2]
  DFTr2[1:nax2, nx - nax2 + (1:nax2)] <- DFTr[1:nax2, nax2 + (1:nax2)]
  DFTr2[nx - nax2 + (1:nax2), nx - nax2 + (1:nax2)] <- DFTr[nax2 + (1:nax2),
                                                            nax2 + (1:nax2)]
  DFTr2[nax2 + 1, 1] <- Conj(DFTr2[nx + 1 - nax2, 1 ])
  DFTr2[1, nax2 + 1] <- Conj(DFTr2[1, nx + 1 - nax2])

  kx <- rep(c(0:(nx / 2), (-nx / 2 + 1):-1), nx)
  dim(kx) <- c(nx, nx)
  k2 <- kx ^ 2 + t(kx) ^ 2  # t(kx) is ky
  kmax2 <- kmax ^ 2
  ddx <- 2 * pi / nax / 2 - 2 * pi / nx / 2

  DFTf <- DFTf * (k2 > kmax2) + DFTr2 * (k2 <= kmax2) * exp(-1i * ddx * kx -
                                                            1i * ddx * t(kx))
  r <- Re(fft(DFTf, inverse = TRUE)) / length(DFTf)
  return(r)
}
