#' Spectral merging of coarse field `ra` and fine field `r` at wavenumber `kmax`.
#' @param ra The coarse field of size `c(nas, nas)`
#' @param r The fine field
#' @param kmax Maximum wavenumber to use (default nas/2)
#' @return The merged field
mergespec <- function(ra, r, kmax = 0) {
  nx <- dim(r)[1]
  nax <- dim(ra)[1]
  nax2 <- nax / 2
  if (kmax == 0) {
    kmax <- nax2
  }
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
