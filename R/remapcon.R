#' Conservative remapping
#'
#' @description Implements conservative remapping weighting with the overlap area between pixels.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param x vector of input longitudes
#' @param y vector of input latitudes
#' @param z matrix of input data
#' @param xo vector of target longitudes
#' @param yo vector of target latitudes
#' @return A remapped matrix of dimensions `c(length(xo), length(yo))`
#' @export
#' @examples
#' z <- 1:(31*51)
#' dim(z) <- c(31, 51)
#' x <- seq(4, 10, 0.2)
#' y <- seq(30, 40, 0.2)
#' xo <- seq(5, 6, 0.5)
#' yo <- seq(35, 37, 0.5)
#' zo <- remapcon(x, y, z, xo, yo)
#' zo
#' #       [,1]  [,2]  [,3]   [,4]   [,5]
#' # [1,] 781.0 858.5 936.0 1013.5 1091.0
#' # [2,] 783.5 861.0 938.5 1016.0 1093.5
#' # [3,] 786.0 863.5 941.0 1018.5 1096.0
remapcon <- function(x, y, z, xo, yo) {
  nx <- length(x)
  nxo <- length(xo)

  xoe <- numeric(nxo + 2)
  xoe[2:(nxo + 1)] <- xo
  xoe[1] <- 2 * xo[1] - xo[2]
  xoe[nxo + 2] <- 2 * xo[nxo] - xo[nxo - 1]
  xol <- (xoe[2:(nxo + 2)] + xoe[1:(nxo + 1)]) / 2
  xo1 <- xol[1:nxo]
  xo2 <- xol[2:(nxo + 1)]

  xe <- numeric(nx + 2)
  xe[2:(nx + 1)] <- x
  xe[1] <- 2 * x[1] - x[2]
  xe[nx + 2] <- 2 * x[nx] - x[nx - 1]
  xl <- (xe[2:(nx + 2)] + xe[1:(nx + 1)]) / 2
  x1 <- xl[1:nx]
  x2 <- xl[2:(nx + 1)]

  ny <- length(y)
  nyo <- length(yo)

  yoe <- numeric(nyo + 2)
  yoe[2:(nyo + 1)] <- yo
  yoe[1] <- 2 * yo[1] - yo[2]
  yoe[nyo + 2] <- 2 * yo[nyo] - yo[nyo - 1]
  yol <- (yoe[2:(nyo + 2)] + yoe[1:(nyo + 1)]) / 2
  yo1 <- yol[1:nyo]
  yo2 <- yol[2:(nyo + 1)]

  ye <- numeric(ny + 2)
  ye[2:(ny + 1)] <- y
  ye[1] <- 2 * y[1] - y[2]
  ye[ny + 2] <- 2 * y[ny] - y[ny - 1]
  yl <- (ye[2:(ny + 2)] + ye[1:(ny + 1)]) / 2
  y1 <- yl[1:ny]
  y2 <- yl[2:(ny + 1)]

  wx <- matrix(0, nx, nxo)
  for (i in 1:nx) {
    for (j in 1:nxo) {
      d <- min(x2[i], xo2[j]) - max(x1[i], xo1[j])
      wx[i, j] <- d * (d > 0)
    }
  }
  wx <- wx / colSums(wx)[col(wx)]

  wy <- matrix(0, ny, nyo)
  for (i in 1:ny) {
    for (j in 1:nyo) {
      d <- min(y2[i], yo2[j]) - max(y1[i], yo1[j])
      wy[i, j] <- d * (d > 0)
    }
  }
  wy <- wy / colSums(wy)[col(wy)]

  zo <- matrix(0., nxo, nyo)
  for (i in 1:nxo) {
    for (j in 1:nyo) {
      ii <- which(wx[ , i] != 0)
      jj <- which(wy[ , j] != 0)
      ww <- sum(z[ii, jj] * (wx[ii, i] %*% t(wy[jj, j])), na.rm = TRUE)
      zo[i, j] <- (sum(z[ii, jj] * ww, na.rm = TRUE) /
                   sum(is.finite(z[ii, jj]) * ww))
      if (sum(!is.finite(z[ii, jj])) == length(z[ii, jj])) {
        zo[i, j] <- NA
      }
    }
  }
  return(zo)
}
