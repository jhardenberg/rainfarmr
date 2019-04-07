#' Interpolate field using nearest neighbors
#'
#' @description Interpolates a square input matrix to a finer
#' resolution `ns` using nearest neighbours.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param z matrix containing the input field at coarse resolution.
#' @param ns the target size.
#' @return The resulting fine-scale field with dimensions `c(ns, ns)`.
#' @export
#' @examples
#' za <- rnorm(4 * 4)
#' dim(za) <- c(4, 4)
#' z <- interpola(za, 16)
#' dim(z)
#' # [1] 16 16
interpola <- function(z, ns) {
  nas <- dim(z)[1]
  sdim <- ns / nas
  rs <- seq(0, ns - sdim, sdim)
  zi <- matrix(0., ns, ns)
  for (i in 1:sdim) {
    for (j in 1:sdim) {
      zi[i + rs, j + rs] <- z
    }
  }
  return(zi)
}
