#' Interpolate field `z` to size `(ns,ns)` using nearest neighbors.
#' @param z The input field at coarse resolution
#' @param ns The target size
#' @return The fine-scale field with dimensions `c(ns, ns)`
#' @export
#' @examples
#' za <- rnorm(4 * 4)
#' dim(za) <- c(4, 4)
#' z <- interpola(z, 16)
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
