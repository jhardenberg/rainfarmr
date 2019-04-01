#' Aggregate fine field to a coarse array of size `(nas,nas)` by box averaging.
#' @param zi The fine field to aggregate
#' @param nas The target dimension for the output field
#' @return The coarse field obtained by box averaging
#' @export
#' @examples
#' z <- rnorm(64 * 64)
#' dim(z) <- c(64, 64)
#' za <- agg(z, 4)
#' dim(za)
#'   [1] 4 4
agg <- function(zi, nas) {
  ns <- dim(zi)[1]
  sdim <- ns / nas
  rs <- seq(0, ns - sdim, sdim)
  z <- matrix(0., nas, nas);
  n <- matrix(0., nas, nas);
  for (i in 1:sdim) {
    for (j in 1:sdim) {
      imask <- is.finite(zi[i + rs, j + rs])
      n <- n + imask
      z <- z + zi[i + rs, j + rs] * imask
    }
  }
  return(z / n)
}
