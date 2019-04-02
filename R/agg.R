#' Aggregate field using box-averaging
#' @description Aggregate fine input field `z` to a coarse array sith dimensions
#' `c(nas,nas)` by box averaging.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param z matrix containing the fine field to aggregate.
#' @param nas target dimension for the output field.
#' @return The coarse field obtained by box averaging.
#' @export
#' @examples
#' z <- rnorm(64 * 64)
#' dim(z) <- c(64, 64)
#' za <- agg(z, 4)
#' dim(za)
#' # [1] 4 4
agg <- function(z, nas) {
  ns <- dim(z)[1]
  sdim <- ns / nas
  rs <- seq(0, ns - sdim, sdim)
  za <- matrix(0., nas, nas);
  n <- matrix(0., nas, nas);
  for (i in 1:sdim) {
    for (j in 1:sdim) {
      zz <- z[i + rs, j + rs]
      imask <- is.finite(zz)
      n <- n + imask
      zz[!imask] <- 0.
      za <- za + zz
    }
  }
  return(za / n)
}
