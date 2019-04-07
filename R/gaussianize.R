#' Gaussianize field using rank ordering
#'
#' @description The amplitudes of the input field are rank ordered and
#' substituted with the corresponding rank-ordered samples from a
#' Normal distribution.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param z matrix containing the input field to be Gaussianized.
#' @return The Gaussianized field with the same dimensions as the input field.
#' @export
#' @examples
#' # Make some sample synthetic rainfall data
#' r <- exp(rnorm(64 * 64))
#' dim(r) <- c(64, 64)
#' # Gaussianize and check standard deviation
#' g <- gaussianize(r)
#' sd(g)
#' # [1] 1
gaussianize <- function(z) {
  m <- dim(z)[1]
  n <- dim(z)[2]
  nn <- m * n
  ii <- order(z)
  g <- matrix(0, m, n)
  g[ii] <- sort(rnorm(nn))
  gs <- sd(g)
  if (gs == 0) {
    gs <- 1
  }
  g <- g / gs
  return(g)
}
