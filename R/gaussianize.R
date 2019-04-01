#' Gaussianize field using a sorted random sequence.
#' @param z The field to be sorted
#' @return The gaussianized field with the same dimensions as the input field
#' @export
#' @examples
#' z <- rnorm(64 * 64)
#' dim(z) <- c(64, 64)
#' g <- gaussianize(z)
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
