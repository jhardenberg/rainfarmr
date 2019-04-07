#' Perform RainFARM downscaling
#'
#' @description The input array is downscaled to finer spatial resolution
#' using the RainFARM stochastic precipitation downscaling method.
#' Orographic correction weights can be applied as described in 
#' Terzago et al. (2018) \doi{10.5194/nhess-18-2825-2018}.
#' Precipitation can be conserved globally (`fglob`), using convolution
#' (`fsmooth`) or over the original coarse-scale boxes. 
#' @references Terzago, S. et al. (2018). NHESS 18(11), 2825–2840
#' \doi{10.5194/nhess-18-2825-2018}; D'Onofrio et al. (2014).
#' J of Hydrometeorology 15, 830-843 \doi{10.1175/JHM-D-13-096.1};
#' Rebora et. al. (2006), JHM 7, 724 \doi{10.1175/JHM517.1}.
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param r matrix or array with large-scale field to downscale. Can be a
#' three-dimensional array with multiple frames at different times.
#' Spatial downscaling is performed separately for each element of the
#' third dimension of `r`.
#' @param slope spatial spectral slope.
#' @param nf refinement factor for spatial downscaling.
#' @param weights matrix with weights for orographic downscaling generated
#' by the [rfweights()] function.
#' @param fglob logical to conserve global average over domain.
#' @param fsmooth logical to use smoothing for conservation. 
#' If neither `fsmooth` or `fglob` is set precipitation is conserved over
#' each coarse pixel of the input field.
#' @param verbose logical to provide some progress report.
#' @return The downscaled array.
#' @export
#' @examples
#' # Make some sample synthetic rainfall data
#' r <- exp(rnorm(4 * 4 * 10))
#' dim(r) <- c(4, 4, 10)
#' r[ , , 1]
#' #           [,1]      [,2]      [,3]      [,4]
#' # [1,] 1.8459816 1.8536550 2.1600665 1.3102116
#' # [2,] 1.3851011 1.4647348 0.2708219 0.4571810
#' # [3,] 0.2492451 0.8227134 0.4790567 1.9320403
#' # [4,] 0.5985922 3.3065678 2.1282795 0.6849944
#' # Downscale with spectral slope=1.7 to size 32x32
#' rd <- rainfarm(r, 1.7, 8, fsmooth=FALSE) 
#' # Verify that downscaled data maintained original box averages
#' agg(rd[ , , 1], 4) 
#' #           [,1]      [,2]      [,3]      [,4]
#' # [1,] 1.8459816 1.8536550 2.1600665 1.3102116
#' # [2,] 1.3851011 1.4647348 0.2708219 0.4571810
#' # [3,] 0.2492451 0.8227134 0.4790567 1.9320403
#' # [4,] 0.5985922 3.3065678 2.1282795 0.6849944
rainfarm <- function(r, slope, nf, weights = 1.,
                     fglob = FALSE, fsmooth = TRUE, verbose = FALSE) {
  drop <- FALSE
  nax <- dim(r)[1]
  nay <- dim(r)[2]
  nt <- dim(r)[3]
  if ( is.na(nt) ) {
    drop <- TRUE
    nt <- 1
  }
  dim(r) <- c(nax, nay, nt)

  if ( nax == nay ) {
    nas <- nax
  } else {
    stop("The spatial dimensions of input array must be square")
  }

  ns <- nas * nf
  f <- initmetagauss(slope, ns)
  rd <- array(0., c(ns, ns, nt))

  for (k in 1:nt) {
    if (verbose) {
       print(paste0("Frame ", k))
    }
    r1 <- r[ , , k]
    if (mean(r1) == 0.) {
        rd[ , , k] <- matrix(0., ns, ns)
    } else {
        fm <- downscale(r1, f, weights, fglob = fglob, fsmooth = fsmooth)
        rd[ , , k] <- fm
    }
  }
  if (drop) {
    rd <- drop(rd)
  }
  return(rd)
}
