#' Lapse-rate temperature correction downscaling
#'
#' @description Downscales a temperature field using a lapse-rate
#' correction based on a reference orography
#' @author Jost von Hardenberg, \email{j.vonhardenberg@isac.cnr.it}
#' @param lon vector of input longitudes
#' @param lat vector of input latitudes
#' @param t matrix of input temperature data
#' @param lono vector of orography longitudes
#' @param lato vector of orography latitudes
#' @param oro matrix of topographical elevations (in meters) 
#' @param xlim vector of longitude bounds
#' @param ylim vector of latitude bounds
#' @param radius smoothing radius expressed in longitude units (default is half a large-scale pixel)
#' @param lapse environmental lapse rate (in K/Km)
#' @return A downscaled temperature matrix
#' @export
#' @examples
#' lon=5:20
#' lat=35:40
#' t=runif(16*6); dim(t)=c(16,6)
#' lono=seq(5,20,0.1)
#' lato=seq(35,40,0.1)
#' o=runif(151*51)*2000; dim(o)=c(151,51)
#' td=downscalet(lon, lat,t,lono,lato,o,c(8,12),c(36,38))

downscalet <- function(lon, lat, t, lono, lato, oro, xlim, ylim, radius = 0, lapse = 6.5) {

    if (length(dim(t))==3) {
        nt = dim(t)[3]
    } else {
        nt = 1
        dim(t) <- c(dim(t), 1)
    }

    dxl <- lon[2] - lon[1]
    if(radius==0) {
        radius <- dxl/2
    }
    #Add buffer
    lon1 <- xlim[1] - radius
    lon2 <- xlim[2] + radius
    lat1 <- ylim[1] - radius
    lat2 <- ylim[2] + radius

    orocut <- oro[(lono <= lon2) & (lono >= lon1),
                  (lato <= lat2) & (lato >= lat1)]
    lonocut <- lono[(lono <= lon2) & (lono >= lon1)]
    latocut <- lato[(lato <= lat2) & (lato >= lat1)]

    if (length(dim(tcut0))<3) {
        dim(tcut0) <- c(dim(tcut0), 1)
    }

    nx <- length(lonocut)
    ny <- length(latocut)
    
    # Interpolate nn assuming a regular grid 
    tcut <- array(0., c(nx, ny, nt))
    dy <- lat[2]-lat[1]
    dx <- lon[2]-lon[1]
    jj <- ((latocut-(lat[1]-dy/2)) %/% dy + 1)  
    ii <- ((lonocut-(lon[1]-dx/2)) %/% dx + 1)  

    if((ii[1]<0)|(ii[length(ii)]>length(lon))|
       (jj[1]<0)|(jj[length(jj)]>length(lat))) { 
           stop("Downscaling area not contained in input data")
    }
    for (k in seq_len(nt)) {
        for (j in seq_along(jj)) {
            tcut[,j,k] <- t[ii,jj[j],k]
        }
    }

    dxol <- lonocut[2] - lonocut[1]
    nrad <- as.integer(radius / abs(dxol) + 0.5)

    ns = length(lonocut)
    orocuts=smooth(orocut, nrad)
    t1 <- -(orocut-orocuts)*lapse/1000.

    tout <- array(0., c(nx, ny, nt))
    for (i in seq_len(nt)) {
        tout[,, i] <- smooth(tcut[,,i], nrad) + t1
    }
    downscalet <- tout[(lonocut <= xlim[2]) & (lonocut >= xlim[1]),
                       (latocut <= ylim[2]) & (latocut >= ylim[1]),]
}
