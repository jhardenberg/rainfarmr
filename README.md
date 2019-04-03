# rainfarmr
![RainFARM logo](https://raw.githubusercontent.com/jhardenberg/RainFARM.jl/master/docs/src/assets/logo-small.png)

Stochastic precipitation downscaling with the RainFARM method.

`rainfarmr` is a R package implementing the RainFARM (Rainfall Filtered Autoregressive Model) stochastic precipitation downscaling method. Adapted for climate downscaling according to (D'Onofrio et al. 2018) and with fine-scale orographic weights (Terzago et al. 2018).

RainFARM (Rebora et al. 2006) is a metagaussian stochastic downscaling procedure based on the extrapolation of the coarse-scale Fourier power spectrum of a spatio-temporal precipitation field to small scales.

## Example

```r
# Make some sample synthetic rainfall data
# 10 fields of 8 by 8 pixel resolution
r <- exp(rnorm(8 * 8 * 10))
# The corresponding latitudes and longitudes
lon <- seq(5, 8.5, 0.5)
lat <- seq(43.5, 47, 0.5)
dim(r) <- c(8, 8, 10)
nf <- 8  # This is the factor by which we will increase resolution
# Downscale with spectral slope=1.7 to size 64x64
rd <- rainfarm(r, 1.7, nf, fsmooth = TRUE) 
# Get the corresponding fine-scale longitudes and latititudes
grid <- lon_lat_fine(lon, lat, nf)
grid$lon[1:4]
# [1] 4.78125 4.84375 4.90625 4.96875
```

## References

- Terzago, S., Palazzi, E., and von Hardenberg, J. (2018). Stochastic downscaling of precipitation in complex orography: a simple method to reproduce a realistic fine-scale climatology, Nat. Hazards Earth Syst. Sci., 18, 2825-2840, doi: <https://doi.org/10.5194/nhess-18-2825-2018>

- D’Onofrio, D., Palazzi, E., von Hardenberg, J., Provenzale, a., & Calmanti, S. (2014). Stochastic Rainfall Downscaling of Climate Models. Journal of Hydrometeorology, 15(2), 830–843. doi: <https://doi.org/10.1175/JHM-D-13-096.1>

- Rebora, N., Ferraris, L., von Hardenberg, J., & Provenzale, A. (2006). RainFARM: Rainfall Downscaling by a Filtered Autoregressive Model. Journal of Hydrometeorology, 7(4), 724–738. doi: <https://doi.org/10.1175/JHM517.1>

## Other languages

A [julia version](https://github.com/jhardenberg/RainFARM.jl) by the same author is available.

## Authors: 

*R version* - J. von Hardenberg, ISAC-CNR (2019)
