#################################
# 3R03: FOCAL RASTER OPERATION  #
#################################


# sp              --> NO
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (257 sec)
# spatial.tools   --> YES (155 sec with 3 CPU cores)
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> YES, but no custom functions (only mean, median, sum) --> not tested

# rgdal           --> NO
# RQGIS           --> YES, e.g. via SAGA or GRASS, but no custom functions (only custom matrices) --> not tested

# fasterize       --> NO
# ----------------------------------



# raster
library(raster)
raster_input <- raster("data_input/raster_s.tif")


fun <- function(x) {
  quantiles <- quantile(x, names = FALSE, na.rm = TRUE)
  return(quantiles[2] - 1.5*(quantiles[4] - quantiles[2]))
}


start_time <- Sys.time()
window <- matrix(data = 1, nrow = 21, ncol = 21)
raster_focal_raster <- focal(x = raster_input, w = window, fun = fun)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

plot(raster_focal_raster)



# spatial.tools
library(spatial.tools)
library(raster)

raster_input <- raster("data_input/raster_s.tif")

fun <- function(inraster) {
  quantiles <- quantile(inraster, names = FALSE, na.rm = TRUE)
  return(quantiles[2] - 1.5*(quantiles[4] - quantiles[2]))
}

start_time <- Sys.time()
sfQuickInit(cpus=3)
raster_focal_spatial_tools <- rasterEngine(
  inraster = raster_input,
  fun = fun,
  window_dims = c(21,21))
sfQuickStop()
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

plot(raster_focal_spatial_tools)
