##################################
# 2V01: MODIFY VECTOR ATTRIBUTES #
##################################

# sp              --> YES (1.1 sec)
# sf              --> YES (0.9 sec)
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES, but not tested because of good performance of other (standard) libraries

# fasterize       --> NO
# ----------------------------------

# sp
library(sp)
load("data_input/poly_msp.RData")
start_time <- Sys.time()
poly_msp@data$Shape_Area <- '/'(as.numeric(as.character(poly_msp@data$Shape_Area)), as.numeric(as.character(poly_msp@data$Shape_Area)) + 1)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# sf
library(sf)
load("data_input/poly_msf.RData")
start_time <- Sys.time()
poly_msf$Shape_Area <- '/'(as.numeric(as.character(poly_msf$Shape_Area)), as.numeric(as.character(poly_msf$Shape_Area)) + 1)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))


