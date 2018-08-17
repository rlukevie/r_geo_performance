###########################
# 1R01: READ RASTER DATA  #
###########################

# =================
# ===  GeoTIFF  ===
# =================
# sp              --> NO (use rgdal)
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (0.1 sec)
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES (1.5 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# rgdal
library(rgdal)
start_time <- Sys.time()
raster <- readGDAL(fname = "data_input/raster_m.tif")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# raster
library(raster)
start_time <- Sys.time()
raster <- raster("data_input/raster_m.tif")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# =================
# ===    asc    ===
# =================
# sp              --> YES (27 sec)
# sf              --> NO
# rgeos           --> NO
# maptools        --> YES (26 sec)
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (0.4 sec)
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES (45 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# sp
library(sp)
start_time <- Sys.time()
raster <- read.asciigrid("data_input/raster_m.asc")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# maptools
library(maptools)
start_time <- Sys.time()
raster <- readAsciiGrid("data_input/raster_m.asc")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# rgdal
library(rgdal)
start_time <- Sys.time()
raster <- readGDAL("data_input/raster_m.asc")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# raster
library(raster)
start_time <- Sys.time()
raster <- raster("data_input/raster_m.asc")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))