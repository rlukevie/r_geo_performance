############################
# 1R01: WRITE RASTER DATA  #
############################

# =================
# ===  GeoTIFF  ===
# =================
# sp              --> NO (use rgdal)
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (1.8 sec)
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES (0.2 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# rgdal
library(rgdal)
raster <- readGDAL("data_input/raster_s.tif")
start_time <- Sys.time()
writeGDAL(dataset = raster,
         fname = "data_output/raster_srgdal.tif")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(raster)

# raster
library(raster)
raster <- raster("data_input/raster_s.tif")
start_time <- Sys.time()
writeRaster(raster, "data_output/raster_sraster.tif")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(raster)

# =================
# ===    asc    ===
# =================
# sp              --> YES (5.5 sec)
# sf              --> NO
# rgeos           --> NO
# maptools        --> Yes (5.5 sec)
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (8 sec)
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES (0.2 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# rgdal
library(rgdal)
raster <- readGDAL("data_input/raster_s.tif")
start_time <- Sys.time()
writeGDAL(dataset = raster,
          fname = "data_output/raster_srgdal.asc")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(raster)

# raster
library(raster)
raster <- raster("data_input/raster_s.tif")
start_time <- Sys.time()
writeRaster(raster, "data_output/raster_sraster.asc")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(raster)

# maptools
library(maptools)
raster <- readAsciiGrid("data_input/raster_s.asc")
start_time <- Sys.time()
writeAsciiGrid(raster, "data_output/raster_smaptools.asc")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(raster)

# sp
library(sp)
raster <- read.asciigrid("data_input/raster_s.asc")
start_time <- Sys.time()
write.asciigrid(x = raster, fname = "data_output/raster_ssp.asc")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(raster)