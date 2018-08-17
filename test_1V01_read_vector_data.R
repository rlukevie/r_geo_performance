##########################
# 1V01: READ VECTOR DATA #
##########################

# =================
# === shapefile ===
# =================
# sp              --> NO (use rgdal)
# sf              --> YES, fast (11 sec)
# rgeos           --> NO
# maptools        --> YES, but deprecated --> not tested
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES, slow (126 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# rgdal
library(rgdal)
start_time <- Sys.time()
poly <- readOGR("data_input/poly_m.shp", "poly_m")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# sf
library(sf)
start_time <- Sys.time()
poly <- st_read("data_input/poly_m.shp")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))


# =================
# ===  GeoJSON  ===
# =================
# sp              --> NO (use rgdal)
# sf              --> YES, slow (83 sec)
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> YES, but depends on rgdal for files > 15 MB --> not tested
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES, very slow (606 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# rgdal
library(rgdal)
start_time <- Sys.time()
poly <- readOGR("data_input/poly_m.geojson", "poly_m")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# sf
library(sf)
start_time <- Sys.time()
poly <- st_read("data_input/poly_m.geojson")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# =================
# ===    KML    ===
# =================
# sp              --> NO (use rgdal)
# sf              --> YES, slow (88 sec)
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> 
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES, very slow (731 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# rgdal
library(rgdal)
start_time <- Sys.time()
poly <- readOGR("data_input/poly_m.kml", "poly_m")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# sf
library(sf)
start_time <- Sys.time()
poly <- st_read("data_input/poly_m.kml")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))


# ==================
# === SERIALIZED ===
# ==================
# Spatial*        --> fast (4 sec)
# sf              --> faster (2.8 sec)

# ----------------------------------
# sp
start_time <- Sys.time()
load("data_input/poly_msp.RData")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# sf
start_time <- Sys.time()
load("data_input/poly_msf.RData")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))