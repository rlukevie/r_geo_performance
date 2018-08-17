###########################
# 1V02: WRITE VECTOR DATA #
###########################

# =================
# === shapefile ===
# =================
# sp              --> NO (use rgdal)
# sf              --> YES (7.8 sec)
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES (4.3 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# rgdal
library(rgdal)
load("data_input/poly_ssp.RData")
start_time <- Sys.time()
writeOGR(obj = poly_ssp,
         dsn = "./data_output", 
         layer = "poly_ssp",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_ssp)

# sf
library(sf)
load("data_input/poly_ssf.RData")
start_time <- Sys.time()
st_write(poly_ssf, "data_output/poly_ssf.shp")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_ssf)


# =================
# ===  GeoJSON  ===
# =================
# sp              --> NO (use rgdal)
# sf              --> YES, slow (55 sec)
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> YES, but depends on rgdal / sf functions --> not tested
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES, slow (55 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------

# rgdal
library(rgdal)
load("data_input/poly_ssp.RData")
start_time <- Sys.time()
writeOGR(obj = poly_ssp,
         dsn = "./data_output/poly_ssp.geojson", 
         layer = "poly_ssp",
         driver = "GeoJSON",
         overwrite_layer = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_ssp)

# sf
library(sf)
load("data_input/poly_ssf.RData")
start_time <- Sys.time()
st_write(poly_ssf, "data_output/poly_ssf.geojson")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_ssf)

# =================
# ===    KML    ===
# =================
# sp              --> NO (use rgdal)
# sf              --> YES, but CORRUPT FILE (24 sec)
# rgeos           --> NO
# maptools        --> YES, fast, but CORRUPT FILE (6.5 sec)
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES (28 sec)
# RQGIS           --> NO

# fasterize       --> NO
# ----------------------------------


# maptools
library(maptools)
load("data_input/poly_ssp.RData")
start_time <- Sys.time()
kmlPolygons(obj = poly_ssp,
           kmlfile = "data_output/poly_smaptools.kml")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_ssp)

# rgdal
library(rgdal)
load("data_input/poly_ssp.RData")
start_time <- Sys.time()
writeOGR(obj = poly_ssp,
         dsn = "./data_output/poly_ssp.kml", 
         layer = "poly_ssp",
         driver = "KML",
         overwrite_layer = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_ssp)

# sf
library(sf)
load("data_input/poly_ssf.RData")
start_time <- Sys.time()
st_write(poly_ssf, "data_output/poly_ssf.kml")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_ssf)


# ==================
# === SERIALIZED ===
# ==================
# Spatial*        --> 1.3 sec
# sf              --> 0.9 sec

# ----------------------------------
# sp

load("data_input/poly_ssp.RData")
start_time <- Sys.time()
save(poly_ssp, file = "data_output/poly_ssp.RData")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# sf
load("data_input/poly_ssf.RData")
start_time <- Sys.time()
save(poly_ssf, file = "data_output/poly_ssf.RData")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))