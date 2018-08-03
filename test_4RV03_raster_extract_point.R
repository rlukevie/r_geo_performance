#################################
# 4RV04: RASTER EXTRACT POINT   #
#################################

# sp              --> YES (0.12 sec)
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (0.6 sec)
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES, but only via SAGA or GRASS (not tested)

# fasterize       --> NO
# ----------------------------------
st_over


# sp
library(sp)
library(rgdal)

points <- readOGR("data_input/point_l.shp", "point_l")
grid <- readGDAL("data_input/raster_landuse200708.tif")
crs(grid) <- crs(points)

start_time <- Sys.time()
extract <- over(y = grid, x = points)

points@data <- cbind(points@data, extract)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

writeOGR(points, "./data_output", "raster_extract_point_sp", driver = "ESRI Shapefile", overwrite_layer = TRUE)

# raster
library(rgdal)
library(raster)

points <- readOGR("data_input/point_l.shp", "point_l")
raster <- raster("data_input/raster_landuse200708.tif")

start_time <- Sys.time()
landuse_p <- extract(raster, points)

points@data <- cbind(points@data, landuse_p)

print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
writeOGR(points, "./data_output", "raster_extract_point_raster", driver = "ESRI Shapefile", overwrite_layer = TRUE)


