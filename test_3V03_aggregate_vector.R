#################################
# 3V03: AGGREGATE VECTOR DATA   #
#################################

# sp              --> YES (0.6 sec)
# sf              --> YES (0.8 sec)
# rgeos           --> NO
# maptools        --> YES, but not tested because of good performance of other (standard) libraries
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (0.7 sec)
# spatial.tools   --> NO
# gdalUtils       --> YES, via ogr2ogr and sqlite (2 sec)
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES, but not tested because of good performance of other (standard) libraries

# fasterize       --> NO
# ----------------------------------


# gdalUtils
library(gdalUtils)
library(raster)

start_time <- Sys.time()
ogr2ogr(src_datasource_name = "data_input/poly_2_s.shp",
        dst_datasource_name = "data_output/poly_agg_gdalutils.shp",
        dialect = "sqlite",
        sql = "SELECT ST_Union(geometry), SUM(area) FROM poly_2_s GROUP BY BEZIRK")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))


# sp
library(rgdal)
library(sp)

poly <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
poly@data <- poly@data[, c("BEZIRK", "area")]
start_time <- Sys.time()
poly@data$area <- as.numeric(as.character(poly@data$area))
poly_agg <- aggregate(poly, by = list(bezirk = poly$BEZIRK), FUN = sum)
writeOGR(poly_agg, "./data_output", "poly_agg_sp", driver = "ESRI Shapefile", overwrite_layer = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# raster
library(rgdal)
library(raster)
poly <- readOGR("data_input/poly_2_s.shp", "poly_2_s")

start_time <- Sys.time()

poly@data$area <- as.numeric(as.character(poly@data$area))
poly_agg <- aggregate(poly, by = "BEZIRK", sums = list(list(sum, c("area"))))
writeOGR(poly_agg, "./data_output", "poly_agg_raster", driver = "ESRI Shapefile", overwrite_layer = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))



# sf
library(sf)
poly <- st_read("data_input/poly_2_s.shp")

start_time <- Sys.time()
poly <- poly[, c("BEZIRK", "area", "geometry")]
poly_agg <- aggregate(poly, by = list(poly$BEZIRK), FUN = sum)
st_write(poly_agg, "data_output/poly_agg_sf.shp", delete_dsn = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

l