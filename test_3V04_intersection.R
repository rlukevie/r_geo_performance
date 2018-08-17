######################
# 3V06: INTERSECTION #
######################

# sp              --> NO
# sf              --> YES, with attributes, fast (11 sec)
# rgeos           --> YES, without attributes
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES, with attributes, slow (670 sec), needs a lot of memory
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES, fast (11 sec) but with some errors

# fasterize       --> NO
# ----------------------------------



# sf
library(sf)
poly_s <- st_read("data_input/poly_s.shp")
poly_2_s <- st_read("data_input/poly_2_s.shp")
start_time <- Sys.time()
intersection <- st_intersection(poly_s, poly_2_s)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
st_write(intersection, "data_output/intersection_sf.shp", delete_layer = TRUE)

# rgeos
library(rgdal)
library(rgeos)
poly_s <- readOGR("data_input/poly_s.shp", "poly_s")
poly_2_s <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
start_time <- Sys.time()
intersection <- gIntersection(poly_s, poly_2_s)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
intersection_spdf <- as(intersection, "SpatialPolygonsDataFrame")
writeOGR(intersection_spdf, "./data_output", "intersection_rgeos", driver = "ESRI Shapefile", overwrite_layer = TRUE)
# only geometry as output, no attributes

# raster
library(rgdal)
library(raster)
library(sf)
poly_s <- readOGR("data_input/poly_s.shp", "poly_s")
poly_2_s <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
start_time <- Sys.time()
intersection <- intersect(poly_s, poly_2_s)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
intersection_sf <-st_as_sf(intersection)
st_write(intersection_sf, "data_output/intersection_raster.shp", delete_layer = TRUE)


# RQGIS
library(RQGIS)
library(raster)
library(sf)
library(rgdal)

poly_s <- readOGR("data_input/poly_s.shp", "poly_s_repaired")
poly_2_s <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
# poly_s <- st_read("data_input/poly_s.shp")
# poly_2_s <- st_read("data_input/poly_2_s.shp")

set_env(root = "C:/OSGeo4W64",
        new = TRUE)
# find_algorithms("intersect")
# get_usage("qgis:intersection")
# get_options("qgis:intersection")

start_time <- Sys.time()
params <- get_args_man(alg = "qgis:intersection")
params$INPUT <- poly_2_s
params$INPUT2 <- poly_s
params$OUTPUT <- "data_output/intersection_rqgis.shp"
intersection <- run_qgis(alg = "qgis:intersection",
                                    params = params,
                                    load_output = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
# Result not complete: many polygons are missing! Maybe because of some invalid geometries?