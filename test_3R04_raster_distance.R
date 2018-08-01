#########################
# 3R04: RASTER DISTANCE #
#########################

# sp              --> NO
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES, very slow
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES, fast (26 sec)

# fasterize       --> NO
# ----------------------------------

# raster
library(raster)
raster_streets <- raster("data_input/raster_streets.tif")
start_time <- Sys.time()
raster_streets_distance_raster <- distance(raster_streets)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
plot(raster_streets_distance_raster)

# RQGIS
library(RQGIS)
library(raster)

raster_streets <- raster("data_input/raster_streets.tif")

set_env(root = "C:/OSGeo4W64",
        new = TRUE)

# find_algorithms("proximity")
# get_usage("gdalogr:proximity")
# get_options("gdalogr:proximity")

start_time <- Sys.time()
params <- get_args_man(alg = "gdalogr:proximity")
print(params)
params$INPUT <- raster_streets
params$OUTPUT <- "data_output/raster_streets_RQGIS.tif"
params$VALUES <- 30201
params$RTYPE <- 5
params$UNITS <- 0

raster_streets_distance <- run_qgis(alg = "gdalogr:proximity",
         params = params,
         load_output = TRUE)

print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
plot(raster_streets_distance)
