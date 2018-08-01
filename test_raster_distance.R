# sp --> NEIN
# sf --> NEIN
# rgeos --> NEIN
# maptools --> NEIN
# geojsonio --> NEIN
# rmapshaper --> NEIN

# raster --> JA, aber sehr langsam
# spatial.tools --> NEIN
# gdalUtils --> NEIN
# landscapetools --> NEIN
# velox --> NEIN

# rgdal --> NEIN
# RQGIS --> JA, schnell

# fasterize --> NEIN
# ----------------------------------

# raster
library(raster)
raster_streets_s <- raster("data_input/raster_streets_s.tif")
raster_streets_distance_raster <- distance(raster_streets_s)
plot(raster_streets_distance_raster)

# RQGIS
library(RQGIS)
library(raster)

raster_streets <- raster("data_input/raster_streets.tif")

set_env(root = "C:/OSGeo4W64",
        new = TRUE)

find_algorithms("proximity")

get_usage("gdalogr:proximity")

get_options("gdalogr:proximity")

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

plot(raster_streets_distance)
