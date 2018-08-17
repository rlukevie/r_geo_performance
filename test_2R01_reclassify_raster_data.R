################################
# 2R01: RECLASSIFY RASTER DATA #
################################

# sp              --> NO
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (4.4 sec)
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES (GRASS7: 3.3 sec)

# fasterize       --> NO
# ----------------------------------

# raster
library(raster)
raster = raster("data_input/raster_m.tif")
start_time <- Sys.time()
rcl = matrix(c(-200, 0, 1, 0, 100, 2, 100, 300, 3), ncol = 3, byrow = TRUE)
raster <- reclassify(raster, rcl = rcl)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# RQGIS
library(RQGIS)
library(raster)
raster = raster("data_input/raster_m.tif")
set_env(root = "C:/OSGeo4W64",
        new = TRUE)
# find_algorithms("recl")
# get_usage("grass7:r.reclass")
# get_options("grass7:r.reclass")
start_time <- Sys.time()
params <- get_args_man(alg = "grass7:r.reclass")
params$input <- "data_input/raster_m.tif"
params$txtrules <- " -200 thru 0 = 1 0.1 thru 100 = 2 100.1 thru 300 = 3 end"
params$output <- "data_output/raster_rec_qgis_grass7.tif"
raster_rec <- run_qgis(alg = "grass7:r.reclass",
                         params = params,
                         load_output = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
