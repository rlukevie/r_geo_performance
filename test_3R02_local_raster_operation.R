#################################
# 3R02: LOCAL RASTER OPERATION  #
#################################


# sp              --> NO
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES (14.2 sec)
# spatial.tools   --> YES (32.9 sec with 1 CPU core, 23.5 sec with 3 CPU cores,  21.3 sec with 4 CPU cores)
# gdalUtils       --> NO (gdal_calc is missing)
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES (15.23 sec)

# fasterize       --> NO
# ----------------------------------


# raster
library(raster)
rasterOptions(maxmemory = 10000)
raster_input <- raster("data_input/raster_m.tif")

fun <- function(x) { # nonsensical test function
  ifelse(x < 80, x / (x + 1) * (x + 2), x / (x + 1) * (x + 2) + 1000)
}

start_time <- Sys.time()
raster_local <- calc(x = raster_input, fun = fun)
writeRaster(raster_local, "data_output/raster_local_raster.tif", overwrite = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

plot(raster_local)

# spatial.tools
library(spatial.tools)

fun <- function(inraster) { # nonsensical test function
  ifelse(inraster < 80, inraster / (inraster + 1) * (inraster + 2), 
         inraster / (inraster + 1) * (inraster + 2) + 1000)
}
raster <- setMinMax(raster("data_input/raster_m.tif"))

start_time <- Sys.time()
sfQuickInit(cpus=4)
raster_local_spatial_tools <- rasterEngine(
  inraster = raster,
  fun = fun)
sfQuickStop()

writeRaster(raster_local_spatial_tools, "data_output/raster_local_spatial_tools.tif", overwrite = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

plot(raster_local_spatial_tools)


# RQGIS
library(RQGIS)
library(raster)
library(sf)
library(rgdal)

raster_input <- raster("data_input/raster_m.tif")

set_env(root = "C:/OSGeo4W64",
        new = TRUE)
find_algorithms("calc")
get_usage("gdalogr:rastercalculator")
get_options("gdalogr:rastercalculator")

start_time <- Sys.time()
params <- get_args_man(alg = "gdalogr:rastercalculator")
params
params$INPUT_A <- raster_input
params$FORMULA <- "where(A<80,A/(A+1)*(A+2),A/(A+1)*(A+2)+1000)"
params$OUTPUT <- "data_output/raster_local_qgis.tif"
params$RTYPE <- 5
raster_local_qgis_gdal <- run_qgis(alg = "gdalogr:rastercalculator",
                         params = params,
                         load_output = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
