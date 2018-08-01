#############
# MASKIEREN #
#############

# sp --> NEIN
# sf --> NEIN
# rgeos --> NEIN
# maptools --> NEIN
# geojsonio --> NEIN
# rmapshaper --> NEIN

# raster --> JA, langsam bei großer Maske
#        --> mit vorhergehendem fasterize deutlich besser
# spatial.tools --> NEIN
# gdalUtils --> eventuell mit gdal_translate
# landscapetools
# velox

# rgdal --> NEIN
# RQGIS --> JA, aber Zellstatistik stimmt nicht

# fasterize --> NEIN (nur als Zwischenschritt mit raster)


# ----------------------------------


# raster mit vorhergehendem crop
library(raster)
library(sf)
start_time <- Sys.time()
raster <- raster("data_input/raster_m.tif")
mask_sp <- as(st_read("data_input/crop_m.shp"), "Spatial")
raster <- crop(raster, mask_sp)
crs(mask_sp) <- crs(raster)
raster_mask_raster <- mask(raster, mask_sp)
print(Sys.time() - start_time)
cellStats(raster_mask_raster, stat = 'sum')
plot(raster_mask_raster)

# raster mit fasterize
library(raster)
library(sf)
start_time <- Sys.time()
raster <- raster("data_input/raster_m.tif")
mask_sf <- st_read("data_input/crop_m.shp")
raster <- crop(raster, extent(as(mask_sf, "Spatial")))
mask_raster <- fasterize::fasterize(sf = mask_sf, raster = raster)
crs(mask_raster) <- crs(raster)
raster_mask_raster_fasterize <- mask(raster, mask_raster)
print(Sys.time() - start_time)
cellStats(raster_mask_raster_fasterize, stat = 'sum')
plot(raster_mask_raster_fasterize)

# RQGIS
library(raster)
library(RQGIS)
library(rgdal)
start_time <- Sys.time()
raster <- raster("data_input/raster_s.tif")
mask_sp <- readOGR("data_input/crop_s.shp", "crop_s")
crs(mask_sp) <- crs(raster)
set_env(root = "C:/OSGeo4W64",
        new = TRUE)
params <- get_args_man(alg = "gdalogr:cliprasterbymasklayer")
print(params)
params$INPUT <- raster
params$MASK <- mask_sp
params$CROP_TO_CUTLINE <- TRUE
params$KEEP_RESOLUTION <- TRUE
# params$RTYPE <- 0
params$OUTPUT <- "data_output/raster_mask_RQGIS.tif"


raster_mask_rqgis <- run_qgis(alg = "gdalogr:cliprasterbymasklayer",
                                    params = params,
                                    load_output = TRUE)
print(Sys.time() - start_time)
cellStats(raster_mask_rqgis, stat = 'sum')
plot(raster_mask_rqgis)
