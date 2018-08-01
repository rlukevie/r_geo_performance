####################
# 4RV01: RASTERIZE #
####################

# sp              --> NO
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES, slow (49 sec)
# spatial.tools   --> NO
# gdalUtils       --> YES, fast (2 sec, including write GeoTIFF)
# landscapetools  --> NO
# velox           --> YES, medium performance (8 sec)

# rgdal           --> NO
# RQGIS           --> YES, slow with grass7 algorithm (30 sec)

# fasterize       --> YES, fast (0.55 sec)
# ----------------------------------

# rqgis (gdalogr:rasterize)
library(RQGIS)
library(raster)
poly_s <- readOGR("data_input/poly_s.shp", "poly_s")
set_env(root = "C:/OSGeo4W64",
        new = TRUE)

find_algorithms("rasterize")
get_usage("grass7:v.to.rast.attribute")
get_options("grass7:v.to.rast.attribute")

start_time <- Sys.time()
params <- get_args_man(alg = "grass7:v.to.rast.attribute")
print(params)
params$input <- poly_s
params$attribute_column = "F_KLASSE"
params$GRASS_REGION_CELLSIZE_PARAMETER = 10
params$output <- "data_output/rasterize_rqgis_grass7.tif"


poly_s_raster <- run_qgis(alg = "grass7:v.to.rast.attribute",
                                    params = params,
                                    load_output = TRUE)

print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
plot(poly_s_raster)


# gdalutils
library(gdalUtils)
library(raster)
start_time <- Sys.time()
poly_s_raster <- gdal_rasterize(src_datasource = "data_input/poly_s.shp",
               dst_filename = "data_output/rasterize_gdalutils.tiff",
               a = "F_KLASSE",
               ot = "Byte",
               tr = "10 10",
               l = "poly_s",
               output_Raster = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
plot(poly_s_raster)


# fasterize
library(fasterize)
library(sf)
library(raster)
res = 10
poly_s <- st_read("data_input/poly_s.shp")
start_time <- Sys.time()
raster <- raster(poly_s, res = res, val = 99999999)
poly_s_raster_fasterize <- fasterize::fasterize(sf = poly_s, raster = raster, field = "F_KLASSE", fun = "first")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
plot(poly_s_raster_fasterize)
writeRaster(poly_s_raster_fasterize, "data_output/rasterize_fasterize.tiff", overwrite = TRUE)

# raster::rasterize
library(rgdal)
library(raster)
res = 10
poly_ssp <- readOGR("data_input/poly_s.shp", "poly_s")
start_time <- Sys.time()
raster <- raster(poly_ssp, res = res, val = 99999999)
poly_s_raster_rasterize <- raster::rasterize(x = poly_ssp, y = raster, field = "F_KLASSE", fun = "first")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
plot(poly_s_raster_rasterize)
writeRaster(poly_s_raster_rasterize, "data_output/rasterize_raster.tiff", overwrite = TRUE)

# velox
# ?VeloxRaster_rasterize
library(sf)
library(raster)
library(velox)
res = 10
poly_ssp <- readOGR("data_input/poly_s.shp", "poly_s")
start_time <- Sys.time()
raster <- raster(poly_ssp, res = res, val = 99999999)
raster_vx <- velox(raster)
raster_vx$rasterize(poly_ssp, field = "F_KLASSE", band = 1, small = FALSE)
extent <- extent(raster)
crs <- crs(raster)
poly_s_raster_velox <- raster(x = raster_vx$rasterbands[[1]],
                              xmn = extent@xmin,
                              xmx = extent@xmax,
                              ymn = extent@ymin,
                              ymx = extent@ymax,
                              crs = crs)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
plot(poly_s_raster_velox)
writeRaster(poly_s_raster_velox, "data_output/rasterize_velox.tiff", NAflag = 99999999, overwrite = TRUE)
