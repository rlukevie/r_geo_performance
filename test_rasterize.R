#################
# RASTERISIEREN #
#################

# sp --> 
# sf --> 
# rgeos --> 
# maptools --> 
# geojsonio --> 
# rmapshaper --> 

# raster --> 
# spatial.tools --> 
# gdalUtils --> 
# landscapetools
# velox

# rgdal --> 
# RQGIS

# fasterize --> 


# ----------------------------------


library(fasterize)
library(sf)
library(raster)
library(velox)

res = 10

# fasterize
poly_s <- st_read("data_input/poly_s.shp")
raster <- raster(poly_s, res = res, val = 99999999)
poly_s_raster_fasterize <- fasterize::fasterize(sf = poly_s, raster = raster, field = "F_KLASSE", fun = "first")
plot(poly_s_raster_fasterize)
writeRaster(poly_s_raster_fasterize, "data_output/rasterize_fasterize.tiff", overwrite = TRUE)

# raster::rasterize
poly_ssp <- as(poly_s, "Spatial")
poly_s_raster_rasterize <- raster::rasterize(x = poly_ssp, y = raster, field = "F_KLASSE", fun = "first")
plot(poly_s_raster_rasterize)
writeRaster(poly_s_raster_rasterize, "data_output/rasterize_raster.tiff", overwrite = TRUE)
# Ergebnis: nicht so fein aufgelöst wie fasterize: bei Polygonen, die gemeinsam geschlossene Flächen bilden --> ausgefüllt

# velox
# ?VeloxRaster_rasterize
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
plot(poly_s_raster_velox)

writeRaster(poly_s_raster_velox, "data_output/rasterize_velox.tiff", NAflag = 99999999, overwrite = TRUE)
