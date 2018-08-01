# sp
# sf
# rgeos
# maptools
# geojsonio
# rmapshaper

# raster --> JA
# spatial.tools
# gdalUtils
# landscapetools
# velox

# rgdal
# RQGIS

# fasterize
# ----------------------------------

# raster
library(raster)
library(sf)
raster <- raster("data_input/raster_l.tif")
crop <- st_read("data_input/crop.shp")
ext <- extent(as(crop, "Spatial"))
raster_crop_raster <- crop(raster, ext)
plot(raster_crop_raster)

# velox
library(raster)
library(velox)
raster <- raster("data_input/raster_l.tif")
crop <- st_read("data_input/crop.shp")
raster_vx <- velox(raster)
raster_vx$crop(crop)
ext <- extent(as(crop, "Spatial"))
crs <- st_crs(crop)
raster_crop_velox <- raster(x = raster_vx$rasterbands[[1]],
                              xmn = ext@xmin,
                              xmx = ext@xmax,
                              ymn = ext@ymin,
                              ymx = ext@ymax,
                              crs = crs)
plot(raster_crop_velox)
