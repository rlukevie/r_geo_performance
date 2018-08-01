#################
# UNION         #
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


# sf
library(sf)

landuse_200708 <- st_read("data_input/REALNUT200708OGD/REALNUT200708OGDPolygon.shp")
landuse_2016 <- st_read("data_input/REALNUT2016GOGD/REALNUT2016GOGDPolygon.shp")
poly_s <- st_read("data_input/poly_s.shp")
poly_2_s <- st_read("data_input/poly_2_s.shp")
point_s <- st_read("data_input/point_s.shp")
start_time <- Sys.time()
union <- st_union(poly_s, poly_2_s, by_feature = FALSE)

print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))


# rgeos
library(rgdal)
library(rgeos)
poly_s <- readOGR("data_input/poly_s.shp", "poly_s")
poly_2_s <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
union <- gUnion(poly_s, poly_2_s)
plot(union)

# raster
library(rgdal)
library(raster)
poly_s <- readOGR("data_input/poly_s.shp", "poly_s")
poly_2_s <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
union <- raster::union(poly_s, poly_2_s)

