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


library(sf)
landuse_2001 <- st_read("data_input/REALNUT2001OGD/REALNUT2001OGDPolygon.shp")
landuse_2001 <- st_read("data_input/REALNUT2016GOGD/REALNUT2016GOGDPolygon.shp")
start_time <- Sys.time()
st_union()