#################################
# 4RV04: RASTER EXTRACT POLYGON #
#################################

# sp              --> 
# sf              --> NO
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> YES, slow (5 min)
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> YES, fast (13 sec sec)

# rgdal           --> NO
# RQGIS           --> NO (QGIS has no extract function for polygons)

# fasterize       --> NO
# ----------------------------------


# raster
# Example: extract landuse in Vienna

library(sf)
library(rgdal)
library(raster)
library(fasterize)
library(velox)

start_time <- Sys.time()
sprengel <- readOGR("data_input/poly_2_s.shp", "poly_2_s")

RESOLUTION = 10

# landuse <- st_read("data_input/REALNUT2016GOGD/REALNUT2016GOGDPolygon.shp")
landuse <- st_read("data_input/REALNUT200708OGD/REALNUT200708OGDPolygon.shp")
landuse <- st_transform(landuse, 31256)
landuse_raster <- raster(landuse, res = RESOLUTION, val = 1)
landuse_raster <- fasterize(sf = landuse, raster = landuse_raster, field = "NUTZUNG_CO")

ext <- extract(landuse_raster, sprengel)
res = res(landuse_raster)
data_sprengel <- data.frame()

for (i in 1:length(ext)) {
  tab <- table(ext[[i]])
  unique_values_one_sprengel <- attr(tab, "names")
  counts_one_sprengel <- as.vector(tab)
  nrows <- nrow(tab)
  for (u in 1:nrows) {
    value_string <- toString(unique_values_one_sprengel[[u]])
    count <- counts_one_sprengel[[u]]
    data_sprengel[i, value_string] <- count
  }
}

sprengel@data <- cbind(sprengel@data, data_sprengel)
writeOGR(sprengel, "./data_output", "sprengel_nutzung_count_200708_raster", driver = "ESRI Shapefile", overwrite_layer = TRUE)

print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

Sys.time()- start_time

# velox
# Example: extract landuse in Vienna

library(sf)
library(rgdal)
library(raster)
library(fasterize)
library(velox)


start_time <- Sys.time()
sprengel <- readOGR("data_input/poly_2_s.shp", "poly_2_s")

res = 10
# landuse <- st_read("data_input/REALNUT2016GOGD/REALNUT2016GOGDPolygon.shp")
landuse <- st_read("data_input/REALNUT200708OGD/REALNUT200708OGDPolygon.shp")
landuse <- st_transform(landuse, 31256)
landuse_raster <- raster(landuse, res = res, val = 1)
landuse_raster <- fasterize(sf = landuse, raster = landuse_raster, field = "NUTZUNG_CO")

landuse_raster_vx <- velox(landuse_raster)
ext <- landuse_raster_vx$extract(sprengel)
data_sprengel <- data.frame()
res_factor = res(landuse_raster)[1]*res(landuse_raster)[2]

for (i in 1:length(ext)) {
  tab <- table(ext[[i]])
  unique_values_one_sprengel <- attr(tab, "names")
  counts_one_sprengel <- as.vector(tab)
  nrows <- nrow(tab)
  for (u in 1:nrows) {
    value_string <- toString(unique_values_one_sprengel[[u]])
    count <- counts_one_sprengel[[u]]
    data_sprengel[i, value_string] <- count * res_factor
  }
}


sprengel@data <- cbind(sprengel@data, data_sprengel)
writeOGR(sprengel, "./data_output", "sprengel_nutzung_count_200708_velox", driver = "ESRI Shapefile", overwrite_layer = TRUE)

print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
      