##################
# RASTER EXTRACT #
##################

# sp --> 
# sf --> 
# rgeos --> 
# maptools --> 
# geojsonio --> 
# rmapshaper --> 

# raster --> JA (Performance testen)
# spatial.tools --> 
# gdalUtils --> 
# landscapetools
# velox --> JA, flott

# rgdal --> 
# RQGIS

# fasterize --> 


# ----------------------------------


# velox
# Anwendungsbeispiel Extraktion Realnutzungskartierung

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
writeOGR(sprengel, "./data_ouput", "sprengel_nutzung_count_200708", driver = "ESRI Shapefile", overwrite_layer = TRUE)

print(paste0("Duration: ", Sys.time() - start_time))

