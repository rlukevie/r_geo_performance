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
landuse_2016 <- st_read("data_input/REALNUT2016GOGD/REALNUT2016GOGDPolygon.shp")
landuse_2016 <- st_transform(landuse_2016, 31256)
landuse_2016_raster <- raster(landuse_2016, res = res, val = 1)
landuse_2016_raster <- fasterize(sf = landuse_2016, raster = landuse_2016_raster, field = "NUTZUNG_CO")

landuse_2016_raster_vx <- velox(landuse_2016_raster)
ext <- landuse_2016_raster_vx$extract(sprengel)
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
writeOGR(sprengel, "./data_ouput", "sprengel_nutzung_count_2016", driver = "ESRI Shapefile", overwrite_layer = TRUE)

print(paste0("Duration: ", Sys.time() - start_time))

