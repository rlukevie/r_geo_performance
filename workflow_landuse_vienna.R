# velox
library(sf)
library(rgdal)
library(raster)
library(fasterize)
library(velox)
library(sp)
library(tibble)
library(dplyr)
library(readxl)

start_time <- Sys.time()

RESOLUTION = 5

# load input data
landuse2001 <- st_read("data_input/REALNUT2001OGD/REALNUT2001OGDPolygon_31256.shp")
# landuse2001 <- st_transform(landuse2001, 31256)
# st_write(landuse2001, "data_input/REALNUT2001OGD/REALNUT2001OGDPolygon_31256.shp")

# landuse2016 <- st_read("data_input/REALNUT2016GOGD/REALNUT2016GOGDPolygon.shp")
# landuse2016 <- st_transform(landuse2016, 31256)
# st_write(landuse2016, "data_input/REALNUT2016GOGD/REALNUT2016GOGDPolygon_31256.shp", delete_dsn = TRUE)

landuse2016 <- readOGR("data_input/REALNUT2016GOGD/REALNUT2016GOGDPolygon_31256.shp", "REALNUT2016GOGDPolygon_31256")

# prepare transformation of nominal landuse classes into integer values (for year 2001)
# this is needed for rasterization
classes_nominal <- sort(unique(na.omit(landuse2001$NUTCD)))
landuse2001$NUTCD_INT <- match(landuse2001$NUTCD, classes_nominal)


# rasterize landuse2001
landuse2001_raster <- raster(landuse2001, res = RESOLUTION, val = 1)
landuse2001_raster <- fasterize(sf = landuse2001, raster = landuse2001_raster, field = "NUTCD_INT")

# convert to velox object and extract into landuse2016
landuse2001_raster_vx <- velox(landuse2001_raster)
ext <- landuse2001_raster_vx$extract(landuse2016)


# write extraction result into Spatial object landuse2016 
start_time_loop <- Sys.time()
data_landuse2016 <- data.frame()
for (i in 1:length(ext)) {
  counts_sum <- 0
  tab <- table(ext[[i]])
  unique_values_one_2016unit <- attr(tab, "names")
  counts_one_2016unit <- as.vector(tab)
  nrows <- nrow(tab)
  if (nrows > 0) {
    for (u in 1:nrows) {
      value = as.numeric(unique_values_one_2016unit[[u]])
      value_string <- toString(classes_nominal[[value]])
      count <- counts_one_2016unit[[u]]
      area <- count * RESOLUTION * RESOLUTION
      # print(paste0("i: ", i, " u: ", u, " value: ", value_string, " area: ", area, " len: ", length(ext[[i]]), " counts: ", counts_one_2016unit[[u]]))
      data_landuse2016[i, value_string] <- area 
      counts_sum <- counts_sum + count
    }
  } else {
    data_landuse2016[i, "none"] <- length(ext[[i]]) * RESOLUTION * RESOLUTION
    counts_sum <- length(ext[[i]])
    }
  
  if (counts_sum < length(ext[[i]])) {
    data_landuse2016[i, "none"] <- (length(ext[[i]]) - counts_sum) * RESOLUTION * RESOLUTION
  }
}
landuse2016@data <- cbind(landuse2016@data, data_landuse2016)
print(paste0("Duration loop: ", difftime(Sys.time(), start_time_loop, units = "secs"), " sec"))

# save result
writeOGR(landuse2016, "./data_output", "landuse_2016_with_2001", driver = "ESRI Shapefile", overwrite_layer = TRUE)

# create and save summary table for comparison between 2001 and 2016
landuse2016_2001 <- as_tibble(landuse2016@data)
columns <- c("none", as.vector(classes_nominal))
landuse2016_2001 <- landuse2016_2001[, c("NUTZUNG_CO", columns)]
landuse_change <- landuse2016_2001 %>%
  group_by(NUTZUNG_CO) %>%
  summarise_all(.funs = c(sum), na.rm = TRUE) %>%
  mutate_at(columns, funs(. / 10000))
write_csv(landuse_change, "data_output/landuse_2016_2001.csv")

print(paste0("Duration TOTAL: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

