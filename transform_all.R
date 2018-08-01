library(sf)
layers <- c("point_s", "line_s", "poly_s",
            "point_m", "line_m", "poly_m",
            "point_l", "line_l", "poly_l", "poly_2_s")
for (layer in layers) {
  shp <- st_read(paste0("data_input/", layer, ".shp"))
  shp <- st_transform(shp, 31256)
  st_write(shp, paste0("data_output/", layer, ".shp"), layer)
  rm(shp)
}