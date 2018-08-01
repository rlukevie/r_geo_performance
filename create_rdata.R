

layers <- list("point_s", "line_s", "poly_s",
               "point_m", "line_m", "poly_m",
               "point_l", "line_l", "poly_l")


layers <- list("poly_l")

library(sf)
for (layer in layers) {
  eval(parse(text = paste0(layer, "sf <- st_read('data_input/", layer, ".shp')")))
  eval(parse(text = paste0("save(", layer, "sf, file = 'data_input/", layer, "sf.RData')")))
  eval(parse(text = paste0("rm(", layer, "sf)")))
}

library(sp)
library(rgdal)
for (layer in layers) {
  eval(parse(text = paste0(layer, "sp <- readOGR(dsn = './data_input/", layer, ".shp', layer = '", layer, "')")))
  eval(parse(text = paste0("save(", layer, "sp, file = 'data_input/", layer, "sp.RData')")))
  eval(parse(text = paste0("rm(", layer, "sp)")))
}


