################################
# 3V02: SPATIAL JOIN VECTOR    #
################################

# === Example: Join polygons to points  ===

# sp              --> YES, in combination with rgeos (0.8 sec)
# sf              --> YES (1.5 sec)
# rgeos           --> YES, in combination with sp (0.8 sec)
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES, but not tested because of good performance of other (standard) packages

# fasterize       --> NO
# ----------------------------------


# sp
library(sp)
library(rgdal)
load("data_input/poly_2_ssp.RData")
load("data_input/point_lsp.RData")
poly_2_ssp <- spTransform(poly_2_ssp, CRSobj = CRS("+init=epsg:31256"))
point_lsp <- spTransform(point_lsp, CRSobj = CRS("+init=epsg:31256"))
start_time <- Sys.time()
point_lsp@data <- cbind(point_lsp@data, over(x = point_lsp, y = poly_2_ssp))
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))


# sf
library(sf)
load("data_input/poly_2_ssf.RData")
load("data_input/point_lsf.RData")

poly_2_ssf <- st_transform(poly_2_ssf, st_crs(point_lsf))
start_time <- Sys.time()
point_lsf_joined <- st_join(x = point_lsf, y = poly_2_ssf)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))



# === Example: Join Points to Polygon and count number of points (custom function) ===

# sp              --> YES, in combination with rgeos (1.1 sec)
# sf              --> YES (4.2 sec)
# rgeos           --> YES, in combination with sp (1.1 sec)
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES, but not tested because of good performance of other (standard) packages

# fasterize       --> NO
# ----------------------------------

# sp

library(rgdal)
library(sp)
library(dplyr)

load("data_input/poly_2_ssp.RData")
load("data_input/point_lsp.RData")
poly_2_ssp <- spTransform(poly_2_ssp, CRSobj = CRS("+init=epsg:31256"))
point_lsp <- spTransform(point_lsp, CRSobj = CRS("+init=epsg:31256"))
join_points_to_poly_sp <- function(point, poly) {
  poly@data <- mutate(poly@data, id_poly = as.numeric(rownames(poly@data)))
  point@data <- mutate(point@data, id_point = as.numeric(rownames(point@data)))
  overlay <- over(point, poly)
  overlay <- mutate(overlay, id_point = as.numeric(rownames(overlay)))
  overlay <- left_join(point@data, overlay, by = c("id_point" = "id_point"))
  overlay_agg <- overlay %>%
    group_by(id_poly) %>%
    summarise(point_count = n()) %>%
    arrange(id_poly)
  poly@data <- left_join(poly@data, overlay_agg, by = c("id_poly" = "id_poly"))
  return(poly)
}
start_time <- Sys.time()
join <- join_points_to_poly_sp(point_lsp, poly_2_ssp)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))


# sf
library(sf)
library(dplyr)

load("data_input/poly_2_ssf.RData")
load("data_input/point_lsf.RData")

poly_2_ssf <- st_transform(poly_2_ssf, st_crs(point_lsf))
join_points_to_poly_sf <- function(point, poly) {
  poly <- mutate(poly, id_poly = as.numeric(rownames(poly)))
  joined = st_join(x = point, y = poly)
  joined_agg <- joined %>%
    group_by(id_poly) %>%
    summarise(point_count = n()) %>%
    arrange(id_poly)
  return(joined_agg)
}

start_time <- Sys.time()
join_sf <- join_points_to_poly_sf(point = point_lsf, poly = poly_2_ssf)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))





