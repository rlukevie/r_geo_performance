library(config)


# --- set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# --- load test functions
source("testing.R")


# --- set rasterOptions
library(raster)
rasterOptions(maxmemory = 1e+07)
rasterOptions(datatype = "FLT4S")


# --- run tests

# 
# #############################################################################################
# #                                                                                           #
# #   3V01: Vektordaten projizieren und transformieren                                        #
# #                                                                                           #
# #############################################################################################
# 
# 
# config <- prepare_test("3V01_reproject_vectorsp_s_m")
# read_rdata(layertype = "sp", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
# 
# 
# config <- prepare_test("3V01_reproject_vectorsp_l")
# read_rdata(layertype = "sp", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("l"))
# 
# 
# config <- prepare_test("3V01_reproject_vectorsf_s_m")
# read_rdata(layertype = "sf", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
# 
# 
# config <- prepare_test("3V01_reproject_vectorsf_l")
# read_rdata(layertype = "sf", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("l"))
# 
# 
# #############################################################################################
# #                                                                                           #
# #   3V02: Vektordaten räumlich zusammenführen (Spatial Join)                                #
# #                                                                                           #
# #############################################################################################
# 
# 
# config <- prepare_test("3V02_spatial_joinsf")
# read_rdata(layertype = "sf", sizes = c("s", "m", "l"), geomtypes = c("point"))
# load("data_input/poly_2_ssf.RData")
# point_lsf <- st_transform(point_lsf, st_crs(poly_2_ssf))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
# rm(poly_2_ssf)
# 
# 
# config <- prepare_test("3V02_spatial_joinsp")
# read_rdata(layertype = "sp", sizes = c("s", "m", "l"), geomtypes = c("point"))
# load("data_input/poly_2_ssp.RData")
# poly_2_ssp <- spTransform(poly_2_ssp, CRSobj = CRS("+init=epsg:31256"))
# point_ssp <- spTransform(point_ssp, CRSobj = CRS("+init=epsg:31256"))
# point_msp <- spTransform(point_msp, CRSobj = CRS("+init=epsg:31256"))
# point_lsp <- spTransform(point_lsp, CRSobj = CRS("+init=epsg:31256"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
# rm(poly_2_ssp)
# 
# 
# config <- prepare_test("3V02_spatial_join_points_to_poly_sp")
# library(rgdal)
# library(dplyr)
# point_lsp <- readOGR("./data_input/point_l.shp", "point_l")
# poly_2_ssp <- readOGR("./data_input/poly_2_s.shp", "poly_2_s")
# join_points_to_poly_sp <- function(point, poly) {
#   poly@data <- mutate(poly@data, id_poly = as.numeric(rownames(poly@data)))
#   point@data <- mutate(point@data, id_point = as.numeric(rownames(point@data)))
#   overlay <- over(point, poly)
#   overlay <- mutate(overlay, id_point = as.numeric(rownames(overlay)))
#   overlay <- left_join(point@data, overlay, by = c("id_point" = "id_point"))
#   overlay_agg <- overlay %>%
#     group_by(id_poly) %>%
#     summarise(point_count = n()) %>%
#     arrange(id_poly)
#   poly@data <- left_join(poly@data, overlay_agg, by = c("id_poly" = "id_poly"))
#   return(poly)
# }
# test_performance_grid(config)
# rm(point_lsp)
# rm(poly_2_ssp)
# 
# 
# config <- prepare_test("3V02_spatial_join_points_to_poly_sf")
# point_lsf <- st_read("./data_input/point_l.shp")
# poly_2_ssf <- st_read("./data_input/poly_2_s.shp")
# join_points_to_poly_sf <- function(point, poly) {
#   poly <- mutate(poly, id_poly = as.numeric(rownames(poly)))
#   joined = st_join(x = point, y = poly)
#   joined_agg <- joined %>%
#     group_by(id_poly) %>%
#     summarise(point_count = n()) %>%
#     arrange(id_poly)
#   return(joined_agg)
# }
# test_performance_grid(config)
# rm(point_lsf)
# rm(poly_2_ssf)
# 
# 
# #############################################################################################
# #                                                                                           #
# #   3V04: Vektordaten verschneiden: Überschneidung (Intersect)                              #
# #                                                                                           #
# #############################################################################################
# 
# 
# config <- prepare_test("3V04_vector_intersect")
# read_rdata(layertype = "sf", sizes = c("s", "m"), geomtypes = "poly")
# read_rdata(layertype = "sp", sizes = c("s"), geomtypes = "poly")
# load("data_input/poly_2_ssf.RData")
# load("data_input/poly_2_ssp.RData")
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
# remove_layer_objects(layertype = "sp", sizes = c("s"))
# rm(poly_2_ssf)
# rm(poly_2_ssp)
# # poly_msp mit rgeos nach ca. 2 Stunden abgebrochen
# 
# config <- prepare_test("3V04_vector_intersect_raster")
# read_rdata(layertype = "sp", sizes = c("s"), geomtypes = "poly")
# load("data_input/poly_2_ssp.RData")
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s"))
# rm(poly_2_ssf)
# rm(poly_2_ssp)
# # Speicherüberlauf
# 

#############################################################################################
#                                                                                           #
#   3R01: Euklidische Distanz in Rasterdaten berechnen                                      #
#                                                                                           #
#############################################################################################


raster_distance_rqigs <- function(ras, ...) {
  raster_name <- deparse(substitute(ras))
  filename <- paste0("data_output/raster_distance_RQIGS_", raster_name, ".tif")
  set_env(root = "C:/OSGeo4W64",
          new = TRUE)
  params <- get_args_man(alg = "gdalogr:proximity")
  params$INPUT <- ras
  params$OUTPUT <- filename
  params$VALUES <- 30201
  params$RTYPE <- 5
  params$UNITS <- 0
  run_qgis(alg = "gdalogr:proximity", params = params)
}
library(raster)
config <- prepare_test("3R01_raster_distance_rqigs")
raster_streets <- raster("data_input/raster_streets_30201.tif")
raster_streets_s <- raster("data_input/raster_streets_30201_s.tif")
test_performance_grid(config)
rm(raster_streets)
rm(raster_streets_s)


raster_distance_raster <- function(ras, ...) {
  raster_name <- deparse(substitute(ras))
  filename <- paste0("data_output/raster_distance_raster_", raster_name, ".tif")
  ras_distance <- distance(ras)
  writeRaster(ras_distance, filename, overwrite = TRUE)
}
library(raster)
config <- prepare_test("3R01_raster_distance_raster")
raster_streets <- raster("data_input/raster_streets_30201.tif")
raster_streets_s <- raster("data_input/raster_streets_30201_s.tif")
test_performance_grid(config)
rm(raster_streets)
rm(raster_streets_s)


#############################################################################################
#                                                                                           #
#   3R02: Rasterdaten mit Map Algebra lokal modifizieren                                    #
#                                                                                           #
#############################################################################################


modify_raster_local_raster <- function(ras){
  fun <- function(x) { # nonsensical test function
    ifelse(x < 80, x / (x + 1) * (x + 2), x / (x + 1) * (x + 2) + 1000)
  }
  raster_local <- calc(x = ras, fun = fun)
  raster_name <- deparse(substitute(ras))
  filename <- paste0("data_output/raster_local_raster_", raster_name, ".tif")
  writeRaster(raster_local, filename, overwrite = TRUE)
}
library(raster)
config <- prepare_test("3R02_raster_local_raster")
read_raster_with_raster(format = "tif")
test_performance_grid(config)
remove_raster_objects()


modify_raster_local_spatial_tools <- function(ras, n_cpu){
  fun <- function(inraster) { # nonsensical test function
    ifelse(inraster < 80, inraster / (inraster + 1) * (inraster + 2),
           inraster / (inraster + 1) * (inraster + 2) + 1000)
  }
  start_time <- Sys.time()
  if (n_cpu > 1) {
    sfQuickInit(cpus=n_cpu)
  }
  raster_local_spatial_tools <- rasterEngine(
    inraster = ras,
    fun = fun)
  if (n_cpu > 1) {
    sfQuickStop()
  }
  raster_name <- deparse(substitute(ras))
  filename <- paste0("data_output/raster_local_spatial_tools_", raster_name, ".tif")
  writeRaster(raster_local_spatial_tools, filename, overwrite = TRUE)
}
library(raster)
config <- prepare_test("3R02_raster_local_spatial_tools")
read_raster_with_raster(format = "tif")
test_performance_grid(config)
remove_raster_objects()


modify_raster_local_rqgis <- function(ras) {
  set_env(root = "C:/OSGeo4W64",
          new = TRUE)
  params <- get_args_man(alg = "gdalogr:rastercalculator")
  params$INPUT_A <- ras
  params$FORMULA <- "where(A<80,A/(A+1)*(A+2),A/(A+1)*(A+2)+1000)"
  params$RTYPE <- 5
  raster_name <- deparse(substitute(ras))
  filename <- paste0("data_output/raster_local_rqgis_", raster_name, ".tif")
  params$OUTPUT <- filename
  raster_local_qgis_gdal <- run_qgis(alg = "gdalogr:rastercalculator",
                                     params = params,
                                     load_output = TRUE)
}
library(raster)
config <- prepare_test("3R02_raster_local_rqgis")
read_raster_with_raster(format = "tif")
test_performance_grid(config)
remove_raster_objects()


#############################################################################################
#                                                                                           #
#   3R03: Rasterdaten mit Map Algebra fokal modifizieren                                    #
#                                                                                           #
#############################################################################################


modify_raster_focal_raster <- function(ras, window_size) {
  fun <- function(x) {
    quantiles <- quantile(x, names = FALSE, na.rm = TRUE)
    return(quantiles[2] - 1.5*(quantiles[4] - quantiles[2]))
  }
  window <- matrix(data = 1, nrow = window_size, ncol = window_size)
  raster_focal_raster <- focal(x = ras, w = window, fun = fun)
  raster_name <- deparse(substitute(ras))
  filename <- paste0("data_output/raster_focal_raster_", raster_name, ".tif")
  writeRaster(raster_focal_raster, filename, overwrite = TRUE)
}
library(raster)
config <- prepare_test("3R03_raster_focal_raster")
read_raster_with_raster(format = "tif")
test_performance_grid(config)
remove_raster_objects()


modify_raster_focal_spatial_tools <- function(ras, window_size, n_cpu) {
  fun <- function(inraster) {
    quantiles <- quantile(inraster, names = FALSE, na.rm = TRUE)
    return(quantiles[2] - 1.5*(quantiles[4] - quantiles[2]))
  }
  if (n_cpu > 1) {
    sfQuickInit(cpus=n_cpu)
  }
  raster_focal_spatial_tools <- rasterEngine(
    inraster = ras,
    fun = fun,
    window_dims = c(window_size, window_size))
  if (n_cpu > 1) {
    sfQuickStop()
  }
  raster_name <- deparse(substitute(ras))
  filename <- paste0("data_output/raster_focal_spatial_tools_", raster_name, ".tif")
  writeRaster(raster_focal_spatial_tools, filename, overwrite = TRUE)
}
library(raster)
config <- prepare_test("3R03_raster_focal_spatial_tools")
read_raster_with_raster(format = "tif")
test_performance_grid(config)
remove_raster_objects()

