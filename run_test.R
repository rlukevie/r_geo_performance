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


#############################################################################################
#                                                                                           #
#   1V01: Vektordaten lesen                                                                 #
#                                                                                           #
#############################################################################################


config <- prepare_test("1V01_read_vector_serialized_sp")
test_performance_grid(config)


config <- prepare_test("1V01_read_vector_serialized_sf")
test_performance_grid(config)


config <- prepare_test("1V01_read_vector_shape")
test_performance_grid(config)


config <- prepare_test("1V01_read_vector_geojson")
test_performance_grid(config)
# poly_m und poly_l: Absturz bei readOGR geojson


config <- prepare_test("1V01_read_vector_kml")
test_performance_grid(config)


#############################################################################################
#                                                                                           #
#   1R01: Rasterdaten lesen                                                                 #
#                                                                                           #
#############################################################################################


config <- prepare_test("1R01_read_raster_geotiff")
test_performance_grid(config)


config <- prepare_test("1R01_read_raster_asc")
test_performance_grid(config)


#############################################################################################
#                                                                                           #
#   1V02: Vektordaten schreiben                                                             #
#                                                                                           #
#############################################################################################


config <- prepare_test("1V02_write_vector_shapesp_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


config <- prepare_test("1V02_write_vector_shapesp_l")
read_rdata(layertype = "sp", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


config <- prepare_test("1V02_write_vector_shapesf_s_m")
read_rdata(layertype = "sf", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


config <- prepare_test("1V02_write_vector_shapesf_l")
read_rdata(layertype = "sf", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("l"))


config <- prepare_test("1V02_write_vector_geojsonsp_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


config <- prepare_test("1V02_write_vector_geojsonsp_l")
read_rdata(layertype = "sp", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


config <- prepare_test("1V02_write_vector_geojsonsf_s_m")
read_rdata(layertype = "sf", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


config <- prepare_test("1V02_write_vector_geojsonsf_l")
read_rdata(layertype = "sf", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("l"))


config <- prepare_test("1V02_write_vector_kmlsp_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


config <- prepare_test("1V02_write_vector_kmlsp_l")
read_rdata(layertype = "sp", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


config <- prepare_test("1V02_write_vector_kmlsf_s_m")
read_rdata(layertype = "sf", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


config <- prepare_test("1V02_write_vector_kmlsf_l")
read_rdata(layertype = "sf", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("l"))

save_serialized <- function(object, fname) {
  save(object, file = fname)
}
config <- prepare_test("1V02_write_vector_serializedsp_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))

save_serialized <- function(object, fname) {
  save(object, file = fname)
}
config <- prepare_test("1V02_write_vector_serializedsp_l")
read_rdata(layertype = "sp", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


save_serialized <- function(object, fname) {
  save(object, file = fname)
}
config <- prepare_test("1V02_write_vector_serializedsf_s_m")
read_rdata(layertype = "sf", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


save_serialized <- function(object, fname) {
  save(object, file = fname)
}
config <- prepare_test("1V02_write_vector_serializedsf_l")
read_rdata(layertype = "sf", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("l"))


#############################################################################################
#                                                                                           #
#   1R02: Rasterdaten schreiben                                                             #
#                                                                                           #
#############################################################################################


config <- prepare_test("1R02_write_raster_geotiff")
read_raster_with_raster(format = "tif")
library(rgdal)
raster_sgdal <- readGDAL("data_input/raster_s.tif")
raster_mgdal <- readGDAL("data_input/raster_m.tif")
raster_lgdal <- readGDAL("data_input/raster_l.tif")
test_performance_grid(config)
remove_raster_objects()
rm(raster_sgdal)
rm(raster_mgdal)
rm(raster_lgdal)


config <- prepare_test("1R02_write_raster_asc")
read_raster_with_raster(format = "tif")
raster_sgdal <- readGDAL("data_input/raster_s.tif")
raster_mgdal <- readGDAL("data_input/raster_m.tif")
raster_lgdal <- readGDAL("data_input/raster_l.tif")
test_performance_grid(config)
remove_raster_objects()
rm(raster_sgdal)
rm(raster_mgdal)
rm(raster_lgdal)

config <- prepare_test("1R02_write_raster_asc_maptools")
raster_smaptools <- readAsciiGrid("data_input/raster_s.asc")
raster_mmaptools <- readAsciiGrid("data_input/raster_m.asc")
raster_lmaptools <- readAsciiGrid("data_input/raster_l.asc")
test_performance_grid(config)
rm(raster_smaptools)
rm(raster_mmaptools)
rm(raster_lmaptools)

config <- prepare_test("1R02_write_raster_asc_sp")
raster_ssp <- read.asciigrid("data_input/raster_s.asc")
raster_msp <- read.asciigrid("data_input/raster_m.asc")
raster_lsp <- read.asciigrid("data_input/raster_l.asc")
test_performance_grid(config)
rm(raster_ssp)
rm(raster_msp)
rm(raster_lsp)



#############################################################################################
#                                                                                           #
#   2V01: Vektorattribute modifizieren                                                      #
#                                                                                           #
#############################################################################################


config <- prepare_test("2V01_modify_vectorsp_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


config <- prepare_test("2V01_modify_vectorsp_l")
read_rdata(layertype = "sp", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


config <- prepare_test("2V01_modify_vectorsf_s_m")
read_rdata(layertype = "sf", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


config <- prepare_test("2V01_modify_vectorsf_l")
read_rdata(layertype = "sf", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("l"))


#############################################################################################
#                                                                                           #
#   2R01: Rasterdaten reklassifizieren                                                      #
#                                                                                           #
#############################################################################################


reclassify_raster_qgis_grass7 <- function(raster_name) {
  set_env(root = "C:/OSGeo4W64",
          new = TRUE)
  params <- get_args_man(alg = "grass7:r.reclass")
  params$input <- paste0("data_input/", raster_name, ".tif")
  params$txtrules <- " -200 thru 0 = 1 0.1 thru 100 = 2 100.1 thru 300 = 3 end"
  params$output <- paste0("data_output/", raster_name, "_rec_qgis_grass7.tif")
  raster_rec <- run_qgis(alg = "grass7:r.reclass",
                         params = params,
                         load_output = TRUE)
}
config <- prepare_test("2R01_reclassify_raster")
read_raster_with_raster(format = "tif")
rcl = matrix(c(-200, 0, 1, 0, 100, 2, 100, 300, 3), ncol = 3, byrow = TRUE)
test_performance_grid(config)
remove_raster_objects()
rm(rcl)



#############################################################################################
#                                                                                           #
#   3V01: Vektordaten projizieren und transformieren                                        #
#                                                                                           #
#############################################################################################


config <- prepare_test("3V01_reproject_vectorsp_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


config <- prepare_test("3V01_reproject_vectorsp_l")
read_rdata(layertype = "sp", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


config <- prepare_test("3V01_reproject_vectorsf_s_m")
read_rdata(layertype = "sf", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


config <- prepare_test("3V01_reproject_vectorsf_l")
read_rdata(layertype = "sf", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("l"))


#############################################################################################
#                                                                                           #
#   3V02: Vektordaten räumlich zusammenführen (Spatial Join)                                #
#                                                                                           #
#############################################################################################


config <- prepare_test("3V02_spatial_joinsf")
read_rdata(layertype = "sf", sizes = c("s", "m", "l"), geomtypes = c("point"))
load("data_input/poly_2_ssf.RData")
point_lsf <- st_transform(point_lsf, st_crs(poly_2_ssf))
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
rm(poly_2_ssf)


config <- prepare_test("3V02_spatial_joinsp")
read_rdata(layertype = "sp", sizes = c("s", "m", "l"), geomtypes = c("point"))
load("data_input/poly_2_ssp.RData")
poly_2_ssp <- spTransform(poly_2_ssp, CRSobj = CRS("+init=epsg:31256"))
point_ssp <- spTransform(point_ssp, CRSobj = CRS("+init=epsg:31256"))
point_msp <- spTransform(point_msp, CRSobj = CRS("+init=epsg:31256"))
point_lsp <- spTransform(point_lsp, CRSobj = CRS("+init=epsg:31256"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
rm(poly_2_ssp)


config <- prepare_test("3V02_spatial_join_points_to_poly_sp")
library(rgdal)
library(dplyr)
point_lsp <- readOGR("./data_input/point_l.shp", "point_l")
poly_2_ssp <- readOGR("./data_input/poly_2_s.shp", "poly_2_s")
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
test_performance_grid(config)
rm(point_lsp)
rm(poly_2_ssp)


config <- prepare_test("3V02_spatial_join_points_to_poly_sf")
point_lsf <- st_read("./data_input/point_l.shp")
poly_2_ssf <- st_read("./data_input/poly_2_s.shp")
join_points_to_poly_sf <- function(point, poly) {
  poly <- mutate(poly, id_poly = as.numeric(rownames(poly)))
  joined = st_join(x = point, y = poly)
  joined_agg <- joined %>%
    group_by(id_poly) %>%
    summarise(point_count = n()) %>%
    arrange(id_poly)
  return(joined_agg)
}
test_performance_grid(config)
rm(point_lsf)
rm(poly_2_ssf)


#############################################################################################
#                                                                                           #
#   3V04: Vektordaten verschneiden: Überschneidung (Intersect)                              #
#                                                                                           #
#############################################################################################


config <- prepare_test("3V04_vector_intersect")
read_rdata(layertype = "sf", sizes = c("s", "m"), geomtypes = "poly")
read_rdata(layertype = "sp", sizes = c("s"), geomtypes = "poly")
load("data_input/poly_2_ssf.RData")
load("data_input/poly_2_ssp.RData")
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
remove_layer_objects(layertype = "sp", sizes = c("s"))
rm(poly_2_ssf)
rm(poly_2_ssp)

config <- prepare_test("3V04_vector_intersect_raster")
read_rdata(layertype = "sp", sizes = c("s"), geomtypes = "poly")
load("data_input/poly_2_ssp.RData")
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s"))
rm(poly_2_ssf)
rm(poly_2_ssp)


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



#############################################################################################
#                                                                                           #
#   4RV01: Vektordaten rasterisieren                                                        #
#                                                                                           #
#############################################################################################


rasterize_raster <- function(vector_sp, resolution, field, ...) {
  vector_sp_name <- paste0(deparse(substitute(vector_sp)), "_", resolution)
  filename <- paste0("data_output/rasterize_raster__", vector_sp_name, ".tif")
  vector_sp@data[[field]] <- as.numeric(as.character(vector_sp@data[[field]]))
  raster <- raster(vector_sp, res = resolution)
  raster <- raster::rasterize(x = vector_sp, y = raster, field = field, filename = filename, overwrite = TRUE)
}
config <- prepare_test("4RV01_rasterize_raster_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
config <- prepare_test("4RV01_rasterize_raster_l")
read_rdata(layertype = "sp", sizes = c("l"))
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


rasterize_fasterize <- function(vector_sf, resolution, field, ...) {
  vector_sf_name <- paste0(deparse(substitute(vector_sf)), "_", resolution)
  filename <- paste0("data_output/rasterize_fasterize__", vector_sf_name, ".tif")
  raster <- raster(vector_sf, res = resolution)
  raster <- fasterize::fasterize(sf = vector_sf, raster = raster, field = field, ...)
  writeRaster(raster, filename, overwrite = TRUE)
}
config <- prepare_test("4RV01_rasterize_fasterize_s_m")
read_rdata(layertype = "sf", sizes = c("s", "m"), geomtypes = c("poly"))
library(raster)
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
config <- prepare_test("4RV01_rasterize_fasterize_l")
read_rdata(layertype = "sf", sizes = c("l"), geomtypes = c("poly"))
library(raster)
test_performance_grid(config)
remove_layer_objects(layertype = "sf", sizes = c("l"))


rasterize_velox <- function(vector_sp, resolution, field, ...) {
  vector_sp_name <- paste0(deparse(substitute(vector_sp)), "_", resolution)
  filename <- paste0("data_output/rasterize_velox__", vector_sp_name, ".tif")
  vector_sp@data[[field]] <- as.numeric(as.character(vector_sp@data[[field]]))
  extent <- extent(vector_sp)
  crs <- crs(vector_sp)
  raster <- raster(vector_sp, res = resolution, val = 99999999)
  raster_vx <- velox(raster)
  raster_vx$rasterize(vector_sp, band = 1, field = field, ...)
  raster <- raster(x = raster_vx$rasterbands[[1]],
                   xmn = extent@xmin,
                   xmx = extent@xmax,
                   ymn = extent@ymin,
                   ymx = extent@ymax,
                   crs = crs)
  writeRaster(raster, filename, NAflag = 99999999, overwrite = TRUE)
}
config <- prepare_test("4RV01_rasterize_velox_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"), geomtypes = c("poly", "line"))
library(raster)
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


config <- prepare_test("4RV01_rasterize_velox_l")
read_rdata(layertype = "sp", sizes = c("l"), geomtypes = c("poly", "line"))
library(raster)
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


config <- prepare_test("4RV01_rasterize_gdalutils")
test_performance_grid(config)




rasterize_rqgis_grass7 <- function(vector_sp, resolution, field) {
  set_env(root = "C:/OSGeo4W64",
          new = TRUE)
  params <- get_args_man(alg = "grass7:v.to.rast.attribute")
  params$input <- vector_sp
  params$attribute_column = field
  params$GRASS_REGION_CELLSIZE_PARAMETER = resolution
  params$output <- "data_output/rasterize_rqgis_grass7.tif"
  print(params)
  rasterized <- run_qgis(alg = "grass7:v.to.rast.attribute",
                         params = params,
                         load_output = TRUE)
}
config <- prepare_test("4RV01_rasterize_rqgis_grass7_s_m")
read_rdata(layertype = "sp", sizes = c("s", "m"), geomtypes = c("poly", "line", "point"))
library(RQGIS)
library(raster)
library(rgdal)
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
# Unfortunately, QGIS did not produce: C:/Users/rolandstudio/AIT/50_performance_messung/data_output/rasterize_rqgis_grass7.tif
config <- prepare_test("4RV01_rasterize_rqgis_grass7_l")
read_rdata(layertype = "sp", sizes = c("l"), geomtypes = c("poly", "line", "point"))
library(RQGIS)
library(raster)
test_performance_grid(config)
remove_layer_objects(layertype = "sp", sizes = c("l"))


#############################################################################################
#                                                                                           #
#   4RV02: Raster zuschneiden (maskieren)                                                   #
#                                                                                           #
#############################################################################################


raster_mask_raster <- function(raster, mask) {
  raster_name <- deparse(substitute(raster))
  mask_name <- deparse(substitute(mask))
  filename <- paste0("data_output/raster_mask_raster_", raster_name, "_", mask_name, ".tif")
  raster <- crop(raster, mask)
  crs(mask) <- crs(raster)
  raster_mask_raster <- mask(raster, mask)
  writeRaster(raster_mask_raster, filename, overwrite = TRUE)
}
library(raster)
library(rgdal)
config <- prepare_test("4RV02_raster_mask_raster")
raster_s <- raster("data_input/raster_s.tif")
raster_m <- raster("data_input/raster_m.tif")
raster_l <- raster("data_input/raster_l.tif")
mask_ssp <- readOGR("data_input/crop_s.shp", "crop_s")
mask_msp <- readOGR("data_input/crop_m.shp", "crop_m")
mask_lsp <- readOGR("data_input/crop_l.shp", "crop_l")
test_performance_grid(config)
rm(raster_s)
rm(raster_m)
rm(raster_l)
rm(mask_ssp)
rm(mask_msp)
rm(mask_lsp)


raster_mask_raster_fasterize <- function(raster, mask) {
  raster_name <- deparse(substitute(raster))
  mask_name <- deparse(substitute(mask_sp))
  filename <- paste0("data_output/raster_mask_raster_", raster_name, "_", mask_name, ".tif")
  raster <- crop(raster, extent(as(mask, "Spatial")))
  mask_raster <- fasterize::fasterize(sf = mask, raster = raster)
  crs(mask_raster) <- crs(raster)
  raster_mask_raster <- mask(raster, mask_raster)
  writeRaster(raster_mask_raster, filename, overwrite = TRUE)
}
library(raster)
library(sf)
config <- prepare_test("4RV02_raster_mask_raster_fasterize")
raster_s <- raster("data_input/raster_s.tif")
raster_m <- raster("data_input/raster_m.tif")
raster_l <- raster("data_input/raster_l.tif")
mask_ssf <- st_read("data_input/crop_s.shp")
mask_msf <- st_read("data_input/crop_m.shp")
mask_lsf <- st_read("data_input/crop_l.shp")
test_performance_grid(config)
rm(raster_s)
rm(raster_m)
rm(raster_l)
rm(mask_ssf)
rm(mask_msf)
rm(mask_lsf)


#############################################################################################
#                                                                                           #
#   4RV03: Rasterdaten auf Basis von Punktdaten extrahieren                                 #
#                                                                                           #
#############################################################################################


config <- prepare_test("4RV03_raster_extract_point")
library(rgdal)
grid <- readGDAL("data_input/raster_landuse200708.tif")
raster <- raster("data_input/raster_landuse200708.tif")
load("data_input/point_lsp.RData")
test_performance_grid(config)
rm(grid)
rm(point_lsp)


#############################################################################################
#                                                                                           #
#   4RV04: Rasterdaten auf Basis von Polygondaten extrahieren                               #
#                                                                                           #
#############################################################################################


raster_extract_polygon_velox <- function(ras, poly) {
  ras_vx <- velox(ras)
  ext <- ras_vx$extract(poly)
  data_poly <- data.frame()
  for (i in 1:length(ext)) {
    tab <- table(ext[[i]])
    unique_values_one_poly <- attr(tab, "names")
    counts_one_poly <- as.vector(tab)
    nrows <- nrow(tab)
    for (u in 1:nrows) {
      value_string <- toString(unique_values_one_poly[[u]])
      count <- counts_one_poly[[u]]
      data_poly[i, value_string] <- count
    }
  }
  poly@data <- cbind(poly@data, data_poly)
  writeOGR(poly, "./data_output", "raster_extract_polygon_velox", driver = "ESRI Shapefile", overwrite_layer = TRUE)
}
library(sf)
library(rgdal)
library(raster)
library(fasterize)
config <- prepare_test("4RV04_raster_extract_polygon_velox")
sprengel <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
res = 10
landuse <- st_read("data_input/REALNUT200708OGD/REALNUT200708OGDPolygon.shp")
landuse <- st_transform(landuse, 31256)
landuse_raster <- raster(landuse, res = res, val = 1)
landuse_raster <- fasterize(sf = landuse, raster = landuse_raster, field = "NUTZUNG_CO")
test_performance_grid(config)
rm(sprengel)
rm(landuse)
rm(landuse_raster)


raster_extract_polygon_raster <- function(ras, poly) {
  ext <- extract(ras, poly)
  data_poly <- data.frame()
  for (i in 1:length(ext)) {
    tab <- table(ext[[i]])
    unique_values_one_poly <- attr(tab, "names")
    counts_one_poly <- as.vector(tab)
    nrows <- nrow(tab)
    for (u in 1:nrows) {
      value_string <- toString(unique_values_one_poly[[u]])
      count <- counts_one_poly[[u]]
      data_poly[i, value_string] <- count
    }
  }
  poly@data <- cbind(poly@data, data_poly)
  writeOGR(poly, "./data_output", "raster_extract_polygon_raster", driver = "ESRI Shapefile", overwrite_layer = TRUE)
}
library(sf)
library(rgdal)
library(raster)
library(fasterize)
config <- prepare_test("4RV04_raster_extract_polygon_raster")
sprengel <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
res = 10
landuse <- st_read("data_input/REALNUT200708OGD/REALNUT200708OGDPolygon.shp")
landuse <- st_transform(landuse, 31256)
landuse_raster <- raster(landuse, res = res, val = 1)
landuse_raster <- fasterize(sf = landuse, raster = landuse_raster, field = "NUTZUNG_CO")
test_performance_grid(config)
rm(sprengel)
rm(landuse)
rm(landuse_raster)