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
# poly_m und poly_l: Absturz bei readOGR geojson


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