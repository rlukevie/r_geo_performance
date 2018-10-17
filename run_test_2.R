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
# #   2V01: Vektorattribute modifizieren                                                      #
# #                                                                                           #
# #############################################################################################
# 
# 
# config <- prepare_test("2V01_modify_vectorsp_s_m")
# read_rdata(layertype = "sp", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
# 
# 
# config <- prepare_test("2V01_modify_vectorsp_l")
# read_rdata(layertype = "sp", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("l"))
# 
# 
# config <- prepare_test("2V01_modify_vectorsf_s_m")
# read_rdata(layertype = "sf", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
# 
# 
# config <- prepare_test("2V01_modify_vectorsf_l")
# read_rdata(layertype = "sf", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("l"))
# 

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

