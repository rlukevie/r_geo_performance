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