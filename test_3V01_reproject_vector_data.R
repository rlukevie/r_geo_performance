################################
# 3V01: REPROJECT VECTOR DATA  #
################################


# sp              --> YES, slow (66 sec)
# sf              --> YES, fast (6.6 sec)
# rgeos           --> NO
# maptools        --> NO
# geojsonio       --> NO
# rmapshaper      --> NO

# raster          --> NO
# spatial.tools   --> YES, slow (62 sec)
# gdalUtils       --> YES, slow, with ogr2ogr, but without loading as R object (66 sec)
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> YES (in combination with sp)
# RQGIS           --> YES, but crashes

# fasterize       --> NO
# ----------------------------------

# sp
library(rgdal)
library(sp)
load("data_input/poly_msp.RData")
start_time <- Sys.time()
poly_msp <- spTransform(x = poly_msp, CRSobj = CRS("+init=epsg:4326"))
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_msp)

# sf
library(sf)
load("data_input/poly_msf.RData")
start_time <- Sys.time()
poly_msf <- st_transform(x = poly_msf, crs = 4326)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_msf)

# spatial.tools
library(spatial.tools)
library(sp)
load("data_input/poly_msp.RData")
start_time <- Sys.time()
ref <- SpatialPoints(matrix(c(1,1), nrow = 1, ncol = 2), proj4string = CRS("+init=epsg:4326"))
poly_msp <- spatial_sync_vector(poly_msp, ref)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
rm(poly_msp)

# gdalUtils
library(gdalUtils)
start_time <- Sys.time()
ogr2ogr(src_datasource_name = "./data_input/poly_m.shp",
        dst_datasource_name = "./data_output/poly_m_transf_gdal.shp",
        t_srs = "EPSG:4326")
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

# RQGIS
library(RQGIS)
load("data_input/poly_msp.RData")
set_env(root = "C:/OSGeo4W64",
        new = TRUE)
find_algorithms("project")
get_usage("qgis:reprojectlayer")
get_options("qgis:reprojectlayer")
start_time <- Sys.time()
params <- get_args_man(alg = "qgis:reprojectlayer")
params$INPUT <- poly_msp
params$TARGET_CRS <- "4326"
params$OUTPUT <- "data_output/poly_msp_transf_qgis.shp"
poly_msp <- run_qgis(alg = "qgis:reprojectlayer",
                       params = params,
                       load_output = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))
  