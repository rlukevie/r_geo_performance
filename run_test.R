library(config)

# GLOBAL VARIABLES
WRITE_OUTPUT = TRUE


# --- define helper functions
load_packages <- function(config) {
  packages <- unlist(sapply(config, `[`, "package"), use.names = FALSE)
  packages <- packages[!is.na(packages)]
  lapply(packages, library, character.only = TRUE)
}

prepare_test <- function(testcase) {
  Sys.setenv(R_CONFIG_ACTIVE = testcase)
  config <- config::get()
  load_packages(config)
  gc()
  Sys.sleep(0.2)
  return(config)
}

read_rdata <- function(layertype, sizes = "all", geomtypes = "all") {
  if (sizes[1] == "all" & layertype == "sf") {
    if (geomtypes == "all") {
      layers <- c("point_s", "line_s", "poly_s",
                  "point_m", "line_m", "poly_m",
                  "point_l", "line_l", "poly_l")
    } else {
      layers <- c()
      for (geomtype in geomtypes) {
        layers <- c(layers, paste0(geomtype, "_s"))
        layers <- c(layers, paste0(geomtype, "_m"))
        layers <- c(layers, paste0(geomtype, "_l"))
      }
    }

  } else if (sizes[1] == "all" & layertype == "sp") {
    if (geomtypes[1] == "all") {
      layers <- c("point_s", "line_s", "poly_s",
                  "point_m", "line_m", "poly_m",
                  "point_l", "line_l")
    } else {
      for (geomtype in geomtypes) {
        layers <- c(layers, paste0(geomtype, "_s"))
        layers <- c(layers, paste0(geomtype, "_m"))
        if (geomtype != "poly") {
          layers <- c(layers, paste0(geomtype, "_l"))
        }
      }
    }

  } else {
    layers <- c()
    for (size in sizes) {
      if (geomtypes[1] == "all") {
        layers <- c(layers, paste0("point_", size))
        layers <- c(layers, paste0("line_", size))
        if (!(size == "l" & layertype == "sp")) {layers <- c(layers, paste0("poly_", size))}
      } else {
        for (geomtype in geomtypes) {
          if (!(size == "l" & layertype == "sp")) {
            layers <- c(layers, paste0(geomtype, "_", size))
          }
        }
      }

    }
  }
  print(layers)
  for (layer in layers) {
    load(paste0("data_input/", layer, layertype, ".RData"))
    eval(parse(text = paste0(layer, layertype, " <<- ", layer, layertype)))
  }
}

remove_layer_objects <- function(layertype, sizes = "all") {
  if (sizes[1] == "all") {
    layers <- c("point_s", "line_s", "poly_s",
                "point_m", "line_m", "poly_m",
                "point_l", "line_l", "poly_l")
  } else {
    layers <- c()
    for (size in sizes) {
      layers <- c(layers, paste0("point_", size))
      layers <- c(layers, paste0("line_", size))
      layers <- c(layers, paste0("poly_", size))
    }
  }
  for (layer in layers) {
    eval(parse(text = paste0("rm(", layer, layertype, ", pos = '.GlobalEnv')")))
  }
}


read_raster_with_raster <- function(format) {
  layers <- list("raster_s", "raster_m", "raster_l")
  for (layer in layers) {
    eval(parse(text = paste0(layer, " <<- raster('data_input/", layer, ".", format,"')")))
  }
}

remove_raster_objects <- function() {
  layers <- list("raster_s", "raster_m", "raster_l")
  for (layer in layers) {
    eval(parse(text = paste0("rm(", layer, ", pos = '.GlobalEnv')")))
  }
}


# --- set working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# --- load test functions
source("test_functions.R")

# --- run tests


#############################################################################################
#                                                                                           #
#   Vektordaten lesen                                                                       #
#                                                                                           #
#############################################################################################


# config <- prepare_test("read_vector_serialized_sp")
# test_performance_grid(config)

# config <- prepare_test("read_vector_serialized_sf")
# test_performance_grid(config)


# config <- prepare_test("read_vector_shape")
# test_performance_grid(config)


# config <- prepare_test("read_vector_geojson")
# test_performance_grid(config)


# config <- prepare_test("read_vector_kml")
# test_performance_grid(config)


#############################################################################################
#                                                                                           #
#   Rasterdaten lesen                                                                       #
#                                                                                           #
#############################################################################################


# config <- prepare_test("read_raster_geotiff")
# test_performance_grid(config)


# config <- prepare_test("read_raster_asc")
# test_performance_grid(config)


#############################################################################################
#                                                                                           #
#   Vektordaten schreiben                                                                   #
#                                                                                           #
#############################################################################################


# config <- prepare_test("write_vector_shapesp_s_m")
# read_rdata(layertype = "sp", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


# config <- prepare_test("write_vector_shapesp_l")
# read_rdata(layertype = "sp", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("l"))
# 

# config <- prepare_test("write_vector_shapesf_s_m")
# read_rdata(layertype = "sf", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
# 

# config <- prepare_test("write_vector_shapesf_l")
# read_rdata(layertype = "sf", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("l"))


#############################################################################################
#                                                                                           #
#   Rasterdaten schreiben                                                                   #
#                                                                                           #
#############################################################################################


# config <- prepare_test("write_raster_geotiff")
# read_raster_with_raster(format = "tif")
# test_performance_grid(config)
# remove_raster_objects()


# config <- prepare_test("write_raster_asc")
# read_raster_with_raster(format = "asc")
# test_performance_grid(config)
# remove_raster_objects()


#############################################################################################
#                                                                                           #
#   Vektorattribute modifizieren                                                            #
#                                                                                           #
#############################################################################################


# config <- prepare_test("modify_vectorsp_s_m")
# read_rdata(layertype = "sp", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


# config <- prepare_test("modify_vectorsp_l")
# read_rdata(layertype = "sp", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("l"))


# config <- prepare_test("modify_vectorsf_s_m")
# read_rdata(layertype = "sf", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


# config <- prepare_test("modify_vectorsf_l")
# read_rdata(layertype = "sf", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("l"))


#############################################################################################
#                                                                                           #
#   Rasterdaten reklassifizieren                                                            #
#                                                                                           #
#############################################################################################


# config <- prepare_test("reclassify_raster_tif")
# read_raster_with_raster(format = "tif")
# rcl = matrix(c(-200, 0, 1, 0, 100, 2, 100, 300, 3), ncol = 3, byrow = TRUE)
# test_performance_grid(config)
# remove_raster_objects()
# rm(rcl)


# config <- prepare_test("reclassify_raster_asc")
# read_raster_with_raster(format = "asc")
# rcl = matrix(c(-200, 0, 1, 0, 100, 2, 100, 300, 3), ncol = 3, byrow = TRUE)
# test_performance_grid(config)
# remove_raster_objects()
# rm(rcl)


#############################################################################################
#                                                                                           #
#   Vektordaten projizieren und transformieren                                              #
#                                                                                           #
#############################################################################################


# config <- prepare_test("reproject_vectorsp_s_m")
# read_rdata(layertype = "sp", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))


# config <- prepare_test("reproject_vectorsp_l")
# read_rdata(layertype = "sp", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("l"))


# config <- prepare_test("reproject_vectorsf_s_m")
# read_rdata(layertype = "sf", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


# config <- prepare_test("reproject_vectorsf_l")
# read_rdata(layertype = "sf", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("l"))


#############################################################################################
#                                                                                           #
#   Vektordaten räumlich zusammenführen (Spatial Join)                                      #
#                                                                                           #
#############################################################################################


# config <- prepare_test("spatial_joinsf")
# read_rdata(layertype = "sf", sizes = c("s", "m"))
# point_lsf <- st_read("data_input/point_l.shp")
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))

# config <- prepare_test("spatial_joinsp")
# read_rdata(layertype = "sp", sizes = c("s", "m"))
# point_lsp <- readOGR("data_input/point_l.shp", "point_l")
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))

# config <- prepare_test("spatial_join_points_to_poly_sp")
# read_rdata(layertype = "sp", sizes = c("s", "m"))
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
#   print(max(poly@data$point_count, na.rm = TRUE))
#   return(poly)
# }
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))

# 
# config <- prepare_test("spatial_join_points_to_poly_sf")
# read_rdata(layertype = "sf", sizes = c("s", "m"))
# point_lsf <- st_read("./data_input/point_l.shp")
# poly_2_ssf <- st_read("./data_input/poly_2_s.shp")
# join_points_to_poly_sf <- function(point, poly) {
#   poly <- mutate(poly, id_poly = as.numeric(rownames(poly)))
#   joined = st_join(x = point, y = poly)
#   joined_agg <- joined %>% 
#     group_by(id_poly) %>%
#     summarise(point_count = n()) %>%
#     arrange(id_poly)
#   print(max(subset(joined_agg, !is.na(joined_agg$id_poly))$point_count))
#   return(joined_agg)
# }
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))


#############################################################################################
#                                                                                           #
#   Vektordaten rasterisieren                                                               #
#                                                                                           #
#############################################################################################


# rasterize_raster <- function(vector_sp, resolution, field, ...) {
#   vector_sp_name <- paste0(deparse(substitute(vector_sp)), "_", resolution)
#   filename <- paste0("data_output/rasterize_raster__", vector_sp_name, ".tif")
#   vector_sp@data[[field]] <- as.numeric(as.character(vector_sp@data[[field]]))
#   raster <- raster(vector_sp, res = resolution)
#   raster <- raster::rasterize(x = vector_sp, y = raster, field = field, filename = filename, overwrite = TRUE)
# }
# config <- prepare_test("rasterize_raster_s_m")
# read_rdata(layertype = "sp", sizes = c("s", "m"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
# config <- prepare_test("rasterize_raster_l")
# read_rdata(layertype = "sp", sizes = c("l"))
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("l"))

# rasterize_fasterize <- function(vector_sf, resolution, field, ...) {
#   vector_sf_name <- paste0(deparse(substitute(vector_sf)), "_", resolution)
#   filename <- paste0("data_output/rasterize_fasterize__", vector_sf_name, ".tif")
#   raster <- raster(vector_sf, res = resolution)
#   raster <- fasterize::fasterize(sf = vector_sf, raster = raster, field = field, ...)
#   writeRaster(raster, filename, overwrite = TRUE)
# }
# config <- prepare_test("rasterize_fasterize_s_m")
# read_rdata(layertype = "sf", sizes = c("s", "m"), geomtypes = c("poly"))
# library(raster)
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
# config <- prepare_test("rasterize_fasterize_l")
# read_rdata(layertype = "sf", sizes = c("l"), geomtypes = c("poly"))
# library(raster)
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("l"))
# 
# 
# rasterize_velox <- function(vector_sp, resolution, field, ...) {
#   vector_sp_name <- paste0(deparse(substitute(vector_sp)), "_", resolution)
#   filename <- paste0("data_output/rasterize_velox__", vector_sp_name, ".tif")
#   vector_sp@data[[field]] <- as.numeric(as.character(vector_sp@data[[field]]))
#   extent <- extent(vector_sp)
#   crs <- crs(vector_sp)
#   raster <- raster(vector_sp, res = resolution, val = 99999999)
#   raster_vx <- velox(raster)
#   raster_vx$rasterize(vector_sp, band = 1, field = field, ...)
#   raster <- raster(x = raster_vx$rasterbands[[1]],
#                                 xmn = extent@xmin,
#                                 xmx = extent@xmax,
#                                 ymn = extent@ymin,
#                                 ymx = extent@ymax,
#                                 crs = crs)
#   writeRaster(raster, filename, NAflag = 99999999, overwrite = TRUE)
# }
# config <- prepare_test("rasterize_velox_s_m")
# read_rdata(layertype = "sp", sizes = c("s", "m"), geomtypes = c("poly", "line"))
# library(raster)
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
# config <- prepare_test("rasterize_velox_l")
# read_rdata(layertype = "sp", sizes = c("l"), geomtypes = c("poly", "line"))
# library(raster)
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("l"))
# 


#############################################################################################
#                                                                                           #
#   Euklidische Distanz in Rasterdaten berechnen                                            #
#                                                                                           #
#############################################################################################


# raster_distance_rqigs <- function(ras, ...) {
#   raster_name <- deparse(substitute(ras))
#   filename <- paste0("data_output/raster_distance_RQIGS_", raster_name, ".tif")
#   print(filename)
#   set_env(root = "C:/OSGeo4W64",
#           new = TRUE)
#   params <- get_args_man(alg = "gdalogr:proximity")
#   params$INPUT <- ras
#   params$OUTPUT <- filename
#   params$VALUES <- 30201
#   params$RTYPE <- 5
#   params$UNITS <- 0
#   run_qgis(alg = "gdalogr:proximity", params = params)
# }
# library(raster)
# config <- prepare_test("raster_distance_rqigs")
# raster_streets <- raster("data_input/raster_streets_30201.tif")
# raster_streets_s <- raster("data_input/raster_streets_30201_s.tif")
# test_performance_grid(config)
# rm(raster_streets)
# rm(raster_streets_s)

# raster_distance_raster <- function(ras, ...) {
#   raster_name <- deparse(substitute(ras))
#   filename <- paste0("data_output/raster_distance_raster_", raster_name, ".tif")
#   ras_distance <- distance(ras)
#   writeRaster(ras_distance, filename)
# }
# library(raster)
# config <- prepare_test("raster_distance_raster")
# raster_streets <- raster("data_input/raster_streets_30201.tif")
# raster_streets_s <- raster("data_input/raster_streets_30201_s.tif")
# test_performance_grid(config)
# rm(raster_streets)
# rm(raster_streets_s)


#############################################################################################
#                                                                                           #
#   Raster zuschneiden (maskieren)                                                          #
#                                                                                           #
#############################################################################################

# raster_mask_raster <- function(raster, mask) {
#   raster_name <- deparse(substitute(raster))
#   mask_name <- deparse(substitute(mask))
#   filename <- paste0("data_output/raster_mask_raster_", raster_name, "_", mask_name, ".tif")
#   raster <- crop(raster, mask)
#   crs(mask) <- crs(raster)
#   raster_mask_raster <- mask(raster, mask)
#   writeRaster(raster_mask_raster, filename, overwrite = TRUE)
# }
# library(raster)
# library(rgdal)
# config <- prepare_test("raster_mask_raster")
# raster_s <- raster("data_input/raster_s.tif")
# raster_m <- raster("data_input/raster_m.tif")
# raster_l <- raster("data_input/raster_l.tif")
# mask_ssp <- readOGR("data_input/crop_s.shp", "crop_s")
# mask_msp <- readOGR("data_input/crop_m.shp", "crop_m")
# mask_lsp <- readOGR("data_input/crop_l.shp", "crop_l")
# test_performance_grid(config)
# rm(raster_s)
# rm(raster_m)
# rm(raster_l)
# rm(mask_ssp)
# rm(mask_msp)
# rm(mask_lsp)


# raster_mask_raster_fasterize <- function(raster, mask) {
#   raster_name <- deparse(substitute(raster))
#   mask_name <- deparse(substitute(mask_sp))
#   filename <- paste0("data_output/raster_mask_raster_", raster_name, "_", mask_name, ".tif")
#   raster <- crop(raster, extent(as(mask, "Spatial")))
#   mask_raster <- fasterize::fasterize(sf = mask, raster = raster)
#   crs(mask_raster) <- crs(raster)
#   raster_mask_raster <- mask(raster, mask_raster)
#   writeRaster(raster_mask_raster, filename, overwrite = TRUE)
# }
# library(raster)
# library(sf)
# config <- prepare_test("raster_mask_raster_fasterize")
# raster_s <- raster("data_input/raster_s.tif")
# raster_m <- raster("data_input/raster_m.tif")
# raster_l <- raster("data_input/raster_l.tif")
# mask_ssf <- st_read("data_input/crop_s.shp")
# mask_msf <- st_read("data_input/crop_m.shp")
# mask_lsf <- st_read("data_input/crop_l.shp")
# test_performance_grid(config)
# rm(raster_s)
# rm(raster_m)
# rm(raster_l)
# rm(mask_ssf)
# rm(mask_msf)
# rm(mask_lsf)


#############################################################################################
#                                                                                           #
#   Rasterdaten auf Basis von Polygondaten extrahieren                                      #
#                                                                                           #
#############################################################################################

# raster_extract_polygon_velox <- function(ras, poly) {
#   ras_vx <- velox(ras)
#   ext <- ras_vx$extract(poly)
#   data_poly <- data.frame()
#   for (i in 1:length(ext)) {
#     tab <- table(ext[[i]])
#     unique_values_one_poly <- attr(tab, "names")
#     counts_one_poly <- as.vector(tab)
#     nrows <- nrow(tab)
#     for (u in 1:nrows) {
#       value_string <- toString(unique_values_one_poly[[u]])
#       count <- counts_one_poly[[u]]
#       data_poly[i, value_string] <- count
#     }
#   }
#   poly@data <- cbind(poly@data, data_poly)
#   writeOGR(poly, "./data_output", "raster_extract_polygon_velox", driver = "ESRI Shapefile", overwrite_layer = TRUE)
# }
# library(sf)
# library(rgdal)
# library(raster)
# library(fasterize)
# config <- prepare_test("raster_extract_polygon_velox")
# sprengel <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
# res = 10
# landuse <- st_read("data_input/REALNUT200708OGD/REALNUT200708OGDPolygon.shp")
# landuse <- st_transform(landuse, 31256)
# landuse_raster <- raster(landuse, res = res, val = 1)
# landuse_raster <- fasterize(sf = landuse, raster = landuse_raster, field = "NUTZUNG_CO")
# test_performance_grid(config)
# rm(sprengel)
# rm(landuse)
# rm(landuse_raster)



# raster_extract_polygon_raster <- function(ras, poly) {
#   ext <- extract(ras, poly)
#   data_poly <- data.frame()
#   for (i in 1:length(ext)) {
#     tab <- table(ext[[i]])
#     unique_values_one_poly <- attr(tab, "names")
#     counts_one_poly <- as.vector(tab)
#     nrows <- nrow(tab)
#     for (u in 1:nrows) {
#       value_string <- toString(unique_values_one_poly[[u]])
#       count <- counts_one_poly[[u]]
#       data_poly[i, value_string] <- count
#     }
#   }
#   poly@data <- cbind(poly@data, data_poly)
   # if (WRITE_OUTPUT) {
   #   writeOGR(poly, "./data_output", "raster_extract_polygon_raster", driver = "ESRI Shapefile", overwrite_layer = TRUE)
   #   }
# library(sf)
# library(rgdal)
# library(raster)
# library(fasterize)
# config <- prepare_test("raster_extract_polygon_raster")
# sprengel <- readOGR("data_input/poly_2_s.shp", "poly_2_s")
# res = 10
# landuse <- st_read("data_input/REALNUT200708OGD/REALNUT200708OGDPolygon.shp")
# landuse <- st_transform(landuse, 31256)
# landuse_raster <- raster(landuse, res = res, val = 1)
# landuse_raster <- fasterize(sf = landuse, raster = landuse_raster, field = "NUTZUNG_CO")
# test_performance_grid(config)
# rm(sprengel)
# rm(landuse)
# rm(landuse_raster)


#############################################################################################
#                                                                                           #
#   Vektordaten verschneiden: Überschneidung (Intersect)                                    #
#                                                                                           #
#############################################################################################

# config <- prepare_test("vector_intersect")
# read_rdata(layertype = "sf", sizes = c("s", "m"), geomtypes = "poly")
# read_rdata(layertype = "sp", sizes = c("s", "m"), geomtypes = "poly")
# load("data_input/poly_2_ssf.RData")
# load("data_input/poly_2_ssp.RData")
# test_performance_grid(config)
# remove_layer_objects(layertype = "sf", sizes = c("s", "m"))
# remove_layer_objects(layertype = "sp", sizes = c("s", "m"))
# rm(poly_2_ssf)
# rm(poly_2_ssp)

# config <- prepare_test("vector_intersect_raster")
# read_rdata(layertype = "sp", sizes = c("s"), geomtypes = "poly")
# load("data_input/poly_2_ssp.RData")
# test_performance_grid(config)
# remove_layer_objects(layertype = "sp", sizes = c("s"))
# rm(poly_2_ssf)
# rm(poly_2_ssp)

