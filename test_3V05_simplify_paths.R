#########################
# 3V05: SIMPLIFY PATHS  #
#########################

# sp              --> NO
# sf              --> YES (0.4 sec, uses rgeos)
# rgeos           --> YES (0.9 sec)
# maptools        --> YES (0.9 sec, uses rgeos)
# geojsonio       --> NO
# rmapshaper      --> YES (4.6 sec)

# raster          --> NO
# spatial.tools   --> NO
# gdalUtils       --> NO
# landscapetools  --> NO
# velox           --> NO

# rgdal           --> NO
# RQGIS           --> YES (uses GEOS)

# fasterize       --> NO
# ----------------------------------

# --> either the same geometry package (rgeos) is used by different packages,
#     or different algorithms are used
#     --> comparison is senseless,
#     --> results are not comparable, thus no performance measurement is conducted

library(spbabel)

# sf
library(sf)
load("data_input/poly_2_ssf.RData")

start_time <- Sys.time()
simplified_sf <- st_simplify(poly_2_ssf, preserveTopology = TRUE, dTolerance = 0.01)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

plot(st_geometry(simplified_sf))

format(object.size(poly_2_ssf), units = "Mb")
nrow(sptable(poly_2_ssf))

format(object.size(simplified_sf), units = "Mb")
nrow(sptable(simplified_sf))


# maptools (uses rgeos)
library(maptools)
load("data_input/poly_2_ssp.RData")

start_time <- Sys.time()
simplified_maptools <- thinnedSpatialPoly(SP = poly_2_ssp, tolerance = 0.01, topologyPreserve = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

plot(simplified_maptools)

format(object.size(poly_2_ssp), units = "Mb")
nrow(sptable(poly_2_ssp))

format(object.size(simplified_maptools), units = "Mb")
nrow(sptable(simplified_maptools))


# rgeos
library(rgeos)
load("data_input/poly_2_ssp.RData")

start_time <- Sys.time()
simplified_rgeos <- gSimplify(poly_2_ssp, tol = 0.01, topologyPreserve = TRUE)
simplified_rgeos <- SpatialPolygonsDataFrame(simplified_rgeos, data = poly_2_ssp@data)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

plot(simplified_rgeos)

format(object.size(poly_2_ssp), units = "Mb")
nrow(sptable(poly_2_ssp))

format(object.size(simplified_rgeos), units = "Mb")
nrow(sptable(simplified_rgeos))


# rmapshaper
library(rmapshaper)
load("data_input/poly_2_ssf.RData")

start_time <- Sys.time()
simplified_rmapshaper <- ms_simplify(poly_2_ssf, keep = 0.01, keep_shapes = TRUE)
print(paste0("Duration: ", difftime(Sys.time(), start_time, units = "secs"), " sec"))

plot(st_geometry(simplified_rmapshaper))

format(object.size(simplified_rmapshaper), units = "Mb")
nrow(sptable(simplified_rmapshaper))

plot(st_geometry(poly_2_ssf))
format(object.size(poly_2_ssf), units = "Mb")
nrow(sptable(poly_2_ssf))
