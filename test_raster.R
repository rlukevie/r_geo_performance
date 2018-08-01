library(raster)
r_s = raster("data_input/raster_s.tif")
rcl = matrix(c(-200, 0, 1, 0, 100, 2, 100, 300, 3), ncol = 3, byrow = TRUE)
rcl
system.time(r_s_rcl <- reclassify(r_s, rcl = rcl))[["elapsed"]]
plot(r_s_rcl)
r_s_rcl
r_s
