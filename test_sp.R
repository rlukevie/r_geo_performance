library(sp)
library(rgdal)
library(rgeos)

# read Data
point_s <- readOGR("./data_input", "point_s")
poly_m <- readOGR("./data_input", "poly_m")
line_s <- readOGR("./data_input", "line_s")
poly_s <- readOGR("./data_input", "poly_s")
poly_2_s <- readOGR("./data_input", "poly_2_s")
print(head(poly_s@data))
writeOGR(poly_s, "./data_output", "poly_s.shp", driver="ESRI Shapefile")
poly <- spTransform(poly, CRS('+init=epsg:31286'))

poly_test1 <- readOGR("./data_input", "poly_test_1")

poly_test2 <- readOGR("./data_input", "poly_test_2")

# gIntersection
plot(poly_test1, border = "red")
plot(poly_test2, border = "blue", add = TRUE)
plot(gIntersection(poly_test1, poly_test2), border = "black", lwd = 3, add = TRUE)

as.numeric(as.character(point_s@data$code)) * 10

# projection
point_s_4326 <- spTransform(point_s, CRS("+init=epsg:4326"))
point_s_4326
point_s


system.time(poly_m_4326 <- spTransform(poly_m, CRS("+init=epsg:4326")))

overlay <- over(x = point_s, y = geometry(poly_2_s))
overlay[!is.na('OBJECTID')]
head(overlay)
View(overlay)
plot(poly_2_s)
plot(point_s, add = TRUE)

