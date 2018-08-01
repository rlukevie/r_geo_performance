library(rgdal)
library(dplyr)

poly <- readOGR("./data_input/poly_2_s.shp", "poly_2_s")
point <- readOGR("./data_input/point_l.shp", "point_l")



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
spplot(poly, "point_count")

print(max(as.numeric(as.character(poly@data$point_count))))
max(poly@data$point_count, na.rm = TRUE)
View(poly@data)
