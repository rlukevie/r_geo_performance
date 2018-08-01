library(sf)


poly_2_s <- st_read("data_input/poly_2_s.shp")
poly_m <- st_read("data_input/poly_m.shp")
point_s <- st_read("data_input/point_s.shp")
point_l <- st_read("data_input/point_l.shp")




# st_join
joined = st_join(x = point_l, y = poly_m)
plot(joined)
joined_agg <- joined %>% 
  group_by(OBJECTID) %>%
  summarise(point_count = n())
max(joined_agg$point_count, na.rm = TRUE)
max(subset(joined_agg, !is.na(joined_agg$OBJECTID))$point_count)

subset(joined_agg, !is.na(joined_agg$point_count)) %>% arrange(desc(point_count))

View(joined_agg)
point_l
joined

joined %>% aggregate(list(.$joined.OBJECTID), count)

st_crs(poly_m)
st_crs(poly_2_s)
poly_m <- st_transform(poly_m, 31256)
st_write(obj = poly_m, dsn = "data_output/poly_m.shp", layer = "poly_m", delete_layer = TRUE)

# Convex Hull

# st_triangulate
plot(st_triangulate(st_union(point_s_geom)), col = 0)
plot(point_s_geom, add = TRUE)

# Voronoi
plot(st_voronoi(st_union(point_s)), col = 0)

# st_polygonize
plot(line_geom)
plot(st_polygonize(st_union(line_geom)))

# st_line_merge
st_line_merge(line_geom)
line_geom

# st_point_on_surface
plot(poly_s_geom)
plot(st_point_on_surface(poly_s_geom), add = TRUE)
plot(st_centroid(poly_s_geom), col = 4, add = TRUE)

# st_node
st_node(line_geom)
line_geom

# st_segmentize
plot(st_segmentize(poly_s_geom, dfMaxLength = 100))

# st_agr
st_set_agr(line_geom)

# st_as_binary
st_as_binary(point_s_geom)

# st_as_grob
# ???

# st_as_text
st_as_text(point_s_geom)

# st_coordinates
st_coordinates(point_s)

# st_interpolate_aw (example makes no sense, only technical test)
g = st_make_grid(poly, n = c(20,20))
plot(st_interpolate_aw(poly["BEZIRK"], g, extensive = TRUE)["BEZIRK"])


# st_intersection
plot(poly_test1_geom, border = "red")
plot(poly_test2_geom, border = "blue", add = TRUE)
plot(st_intersection(poly_test1_geom, poly_test2_geom), border = "black", add = TRUE, lwd = 3)
