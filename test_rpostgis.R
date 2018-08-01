library(rpostgis)
library(rgdal)
library(raster)
library(postGIStools)


# read login data 
# postgis_account.csv is in gitignore! needs to be created locally
pg_login <- read.csv("etc_local/postgis_account.csv", stringsAsFactors = FALSE)

# create pg connection
con <- dbConnect("PostgreSQL", 
                 dbname = "rperf",
                 host = "localhost",
                 user = pg_login[1, 1],
                 password = pg_login[1, 2])

# test if PostGIS is installed
pgPostGIS(con)

# load data from filesystem
poly <- readOGR("./data_input", "WAHLSPRNR2017OGDPolygon")
poly <- spTransform(poly, CRS('+init=epsg:31286'))
line <- readOGR("./data_input", "FLIESSGEWOGDLine")
line <- spTransform(line, CRS('+init=epsg:31286'))
poly_test1 <- readOGR("./data_input", "poly_test_1")
poly_test2 <- readOGR("./data_input", "poly_test_2")

# insert data into pg
dbDrop(con, "line", type = "table", ifexists = TRUE)
dbDrop(con, "poly_test1", type = "table", ifexists = TRUE)
dbDrop(con, "poly_test2", type = "table", ifexists = TRUE)
pgInsert(con, name = c("public", "poly"), data.obj = poly)
pgInsert(con, name = c("public", "line"), data.obj = line, encoding = c("latin1", "UTF-8"))
pgInsert(con, name = c("public", "poly_test1"), data.obj = poly_test1)
pgInsert(con, name = c("public", "poly_test2"), data.obj = poly_test2)

# retrieve geoms from pg, intersect, and plot
plot(pgGetGeom(con, query = "SELECT p1.geom FROM poly_test1 p1;"), 
     border = "red")
plot(pgGetGeom(con, query = "SELECT p2.geom  FROM poly_test2 p2;"), 
     border = "blue", 
     add = TRUE)
plot(pgGetGeom(con, query = "SELECT ST_Intersection(p1.geom, p2.geom) as geom FROM poly_test1 p1, poly_test2 p2;"),
     border = "black",
     lwd = 3,
     add = TRUE)

# load raster
ras <- raster("./data_input/domhbf_orig.tif")
pgWriteRast(con, c("public", "domhbf"), ras, blocks = 10, overwrite = TRUE)
pgWriteRast(con, c("public", "domhbf_8bui"), ras, bit.depth = "8BUI", blocks = 1, overwrite = TRUE)


# rasterize
dbSendQuery(con, "DROP TABLE IF EXISTS rasterized")
dbSendQuery(con, "CREATE TABLE rasterized AS 
                   SELECT 1 AS rid, ST_AsRaster((SELECT ST_Union(geom) FROM line), 100.0, 100.0) AS rast;")

plot(pgGetRast(con, c("public", "rasterized")))

plot(dbSendQuery(con, "SELECT rast FROM rasterized"))
# disconnect pg
dbDisconnect(con)
