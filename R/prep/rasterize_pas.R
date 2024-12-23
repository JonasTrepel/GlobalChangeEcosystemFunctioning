### rasterize Protected Areas 
library(sf)
library(data.table)
library(tidyverse)
library(terra)

pas <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg")
table(pas$iucnCatOrd)

pas_trans <- st_transform(pas, crs = "epsg:3857")
temp_r <- rast(ext(pas_trans), resolution = 500, crs = "epsg:3857")

pas_r <- rasterize(pas_trans, temp_r, field = "iucnCatOrd", fun = "max")

pas_r_wsg <- project(pas_r, "EPSG:4326")
plot(pas_r_wsg)

writeRaster(pas_r_wsg, )

