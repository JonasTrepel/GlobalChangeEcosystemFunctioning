### rasterize Protected Areas 
library(sf)
library(data.table)
library(tidyverse)
library(terra)
library(exactextractr)

#Rasterize Biomes ----------

wwf_biome <- read_sf("data/spatial_data/covariates/WWF_BIOMES.gpkg") %>% 
  st_transform(crs = "ESRI:54009")


r_tmp <- terra::rast(ext(wwf_biome), resolution = 1000, crs = "ESRI:54009") # 1 km
r_biome <- terra::rasterize(wwf_biome, r_tmp, field = "BIOME")
plot(r_biome)


writeRaster(r_biome, "data/spatial_data/covariates/wwf_olson_biome.tif", overwrite = T)


#Rasterize Pas ----------


pas <- read_sf("data/spatial_data/protected_areas/pa_shapes.gpkg")
table(pas$iucn_cat_ord)
st_crs(pas)
pas_t <- st_transform(pas, crs = "ESRI:54009")

temp_r <- rast(ext(pas_t), resolution = 250, crs = "ESRI:54009")

pas_r <- rasterize(pas_t, temp_r, field = "iucn_cat_ord", fun = "max")

#pas_r_wsg <- project(pas_r, "EPSG:4326")
plot(pas_r)

writeRaster(pas_r, "data/spatial_data/covariates/pa_raster.tif")

#pas_r <- rast("data/spatial_data/protectedAreas/pa_raster.tif")

pas_bi <- pas_r > 0
pas_bi <- as.numeric(pas_bi)
pas_bi <- ifel(is.na(pas_bi), 0, pas_bi)

plot(pas_bi)
writeRaster(pas_bi, "data/spatial_data/covariates/pa_raster_binary.tif")

### get raster of pa maps 
# 
# temp_r_id <- rast(ext(pas_t), resolution = 250, crs = "ESRI:54009")
# 
# 
# pa_ids_r <- exactextractr::rasterize_polygons(pas_t,
#                       temp_r_id)
# 
# pa_ids_r[] <- pas_t$own_pa_id[pa_ids_r[]]
# 
# plot(pa_ids_r)
# 
# writeRaster(pa_ids_r, "data/spatial_data/protected_areas/pa_id_raster.tif")
