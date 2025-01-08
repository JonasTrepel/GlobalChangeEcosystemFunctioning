### rasterize Protected Areas 
library(sf)
library(data.table)
library(tidyverse)
library(terra)
library(exactextractr)

pas <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg")
table(pas$iucnCatOrd)

be_crs <- st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

pas_t <- st_transform(pas, crs = be_crs)

temp_r <- rast(ext(pas_t), resolution = 500, crs = be_crs$wkt)

pas_r <- rasterize(pas_t, temp_r, field = "iucnCatOrd", fun = "max")

#pas_r_wsg <- project(pas_r, "EPSG:4326")
plot(pas_r)

writeRaster(pas_r, "data/spatialData/protectedAreas/pa_raster.tif")

#pas_r <- rast("data/spatialData/protectedAreas/pa_raster.tif")

pas_bi <- pas_r > 0
pas_bi <- as.numeric(pas_bi)
pas_bi <- ifel(is.na(pas_bi), 0, pas_bi)

plot(pas_bi)
writeRaster(pas_bi, "data/spatialData/protectedAreas/pa_raster_binary.tif")

### get raster of pa maps 

temp_r_id <- rast(ext(pas_t), resolution = 500, crs = be_crs$wkt)


Mode <- function(x, na.rm = TRUE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

pa_ids_r <- rasterize_polygons(pas_t,
                      temp_r_id)

pa_ids_r[] <- pas_t$own_pa_id[pa_ids_r[]]

plot(pa_ids_r)

writeRaster(pa_ids_r, "data/spatialData/protectedAreas/pa_id_raster.tif")
