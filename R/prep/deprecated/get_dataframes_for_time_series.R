### transform all the raster data to a dataframe 

library(terra)
library(tidyverse)
library(data.table)
library(exactextractr)
library(sf)

### MAT ----------------
mat_files <- data.table(filepath = list.files("data/rawData/raw_time_series/era5/era5_mat/",
                        pattern = ".tif", 
                        full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/era5/era5_mat/",
                                                          pattern = ".tif", 
                                                          full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename))
  

for(i in 1:nrow(mat_files)){
  
  mat_r <- rast(mat_files[i,]$filepath)
  plot(mat_r, main = gsub("era_5_", "", mat_files[i,]$filename))
  
  mat_tmp_df <- as.data.frame(mat_r, xy = T)
  
  names(mat_tmp_df) <- c("x", "y", gsub("era_5_", "", mat_files[i,]$filename))
  
  if(i == 1){
    mat_df <- mat_tmp_df
  }else{
    mat_df <- left_join(mat_df, mat_tmp_df)
  }

  print(paste0(gsub("era_5_", "", mat_files[i,]$filename), " done"))
}

fwrite(mat_df, "data/rawData/raw_time_series/data_frames/era5_mat_df.csv")

### Max Temp -----------------

max_temp_files <- data.table(filepath = list.files("data/rawData/raw_time_series/era5/era5_max_temp/",
                                              pattern = ".tif", 
                                              full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/era5/era5_max_temp/",
                                              pattern = ".tif", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename))


for(i in 1:nrow(max_temp_files)){
  
  max_temp_r <- rast(max_temp_files[i,]$filepath)
  plot(max_temp_r, main = gsub("era_5_", "", max_temp_files[i,]$filename))
  
  max_temp_tmp_df <- as.data.frame(max_temp_r, xy = T)
  
  names(max_temp_tmp_df) <- c("x", "y", gsub("era_5_", "", max_temp_files[i,]$filename))
  
  if(i == 1){
    max_temp_df <- max_temp_tmp_df
  }else{
    max_temp_df <- left_join(max_temp_df, max_temp_tmp_df)
  }
  
  print(paste0(gsub("era_5_", "", max_temp_files[i,]$filename), " done"))
}

fwrite(max_temp_df, "data/rawData/raw_time_series/data_frames/era5_max_temp_df.csv")


### Prec --------------------

map_files <- data.table(filepath = list.files("data/rawData/raw_time_series/era5/era5_map/",
                                                   pattern = ".tif", 
                                                   full.names = TRUE), 
                             filename = list.files("data/rawData/raw_time_series/era5/era5_map/",
                                                   pattern = ".tif", 
                                                   full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename))


for(i in 1:nrow(map_files)){
  
  map_r <- rast(map_files[i,]$filepath)
  plot(map_r, main = gsub("era_5_", "", map_files[i,]$filename))
  
  map_tmp_df <- as.data.frame(map_r, xy = T)
  
  names(map_tmp_df) <- c("x", "y", gsub("era_5_", "", map_files[i,]$filename))
  
  if(i == 1){
    map_df <- map_tmp_df
  }else{
    map_df <- left_join(map_df, map_tmp_df)
  }
  
  print(paste0(gsub("era_5_", "", map_files[i,]$filename), " done"))
}

fwrite(map_df, "data/rawData/raw_time_series/data_frames/era5_map_df.csv")


### EVI ---------------------

#get evi files 
evi_files <- data.table(filepath = list.files("data/rawData/raw_time_series/evi/",
                                                   pattern = ".tif", 
                                                   full.names = TRUE), 
                             filename = list.files("data/rawData/raw_time_series/evi/",
                                                   pattern = ".tif", 
                                                   full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename))

#create landcover mask 
# lc_r <- rast("data/spatialData/otherCovariates/GlobalLandCoverCopernicus2019.tif")

# evi_temp_r <- rast(evi_files[12,]$filepath) #load evi raster as template

# lc_res <- exact_resample(lc_r, evi_temp_r, fun = "mode")

# lc_mask <- (lc_res %in% c(40, 50, 80, 200))*1

# plot(lc_mask)

#writeRaster(lc_mask, "data/spatialData/masks/land_cover_mask_modis_evi.tif")

lc_mask <- rast("data/spatialData/masks/land_cover_mask_modis_evi.tif")

roi <- read_sf("data/spatialData/regions_of_interest.gpkg")
roi_v <- vect(roi)

for(i in 1:nrow(evi_files)){
  

  evi_r <- rast(evi_files[i,]$filepath)
  
  crs(evi_r) == crs(roi_v)
  
  evi_masked <- mask(evi_r, roi_v, updatevalue = NA)
  
 # evi_masked <- mask(evi_r, lc_mask, maskvalue = 1, updatevalue = NA)
  
  plot(evi_masked, main = gsub("evi_modis_median_500m", "evi", evi_files[i,]$filename))
  
  tilename <- paste0("data/rawData/raw_time_series/tmp_tiles/",evi_files[i,]$filename, "_tile_.tif")
  
  evi_tiles <- makeTiles(x = evi_masked, y = 25000,
                         filename= tilename, 
                         overwrite = TRUE) 
  
  evi_tmp_df <- data.table()
  
  for(j in 1:n_distinct(evi_tiles)){
    
    evi_t <- rast(evi_tiles[j])
    plot(evi_t)
    
    evi_tile_df <- as.data.frame(evi_t, xy = T)
    
    evi_tmp_df <- rbind(evi_tile_df, evi_tmp_df) %>% as.data.table()
    
    print(paste(j, "done"))
    
  }

  names(evi_tmp_df) <- c("x", "y", gsub("evi_modis_median_500m", "evi", evi_files[i,]$filename))
  
  if(i == 1){
    evi_df <- evi_tmp_df
  }else{
    evi_df <- left_join(evi_df, evi_tmp_df)
  }
  
  setDT(evi_df)
  print(paste0(gsub("evi_modis_median_500m", "evi", evi_files[i,]$filename), " done"))
}

fwrite(evi_df, "data/rawData/raw_time_series/data_frames/modis_evi_df.csv")


### Burned area -------------------

#get burned_area files 
burned_area_files <- data.table(filepath = list.files("data/rawData/raw_time_series/fire/burned_area/",
                                              pattern = ".tif", 
                                              full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/fire/burned_area/",
                                              pattern = ".tif", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename))

#create landcover mask 
# lc_r <- rast("data/spatialData/otherCovariates/GlobalLandCoverCopernicus2019.tif")

# burned_area_temp_r <- rast(burned_area_files[12,]$filepath) #load burned_area raster as template

# lc_res <- exact_resample(lc_r, burned_area_temp_r, fun = "mode")

# lc_mask <- (lc_res %in% c(40, 50, 80, 200))*1

# plot(lc_mask)

#writeRaster(lc_mask, "data/spatialData/masks/land_cover_mask_modis_burned_area.tif")

roi <- read_sf("data/spatialData/regions_of_interest.gpkg")
roi_v <- vect(roi)

for(i in 1:nrow(burned_area_files)){
  
  
  burned_area_r <- rast(burned_area_files[i,]$filepath)
  
  crs(burned_area_r) == crs(roi_v)
  
  burned_area_masked <- mask(burned_area_r, roi_v, updatevalue = NA)
  
  plot(burned_area_masked, main = gsub("burned_area_modis_median_500m", "burned_area", burned_area_files[i,]$filename))
  
  burned_area_tmp_df <- as.data.frame(burned_area_masked, xy = T)
  
  names(burned_area_tmp_df) <- c("x", "y", gsub("burned_area_modis_5000m", "burned_area", burned_area_files[i,]$filename))
  
  if(i == 1){
    burned_area_df <- burned_area_tmp_df
  }else{
    burned_area_df <- left_join(burned_area_df, burned_area_tmp_df)
  }
  
  setDT(burned_area_df)
  print(paste0(gsub("burned_area_modis_5000m", "burned_area", burned_area_files[i,]$filename), " done"))
}

fwrite(burned_area_df, "data/rawData/raw_time_series/data_frames/modis_burned_area_df.csv")


### Start of growing season ---------------------


#get sos files 
sos_files <- data.table(filepath = list.files("data/rawData/raw_time_series/phenology/",
                                              pattern = "doy_greenup_1", 
                                              full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/phenology/",
                                              pattern = "doy_greenup_1", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename))

#create landcover mask 
# lc_r <- rast("data/spatialData/otherCovariates/GlobalLandCoverCopernicus2019.tif")

# sos_temp_r <- rast(sos_files[12,]$filepath) #load sos raster as template

# lc_res <- exact_resample(lc_r, sos_temp_r, fun = "mode")

# lc_mask <- (lc_res %in% c(20, 30, 60, 70, 90, 100,
#                           111, 112, 113, 114, 115, 116, 
#                           121, 122, 123, 124, 125, 126))*1

# plot(lc_mask)

# writeRaster(lc_mask, "data/spatialData/masks/land_cover_mask_modis_sos.tif")

lc_mask_sos <- rast("data/spatialData/masks/land_cover_mask_modis_sos.tif")

roi <- read_sf("data/spatialData/regions_of_interest.gpkg")
roi_v <- vect(roi)


for(i in 1:nrow(sos_files)){
  
  sos_r <- rast(sos_files[i,]$filepath)
  sos_masked <- mask(sos_r, roi_v, updatevalue = NA)
 # sos_masked <- mask(sos_r, lc_mask_sos, maskvalue = 0, updatevalue = NA)
 
  plot(sos_masked, main = gsub("sos_modis_median_500m", "sos", sos_files[i,]$filename))
  
  tilename <- paste0("data/rawData/raw_time_series/tmp_tiles/",sos_files[i,]$filename, "_tile_.tif")
  
  sos_tiles <- makeTiles(x = sos_masked, y = 25000,
                         filename = tilename, 
                         overwrite = TRUE) 
  
  sos_tmp_df <- data.table()
  
  for(j in 1:n_distinct(sos_tiles)){
    
    sos_t <- rast(sos_tiles[j])
    plot(sos_t)
    
    sos_tile_df <- as.data.frame(sos_t, xy = T)
    
    sos_tmp_df <- rbind(sos_tile_df, sos_tmp_df) %>% as.data.table()
    
    print(paste(j, "done"))
    
  }
  
  names(sos_tmp_df) <- c("x", "y", gsub("doy_greenup_1_modis_500m", "greenup_1", sos_files[i,]$filename))
  
  if(i == 1){ 
    sos_df <- sos_tmp_df 
  }else{
    sos_df <- left_join(sos_df, sos_tmp_df) 
    }
  
  setDT(sos_tmp_df)
  
  fwrite(sos_df, "data/rawData/raw_time_series/data_frames/modis_greenup_1_df.csv")
  
  print(paste0(gsub("sos_modis_median_500m", "sos", sos_files[i,]$filename), " done"))
}


fwrite(sos_df, "data/rawData/raw_time_series/data_frames/modis_greenup_1_df.csv")
