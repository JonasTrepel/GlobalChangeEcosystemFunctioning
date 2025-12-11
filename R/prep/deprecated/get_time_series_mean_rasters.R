# get mean time_series values 

library(terra)

#mean MAT

mat_files <- list.files("data/raw_data/time_series/", pattern = "evi", full.names = T)

mat_stack <- rast(mat_files)

mean_mat <- mean(mat_stack, na.rm = T)
plot(mean_mat)
writeRaster(mean_mat, "data/spatial_data/covariates/evi.tif", overwrite = T)

#mean Mx tempt

max_temp_files <- list.files("data/raw_data/time_series/", pattern = "era5_max_temp", full.names = T)

max_temp_stack <- rast(max_temp_files)

mean_max_temp <- mean(max_temp_stack, na.rm = T)
plot(mean_max_temp)
writeRaster(mean_max_temp, "data/spatial_data/covariates/era5_max_temp.tif", overwrite = T)

#mean Prec 

prec_files <- list.files("data/raw_data/time_series/", pattern = "era5_prec", full.names = T)

prec_stack <- rast(prec_files)

mean_prec <- mean(prec_stack, na.rm = T)
plot(mean_prec)
writeRaster(mean_prec, "data/spatial_data/covariates/era5_prec.tif", overwrite = T)

#mean burned area 

burned_area_files <- list.files("data/raw_data/time_series/", pattern = "burned_area", full.names = T)

burned_area_stack <- rast(burned_area_files)

mean_burned_area <- mean(burned_area_stack, na.rm = T)
plot(mean_burned_area)
writeRaster(mean_burned_area, "data/spatial_data/covariates/modis_burned_area.tif", overwrite = T)

#mean EVI 

evi_files <- list.files("data/raw_data/time_series/", pattern = "evi", full.names = T)

evi_stack <- rast(evi_files)

mean_evi <- mean(evi_stack, na.rm = T)
plot(mean_evi)
writeRaster(mean_evi, "data/spatial_data/covariates/modis_evi.tif", overwrite = T)