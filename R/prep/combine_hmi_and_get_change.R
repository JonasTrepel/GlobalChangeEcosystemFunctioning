library(terra)

## change -----
r_2000 <- rast("data/spatial_data/covariates/HMv20240801_2000c_AA.tif")
r_2020 <- rast("data/spatial_data/covariates/HMv20240801_2020c_AA.tif")

r_hmi_change <- (r_2020 - r_2000)
plot(r_hmi_change)
writeRaster(r_hmi_change,"data/spatial_data/covariates/hmi_change.tif", overwrite = T)




r_2000 <- rast("data/spatial_data/covariates/HMv20240801_2000c_AA.tif")
r_2020 <- rast("data/spatial_data/covariates/HMv20240801_2020c_AA.tif")
r_change <- rast("data/spatial_data/covariates/hmi_change.tif")
plot(r_2000)
plot(r_2015)
plot(r_change)

#hist(values(r_change))


### Deprecated because now we have an updated dataset: https://zenodo.org/records/16907328
# #2000----
# (files_2000 <- list.files("data/spatial_data/hmi_tiles/gHMv1_300m_2000_change/", full.names = T))
# 
# rl_2000 <- lapply(files_2000, rast)
# 
# dtype_2000 <- datatype(rl_2000[[1]])
# 
# r_2000_sc <- merge(sprc(rl_2000), 
#                 filename = "data/spatial_data/covariates/hmi_2000_scaled.tif", 
#                 overwrite = T, 
#                 datatype = dtype_2000)
# 
# plot(r_2000_sc)
# 
# r_2000 <- round(r_2000_sc/1000, 2)
# plot(r_2000)
# writeRaster(r_2000,"data/spatial_data/covariates/hmi_2000.tif", overwrite = T)
# 
# #2015 -----
# (files_2015 <- list.files("data/spatial_data/hmi_tiles/gHMv1_300m_2015_change/", full.names = T))
# 
# rl_2015 <- lapply(files_2015, rast)
# 
# dtype_2015 <- datatype(rl_2015[[1]])
# 
# r_2015_sc <- merge(sprc(rl_2015), 
#                    filename = "data/spatial_data/covariates/hmi_2015_scaled.tif", 
#                    overwrite = T, 
#                    datatype = dtype_2015)
# 
# plot(r_2015_sc)
# 
# r_2015 <- round(r_2015_sc/1000, 2)
# plot(r_2015)
# writeRaster(r_2015,"data/spatial_data/covariates/hmi_2015.tif", overwrite = T)
# 
# ## change -----
# r_hmi_change <- (r_2015 - r_2000)
# plot(r_hmi_change)
# writeRaster(r_hmi_change,"data/spatial_data/covariates/hmi_change.tif", overwrite = T)
# 
# 
# 
# 
# r_2015 = rast("data/spatial_data/covariates/hmi_2015.tif")
# r_2000 = rast("data/spatial_data/covariates/hmi_2000.tif")
# r_change = rast("data/spatial_data/covariates/hmi_change.tif")
# plot(r_2000)
# plot(r_2015)
# plot(r_change)
# 
# hist(values(r_change))
