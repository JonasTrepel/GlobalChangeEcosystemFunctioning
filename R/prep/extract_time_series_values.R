## extract time series for vects 

library(sf)
library(terra)
library(exactextractr)
library(data.table)
library(tidyverse)

#read vectors 
#param = "grid"
#param = "pas"
#param = "usa"
param = "europe"

if(param == "grid"){
  vect <- st_read("data/spatialData/grid_sample.gpkg") 
}else if(param == "pas"){
  vect <- st_read("data/spatialData/protectedAreas/pa_and_control_grid_1km_with_covs.gpkg")
} else if(param == "usa"){
  vect <- st_read("data/spatialData/grid_usa.gpkg")
} else if(param == "europe"){
  vect <- st_read("data/spatialData/grid_europe.gpkg")
}


## get file paths sorted 

mat_files <- data.table(filepath = list.files("data/rawData/raw_time_series/era5/era5_mat/",
                                              pattern = ".tif", 
                                              full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/era5/era5_mat/",
                                              pattern = ".tif", 
                                              full.names = FALSE)
                        ) %>% 
  mutate(filename = gsub(".tif", "", filename), 
         colname = gsub("era_5_", "", filename))

max_temp_files <- data.table(filepath = list.files("data/rawData/raw_time_series/era5/era5_max_temp/",
                                                   pattern = ".tif", 
                                                   full.names = TRUE), 
                             filename = list.files("data/rawData/raw_time_series/era5/era5_max_temp/",
                                                   pattern = ".tif", 
                                                   full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename), 
         colname = gsub("era_5_", "", filename))

map_files <- data.table(filepath = list.files("data/rawData/raw_time_series/era5/era5_map/",
                                              pattern = ".tif", 
                                              full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/era5/era5_map/",
                                              pattern = ".tif", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename), 
         colname = gsub("era_5_", "", filename))


evi_files <- data.table(filepath = list.files("data/rawData/raw_time_series/evi/",
                                              pattern = ".tif", 
                                              full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/evi/",
                                              pattern = ".tif", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename), 
         colname = gsub("evi_modis_median_500m", "median_evi", filename))

burned_area_files <- data.table(filepath = list.files("data/rawData/raw_time_series/fire/burned_area/",
                                                      pattern = ".tif", 
                                                      full.names = TRUE), 
                                filename = list.files("data/rawData/raw_time_series/fire/burned_area/",
                                                      pattern = ".tif", 
                                                      full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         colname =  gsub("burned_area_modis_median_500m", "burned_area", filename))


sos_files <- data.table(filepath = list.files("data/rawData/raw_time_series/phenology/",
                                              pattern = "doy_greenup_1", 
                                              full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/phenology/",
                                              pattern = "doy_greenup_1", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         colname =  gsub("_modis_500m", "", filename))


mean_evi_files <- data.table(filepath = list.files("data/rawData/raw_time_series/mean_evi/",
                                              pattern = "envi_", 
                                              full.names = TRUE), 
                        filename = list.files("data/rawData/raw_time_series/mean_evi/",
                                              pattern = "envi_", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         colname =  gsub("envi_mean_500m", "mean_evi", filename))


n_depo_zhu_files <- data.table(filepath = list.files("data/rawData/raw_time_series/zhu_2025_nitrogen_depo/Global_N_deposition_grid_dataset_2008_2020/",
                                                   pattern = "mean_totN_", 
                                                   full.names = TRUE), 
                             filename = list.files("data/rawData/raw_time_series/zhu_2025_nitrogen_depo/Global_N_deposition_grid_dataset_2008_2020/",
                                                   pattern = "mean_totN_", 
                                                   full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         filename = gsub("_hm", "", filename),
         colname =  gsub("mean_totN", "n_depo_zhu", filename))

n_depo_usa_files <- data.table(filepath = list.files("data/rawData/raw_time_series/usa_n_depo/",
                                                     pattern = "n_tw", 
                                                     full.names = TRUE), 
                               filename = list.files("data/rawData/raw_time_series/usa_n_depo/",
                                                     pattern = "n_tw", 
                                                     full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         colname =  gsub("n_tw-", "n_depo_usa_", filename))

n_depo_europe_files <- data.table(filepath = list.files("data/rawData/raw_time_series/europe_n_depo/",
                                                     pattern = "europe_n_depo", 
                                                     full.names = TRUE), 
                               filename = list.files("data/rawData/raw_time_series/europe_n_depo/",
                                                     pattern = "europe_n_depo", 
                                                     full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         colname =  gsub("europe_n_depo_", "n_depo_europe_", filename))


covs <- rbind(mat_files, 
              max_temp_files, 
              map_files,
              evi_files,
              burned_area_files,
              sos_files, 
              mean_evi_files, 
              n_depo_zhu_files)

if(param == "usa"){
  covs <- rbind(covs, n_depo_usa_files)
}
if(param == "europe"){
  covs <- rbind(covs, n_depo_europe_files)
}

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

# Create and register a cluster
clust <- makeCluster(56)
registerDoSNOW(clust)

## progress bar 
iterations <- nrow(covs)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

##############################################################################            
################################## LOOOOOOOOOOOOP ############################            
##############################################################################    
#dt.tier <- dt.tier[3,]
tic()

dt_covs <- foreach(i = 1:nrow(covs),
                  .packages = c('tidyverse', 'exactextractr', 'data.table', 'terra', 'sf'),
                  .options.snow = opts,
                  .inorder = TRUE,
                  .verbose = TRUE,
                  .combine = left_join) %dopar% {
                    
                    #for(i in 1:nrow(covs)){
                    
                    cov_r <- rast(covs[i, ]$filepath)
                    
                    vect_trans <- st_transform(vect, crs = st_crs(cov_r))
                    
                    extr <- exactextractr::exact_extract(cov_r, 
                                                         append_cols = c("unique_id"),
                                                         vect_trans, 
                                                         fun = "mean")
     
                    setnames(extr, "mean", covs[i, ]$colname)
                    
                    dt_extr_fin <- extr %>% 
                      as.data.table() %>%
                      mutate(geom = NULL) %>% 
                      unique()
                   
                    return(dt_extr_fin)
                    
    }

toc()
print(Sys.time())
stopCluster(clust)

#combine
vect_covs <- vect %>% 
  as.data.table() %>% 
  mutate(x = NULL, geom = NULL, geometry = NULL) %>% 
  left_join(dt_covs) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL,
         geometry = NULL) %>% 
  mutate(mean_burned_area = rowMeans(select(., contains("burned_area")), na.rm = TRUE), 
         max_burned_area = apply(select(., contains("burned_area")), 1, max, na.rm = TRUE),
         map_era = rowMeans(select(., contains("map_")), na.rm = TRUE), 
         mat_era = rowMeans(select(., contains("mat_")), na.rm = TRUE), 
         max_temp_era = apply(select(., contains("max_temp")), 1, max, na.rm = TRUE), 
         mean_evi = rowMeans(select(., contains("mean_evi")), na.rm = TRUE), 
         mean_greenup = rowMeans(select(., contains("greenup")), na.rm = TRUE), 
         mean_n_depo_zhu = rowMeans(select(., contains("n_depo_zhu")), na.rm = TRUE)) 


if(param == "grid"){
  fwrite(vect_covs, "data/processedData/dataFragments/grid_sample_with_raw_timeseries.csv")
}else if(param == "pas"){
  fwrite(vect_covs, "data/processedData/dataFragments/pa_and_controls_with_raw_timeseries.csv")
} else if(param == "usa"){
  
  vect_covs <- vect_covs %>% 
    mutate(mean_n_depo_usa = rowMeans(select(., contains("n_depo_usa")), na.rm = TRUE))
  
  fwrite(vect_covs, "data/processedData/dataFragments/grid_usa_with_raw_timeseries.csv")
} else if(param == "europe"){
  
  vect_covs <- vect_covs %>% 
    mutate(mean_n_depo_europe = rowMeans(select(., contains("n_depo_europe")), na.rm = TRUE))
  
  fwrite(vect_covs, "data/processedData/dataFragments/grid_europe_with_raw_timeseries.csv")
}

