## extract time series for grids 

library(sf)
library(terra)
library(exactextractr)
library(data.table)
library(tidyverse)

#read vectors 

grid <- st_read("data/spatialData/grid_sample.gpkg")

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
         colname = gsub("evi_modis_median_500m", "evi", filename))

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


covs <- rbind(mat_files, 
              max_temp_files, 
              map_files,
              evi_files,
              burned_area_files,
              sos_files)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

# Create and register a cluster
clust <- makeCluster(50)
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
                  .inorder = FALSE,
                  .combine = left_join) %dopar% {
                    
                    #for(i in 1:nrow(covs)){
                    
                    cov_r <- rast(covs[i, ]$filepath)
                    
                    grid_trans <- st_transform(grid, crs = st_crs(cov_r))
                    
                    extr <- exactextractr::exact_extract(cov_r, 
                                                         grid_trans, 
                                                         fun = "mean")
                    dt_extr <- data.table(
                      extr_col = extr
                    )
                   
                    setnames(dt_extr, "extr_col", covs[i, ]$colname)
                    
                    dt_extr2 <- cbind(grid[, "gridID"], dt_extr) %>%
                      as.data.table() %>%
                      mutate(geom = NULL) %>% 
                      unique()
                    
                    dt_extr_fin <- grid[, "gridID"] %>% 
                      left_join(dt_extr2) %>% 
                      as.data.table() %>%
                      mutate(geom = NULL) %>% 
                      unique()
                   
                    return(dt_extr_fin)
                    
    }

toc()


#combine
grid_covs <- grid %>% 
  as.data.table() %>% 
  mutate(x = NULL, geom = NULL, geometry = NULL) %>% 
  left_join(dt_covs) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL,
         geometry = NULL)

fwrite(grid_covs, "data/processedData/dataFragments/grid_sample_with_raw_timeseries")

# Identify duplicated gridIDs
duplicated_ids <- grid %>%
  filter(gridID %in% grid$gridID[duplicated(grid$gridID)])

# Display rows with duplicated gridIDs
print(duplicated_ids)

# Optional: Count how many times each gridID is duplicated
duplicates_summary <- duplicated_ids %>%
  count(gridID, sort = TRUE)

print(duplicates_summary)