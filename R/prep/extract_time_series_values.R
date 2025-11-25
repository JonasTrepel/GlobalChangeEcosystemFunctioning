## extract time series for vects 

library(sf)
library(terra)
library(exactextractr)
library(data.table)
library(tidyverse)

#read vectors 
#param = "grid"
#param = "pas"
param = "pa_controls"
#param = "usa"
#param = "europe"

if(param == "grid"){
  
  vect <- st_read("data/spatial_data/grids/world_grid.gpkg")
  
}else if(param == "pas"){
  
  vect <- st_read("data/spatial_data/protected_areas/pa_shapes.gpkg")
  
  coords <- vect %>% st_transform(crs = "ESRI:54009") %>%
    st_centroid() %>% st_coordinates()
  
  sf_use_s2(FALSE)
  coords_lat_lon <- vect %>% st_centroid() %>% st_coordinates()
  sf_use_s2(TRUE)
  
  
  vect$x_moll <- coords[,1]
  vect$y_moll <- coords[,2]
  
  vect$lon <- coords_lat_lon[,1]
  vect$lat <- coords_lat_lon[,2]
  
}  else if(param == "pa_controls"){
  
  vect <- read_sf("data/spatial_data/protected_areas/controls_for_pas.gpkg")
  
  coords <- vect %>% st_transform(crs = "ESRI:54009") %>% st_centroid() %>% st_coordinates()
 
  sf_use_s2(FALSE)
  coords_lat_lon <- vect %>% st_centroid() %>% st_coordinates()
  sf_use_s2(TRUE)
  
  
  vect$x_moll <- coords[,1]
  vect$y_moll <- coords[,2]
  
  vect$lon <- coords_lat_lon[,1]
  vect$lat <- coords_lat_lon[,2]
  
} else if(param == "usa"){
  
  vect_raw <- st_read("data/spatial_data/grids/world_grid.gpkg")
  
  #outline from https://public.opendatasoft.com/explore/dataset/us-state-boundaries/export/
  usa_all <- st_read("data/spatial_data/random_shapefiles//usa_boundaries/us-state-boundaries.shp")
  
  usa <- usa_all %>% 
    filter(name %in% c(
      "Alabama", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", 
      "Delaware", "Florida", "Georgia", "Idaho", "Illinois", "Indiana", "Iowa", 
      "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", 
      "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", 
      "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York", 
      "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", 
      "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", 
      "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", 
      "West Virginia", "Wisconsin", "Wyoming"
    )) %>% summarize() %>% st_transform(crs = st_crs(vect_raw))
  
  vect <- vect_raw %>% 
    filter(lengths(st_intersects(., usa)) > 0)
  
  #mapview::mapview(grid_usa)
  
  write_sf(vect, "data/spatial_data/grids/grid_usa.gpkg", append = FALSE)
  
  
  
} else if(param == "europe"){
  
  vect_raw <- st_read("data/spatial_data/grids/world_grid.gpkg")
  
  r <- rast("data/raw_data/time_series/europe_n_depo_2020.tif")
  
  extent_shp <- as.polygons(ext(r), crs = crs(r)) %>% 
    st_as_sf() %>% 
    st_transform(crs = st_crs(vect_raw))
  
  set.seed(161)
  vect <- vect_raw %>% 
    filter(lengths(st_intersects(., extent_shp)) > 0)
  
  write_sf(vect, "data/spatial_data/grids/grid_europe.gpkg", append = FALSE)
  
}


## get file paths sorted 

(mat_files <- data.table(filepath = list.files("data/raw_data/time_series/",
                                              pattern = "era5_mat", 
                                              full.names = TRUE), 
                        filename = list.files("data/raw_data/time_series/",
                                              pattern = "era5_mat", 
                                              full.names = FALSE)
                        ) %>% 
  mutate(filename = gsub(".tif", "", filename), 
         colname = gsub("era5_", "", filename)))

(max_temp_files <- data.table(filepath = list.files("data/raw_data/time_series/",
                                                   pattern = "era5_max_temp", 
                                                   full.names = TRUE), 
                             filename = list.files("data/raw_data/time_series/",
                                                   pattern = "era5_max_temp", 
                                                   full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename), 
         colname = gsub("era5_", "", filename)))

(map_files <- data.table(filepath = list.files("data/raw_data/time_series/",
                                              pattern = "prec", 
                                              full.names = TRUE), 
                        filename = list.files("data/raw_data/time_series/",
                                              pattern = "prec", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename), 
         colname = gsub("era5_", "", filename)))


(evi_files <- data.table(filepath = list.files("data/raw_data/time_series/",
                                              pattern = "evi", 
                                              full.names = TRUE), 
                        filename = list.files("data/raw_data/time_series/",
                                              pattern = "evi", 
                                              full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename), 
         colname = gsub("modis_", "", filename),
         colname = gsub("500m_", "", colname)))

(burned_area_files <- data.table(filepath = list.files("data/raw_data/time_series/",
                                                      pattern = "burned_area", 
                                                      full.names = TRUE), 
                                filename = list.files("data/raw_data/time_series/",
                                                      pattern = "burned_area", 
                                                      full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         colname =  gsub("modis_5000m_", "", filename)))




(n_depo_zhu_files <- data.table(filepath = list.files("data/raw_data/time_series/",
                                                   pattern = "mean_totN_", 
                                                   full.names = TRUE), 
                             filename = list.files("data/raw_data/time_series/",
                                                   pattern = "mean_totN_", 
                                                   full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         filename = gsub("_hm", "", filename),
         colname =  gsub("mean_totN", "n_depo_zhu", filename)))

(n_depo_usa_files <- data.table(filepath = list.files("data/raw_data/time_series/",
                                                     pattern = "usa_n_depo", 
                                                     full.names = TRUE), 
                               filename = list.files("data/raw_data/time_series/",
                                                     pattern = "usa_n_depo", 
                                                     full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         colname =  gsub("n_tw-", "n_depo_usa_", filename)))

(n_depo_europe_files <- data.table(filepath = list.files("data/raw_data/time_series/",
                                                     pattern = "europe_n_depo", 
                                                     full.names = TRUE), 
                               filename = list.files("data/raw_data/time_series/",
                                                     pattern = "europe_n_depo", 
                                                     full.names = FALSE)) %>% 
  mutate(filename = gsub(".tif", "", filename),
         colname =  gsub("europe_n_depo_", "n_depo_europe_", filename)))


covs <- rbind(mat_files, 
              max_temp_files, 
              map_files,
              evi_files,
              burned_area_files,
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
clust <- makeCluster(40)
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
                    
                    if(param == "pa_controls"){
                      vect_trans <- vect_trans %>% st_make_valid()
                      vect_trans <- vect_trans[!st_is_empty(vect_trans) & !is.na(st_dimension(vect_trans)), ]
                    }
                    
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
                    
                    rm(cov_r)
                    rm(vect_trans) 
                    gc()
                    
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
  mutate(mean_burned_area = rowMeans(select(., contains("burned_area_")), na.rm = TRUE), 
         mean_prec = rowMeans(select(., contains("prec_")), na.rm = TRUE), 
         mean_mat = rowMeans(select(., contains("mat_")), na.rm = TRUE), 
         mean_max_temp = rowMeans(select(., contains("max_temp_")), na.rm = TRUE), 
         mean_evi = rowMeans(select(., contains("evi_")), na.rm = TRUE), 
         mean_n_depo_zhu = rowMeans(select(., contains("n_depo_zhu_")), na.rm = TRUE)) 


if(param == "grid"){
  
  fwrite(vect_covs, "data/processed_data/grid_with_timeseries.csv")
  
}else if(param == "pas"){
  
  fwrite(vect_covs, "data/processed_data/pas_with_timeseries.csv")
  
} else if(param == "pa_controls"){
  
  fwrite(vect_covs, "data/processed_data/pa_controls_with_timeseries.csv")
  
} else if(param == "usa"){
  
  vect_covs <- vect_covs %>% 
    mutate(mean_n_depo_usa = rowMeans(select(., contains("usa_n_depo_")), na.rm = TRUE))
  
  fwrite(vect_covs, "data/processed_data/grid_usa_with_timeseries.csv")
} else if(param == "europe"){
  
  vect_covs <- vect_covs %>% 
    mutate(mean_n_depo_europe = rowMeans(select(., contains("n_depo_europe_")), na.rm = TRUE))
  
  fwrite(vect_covs, "data/processed_data/grid_europe_with_timeseries.csv")
}

