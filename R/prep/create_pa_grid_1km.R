
library(data.table)
library(tidyverse)
library(tidylog)
library(tictoc)
library(sf)
library(mapview)


#source("R/functions/move_polygon.R")
#source("R/functions/resolve_overlaps.R")

#make small PA grid 

shapes_all <- read_sf("data/spatialData/pa_and_control_shapes.gpkg") 

be_crs <- st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

pa_shapes <- shapes_all %>% 
  st_transform(crs = be_crs)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

nCores <- parallel::detectCores()
# Create and register a cluster
clust <- makeCluster(40)
registerDoSNOW(clust)

## progress bar 
iterations <- nrow(pa_shapes)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


pa_controls_raw <- NULL

print(Sys.time())

pa_grid <- foreach(i = 1:nrow(pa_shapes),
                           .packages = c('sf', 'ggplot2', 'data.table', 'tidyverse'),
                           .options.snow = opts,
                           .inorder = FALSE,
                           .verbose = TRUE,
                           .errorhandling = "remove",
                           .combine = rbind) %dopar% {

#for(i in 1:nrow(pa_shapes[1:10, ])){
  
  pa <- pa_shapes[i, ]
  
  tmp_grid <- st_make_grid(pa, cellsize = c(1000, 1000)) %>%
    st_as_sf() %>% 
    mutate(id = 1:nrow(.))
  
  pa <- st_make_valid(pa)
  
  int <- st_intersection(pa, tmp_grid) 
  
  int$area <- st_area(int)
  grid_area <- st_area(tmp_grid) %>% unique()
  
  int <- int %>% 
    mutate(fract_prot = as.numeric(area/grid_area))
  
  ids <- int %>% as.data.table() %>% mutate(geom = NULL) %>% filter(fract_prot == 1) %>% dplyr::select(id) %>% pull()
  
  meta_d <- int %>%
    as.data.table() %>%
    mutate(geom = NULL)
  
  tmp_grid <- tmp_grid %>%
    filter(id %in% ids)%>%
    left_join(meta_d) %>% 
    dplyr::select(-id)
  
  return(tmp_grid)
  
  #if(i == 1){
  #  pa_grid <- tmp_grid
  #} else {
  #  pa_grid <- rbind(tmp_grid, pa_grid)
  # }
  
  print(i)
}

print(paste0("done at: ", Sys.time()))
stopCluster(clust)

n_distinct(pa_grid$unique_id)

write_sf(pa_grid, "data/spatialData/protectedAreas/pa_and_control_grid_1km.gpkg")

