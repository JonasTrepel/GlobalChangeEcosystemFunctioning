### try to find spatial equivalents outside of protected areas but in the same climatic range 

library(data.table)
library(tidyverse)
library(ggridges)
library(scico)
library(mgcv)
library(MuMIn)
library(broom)
library(tictoc)
library(gridExtra)
library(sf)
library(mapview)

source("R/functions/move_polygon.R")
source("R/functions/resolve_overlaps.R")

# just gotta do two data sets. One strict, one not so strict. 
dt_pas <- fread("data/processed_data/pas_with_covariates.csv") %>% 
  mutate(STATUS_YR = ifelse(STATUS_YR == 0, NA, STATUS_YR), 
         pa_age = 2024-STATUS_YR) %>% 
  filter(!land_cover %in% c("water", "ocean", "urban", "cultivated")) %>% unique()

dt_strict_pas <- dt_pas %>% 
  filter(iucn_cat %in% c("Ia", "Ib", "II"))

dt_loose_pas <- dt_pas %>% 
  filter(!iucn_cat %in% c("Ia", "Ib", "II"))

strict_pa_shapes_raw <- read_sf("data/spatial_data/protected_areas/pa_shapes.gpkg") %>% 
  dplyr::select(unique_id) %>% 
  filter(unique_id %in% unique(dt_strict_pas$unique_id)) %>%
  left_join(dt_pas) %>% 
  mutate(seed = 1:nrow(.), 
         protection_cat_broad = "strict") %>% 
  st_transform(., crs = "ESRI:54009")

loose_pa_shapes_raw <- read_sf("data/spatial_data/protected_areas/pa_shapes.gpkg") %>% 
  dplyr::select(unique_id) %>% 
  filter(unique_id %in% unique(dt_loose_pas$unique_id)) %>%
  left_join(dt_pas) %>% 
  mutate(seed = 1:nrow(.), 
         protection_cat_broad = "loose") %>% 
  st_transform(., crs = "ESRI:54009")

#remove overlapping polygons (keep only the largest)
pa_shapes_strict <- resolve_overlaps(polygons = strict_pa_shapes_raw, keep_largest = TRUE)
pa_shapes_loose <- resolve_overlaps(polygons = loose_pa_shapes_raw, keep_largest = TRUE)

pa_shapes <- rbind(pa_shapes_strict, pa_shapes_loose) %>% 
  dplyr::select(-IUCN_CAT)

names(pa_shapes)

write_sf(pa_shapes, "data/spatial_data/protected_areas/pa_shapes_no_overlap.gpkg", append = FALSE)


dt_grid <- fread("data/processed_data/grid_with_covariates.csv")
grid_raw <- read_sf("data/spatial_data/grids/world_grid.gpkg") %>% 
  left_join(dt_grid %>% dplyr::select(-c(x_moll, y_moll,lon, lat )))

grid_not_prot <- grid_raw %>% 
  mutate(iucn_cat_ord = ifelse(is.na(iucn_cat_ord), "not_protected", iucn_cat_ord)) %>% 
  filter(iucn_cat_ord == "not_protected")

st_crs(grid_raw)==st_crs(pa_shapes) 
#pas_t <- pa_shapes %>% st_transform(st_crs(grid_not_prot))

pas_t <- pa_shapes #%>% sample_n(10)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

n_cores <- parallel::detectCores()
# Create and register a cluster
clust <- makeCluster(20)
registerDoSNOW(clust)

## progress bar 
iterations <- nrow(pas_t)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


pa_controls_raw <- NULL

pa_controls_raw <- foreach(i = 1:nrow(pas_t),
               .packages = c('sf', 'ggplot2', 'data.table', 'tidyverse'),
               .options.snow = opts,
               .inorder = FALSE,
               .verbose = TRUE,
               .errorhandling = "remove",
               .combine = rbind) %dopar% {

  
  pa <- pas_t[i,]
  
  mat <- pa$mat
  map <- pa$map
  biome <- pa$functional_biome
  cont <- pa$continent 
  Seed <- pa$seed 

  potential_space <- grid_not_prot %>% 
    filter(!mat > (mat+2) & !mat < (mat-2)) %>% 
    filter(!map > (map+250) & !map < (map-250)) %>% 
    filter(functional_biome == biome & continent == cont) 
  
  if(nrow(potential_space) == 0){return(NULL)}
  
  potential_space <- potential_space %>% 
    summarize()
  
  ## skip if there is no equivalent at all
  bBox <- st_bbox(potential_space)
  if(is.na(bBox$xmin)){return(NULL)}
 
  # create 100k buffer around pa
  sf_use_s2(FALSE)
  
  pa_buff_100k <- st_buffer(pa, dist = 100000)
  potential_space_100k <- st_intersection(potential_space, pa_buff_100k)

  new_poly <- NULL
  # try to place it within 100k 
  if(nrow(potential_space_100k) >= 1){
    
  new_poly <- move_polygon(og_poly = pa,
                           pot_space = potential_space_100k,
                           max_attempts = 1000, 
                           seed = Seed)

  dist_cat <- "100k"
  
  }
  
  # try to place it within 500k if it hasn't found anything withing 100k 
  
  if(is.null(new_poly)){
   
    pa_buff_500k <- st_buffer(pa, dist = 500000)
    potential_space_500k <- st_intersection(potential_space, pa_buff_500k)
    
    if(nrow(potential_space_500k) >= 1){
      
      new_poly <- move_polygon(og_poly = pa,
                               pot_space = potential_space_500k,
                               max_attempts = 1000, 
                               seed = Seed)
      dist_cat <- "500k"
    }
  }
  
  # try to place it within 1000k if it hasn't found anything withing 500k 
  
  if(is.null(new_poly)){
    
    pa_buff_1000k <- st_buffer(pa, dist = 1000000)
    potential_space_1000k <- st_intersection(potential_space, pa_buff_1000k)
    
    if(nrow(potential_space_1000k) >= 1){
      
      new_poly <- move_polygon(og_poly = pa,
                               pot_space = potential_space_1000k,
                               max_attempts = 1000, 
                               seed = Seed)
      dist_cat <- "1000k"
    }
  }
  
  #try to place anywhere if it hasn't found something within 1000k
  if(is.null(new_poly)){
      
      new_poly <- move_polygon(og_poly = pa,
                               pot_space = potential_space,
                               max_attempts = 1000, 
                               seed = Seed)
      dist_cat <- "over_1000k"
      
    }
  sf_use_s2(TRUE)
  
  
  # mapview(potential_space) + mapview(new_poly, col.regions = "red") + mapview(pa, col.regions = "green")
  
  if(!is.null(new_poly)){
     new_poly <- new_poly %>% 
       mutate(control_for = pas_t[i,]$unique_id, 
              continent = cont, 
              functional_biome = biome, 
              control_within_dist = dist_cat)
     
  }
  
  return(new_poly)
  
}

toc()
stopCluster(clust)
print(paste0("Loop done! ", Sys.time()))

pa_controls_raw_strict <- pa_controls_raw %>% 
  filter(control_for %in% unique(dt_strict_pas$unique_id))

pa_controls_raw_loose <- pa_controls_raw %>% 
  filter(control_for %in% unique(dt_loose_pas$unique_id))
  
pa_controls_strict <- resolve_overlaps(polygons = pa_controls_raw_strict, keep_largest = TRUE)
pa_controls_loose <- resolve_overlaps(polygons = pa_controls_raw_loose, keep_largest = TRUE)

pa_controls <- rbind(pa_controls_strict, pa_controls_loose)

print(Sys.time())
print(paste0("Found controls for ", 
             round((nrow(pa_controls)/nrow(pas_t)*100), 1), "% of the PAs (", nrow(pa_controls), " in total)"))
#"Found controls for 52.6% of the PAs (25461 in total)"


print(paste0("Found controls for ", 
             round((nrow(pa_controls %>% 
                           filter(control_for %in% unique(dt_strict_pas$unique_id)))/
                      nrow(pas_t %>% 
                             filter(unique_id %in% unique(dt_strict_pas$unique_id)))*100), 1),
             "% of the PAs (", nrow(pa_controls %>% 
                                      filter(control_for %in% unique(dt_strict_pas$unique_id))), " in total)"))
# "Found controls for 50.5% of the PAs (5345 in total)"
table(pa_controls$control_within_dist)

pa_controls_final <- pa_controls %>% 
  mutate(unique_id = paste0(control_for, "_control"))

write_sf(pa_controls_final, "data/spatial_data/protected_areas/controls_for_pas.gpkg", append = FALSE)
n_distinct(pa_controls_final$unique_id)
