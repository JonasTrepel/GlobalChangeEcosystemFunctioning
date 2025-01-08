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

pasCovsDTRaw <- fread("data/processedData/cleanData/pasWithCovs.csv") %>% 
  mutate(STATUS_YR = ifelse(STATUS_YR == 0, NA, STATUS_YR), 
         PaAge = 2023-STATUS_YR) %>% 
  filter(!LandCover %in% c("Water", "Ocean", "Urban", "Cultivated") & 
           IUCN_CAT %in% c("Ib", "Ia", "II")) %>% unique()


paShapesRaw <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg") %>% 
  dplyr::select(unique_id) %>% 
  filter(unique_id %in% unique(pasCovsDTRaw$unique_id)) %>%
  left_join(pasCovsDTRaw) %>% 
  mutate(seed = 1:nrow(.))

#remove overlapping polygons (keep only the largest)
paShapes <- resolve_overlaps(polygons = paShapesRaw, keep_largest = TRUE)

write_sf(paShapes, "data/spatialData/protectedAreas/strict_pas.gpkg")

pasCovsDT <- pasCovsDTRaw %>% filter(unique_id %in% unique(paShapes$unique_id))

fwrite(pasCovsDT, "data/processedData/cleanData/strict_pas_with_covs.csv")

#paShapes <- st_read("data/spatialData/protectedAreas/strict_pas.gpkg")


file.exists("data/spatialData/protectedAreas/gridWithCovs.gpkg")

gridNotProt_raw <- read_sf("data/spatialData/protectedAreas/gridWithCovs.gpkg") 

gridNotProt <- gridNotProt_raw %>% 
  mutate(iucn_cat_ord = ifelse(is.na(iucn_cat_ord), "not_protected", iucn_cat_ord)) %>% 
  filter(iucn_cat_ord == "not_protected")

pas_t <- paShapes %>% st_transform(st_crs(gridNotProt))

#pas_t <- pas_t %>% sample_n(10)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

nCores <- parallel::detectCores()
# Create and register a cluster
clust <- makeCluster(25)
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
  
  mat <- pa$MAT
  map <- pa$MAP
  biome <- pa$FunctionalBiome
  cont <- pa$Continent 
  Seed <- pa$seed 
  
  potential_space <- gridNotProt %>% 
    filter(!MAT > (mat+2) & !MAT < (mat-2)) %>% 
    filter(!MAP > (map+250) & !MAP < (map-250)) %>% 
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
              Continent = cont, 
              FunctionalBiome = biome, 
              control_within_dist = dist_cat)
  }
  
  return(new_poly)
  
}

toc()
stopCluster(clust)
print(paste0("Loop done!"))

      
pa_controls <- resolve_overlaps(polygons = pa_controls_raw, keep_largest = TRUE)


print(paste0("Found controls for ", 
             round((nrow(pa_controls)/nrow(pas_t)*100), 1), "% of the PAs (", nrow(pa_controls), " in total)"))

#"Found controls for 74% of the PAs (7751 in total)"

pa_controls_final <- pa_controls %>% 
 # rename(control_for = controlFor) %>%
  mutate(unique_id = paste0(control_for, "_control"))

write_sf(pa_controls_final, "data/spatialData/protectedAreas/controls_for_strict_pas.gpkg")
