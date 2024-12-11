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

source("R/functions/movePolygon.R")
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


file.exists("data/spatialData/protectedAreas/gridWithCovs.gpkg")

gridNotProt_raw <- read_sf("data/spatialData/protectedAreas/gridWithCovs.gpkg") 

gridNotProt <- gridNotProt_raw %>% 
  mutate(iucnCat = ifelse(is.na(iucnCat), "NotProtected", iucnCat)) %>% 
  filter(iucnCat == "NotProtected")

#paShapes <- paShapes %>% sample_n(25)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

nCores <- parallel::detectCores()/2
# Create and register a cluster
clust <- makeCluster(nCores)
registerDoSNOW(clust)

## progress bar 
iterations <- nrow(paShapes)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


pa_controls_raw <- NULL

pa_controls_raw <- foreach(i = 1:nrow(paShapes),
               .packages = c('sf', 'ggplot2', 'data.table', 'tidyverse'),
               .options.snow = opts,
               .inorder = FALSE,
               .combine = rbind) %dopar% {

  
  pa <- paShapes[i,]
  
  mat <- pa$MAT
  map <- pa$MAP
  biome <- pa$FunctionalBiome
  cont <- pa$Continent 
  Seed <- pa$seed 
  
  potentialSpace <- gridNotProt %>% 
    filter(!MAT > (mat+2) & !MAT < (mat-2)) %>% 
    filter(!MAP > (map+250) & !MAP < (map-250)) %>% 
    filter(FunctionalBiomeShort == biome & continent == cont) %>% 
    summarize()
  
  ## skip if there is no equivalent at all
  bBox <- st_bbox(potentialSpace)
  if(is.na(bBox$xmin)){return(NULL)}
 
  # mapview(potentialSpace) + mapview(pa, col.regions = "red")
  sf_use_s2(FALSE)
  newPoly <- movePolygon(ogPoly = pa,
                         potSpace = potentialSpace,
                         maxAttempts = 1000, 
                         seed = Seed)
  sf_use_s2(TRUE)
  
  # mapview(potentialSpace) + mapview(newPoly, col.regions = "red") + mapview(pa, col.regions = "green")
  
  if(!is.null(newPoly)){
   
     newPoly <- newPoly %>% 
       mutate(control_for = paShapes[i,]$unique_id, 
              Continent = cont, 
              FunctionalBiome = biome)
  }
  
  return(newPoly)
  
}

toc()
stopCluster(clust)
print(paste0("Loop done!"))

      
pa_controls <- resolve_overlaps(polygons = pa_controls_raw, keep_largest = TRUE)


print(paste0("Found controls for ", 
             round((nrow(pa_controls)/nrow(paShapes)*100), 1), "% of the PAs (", nrow(pa_controls), " in total)"))

#"Found controls for 74% of the PAs (7751 in total)"

pa_controls_final <- pa_controls %>% 
 # rename(control_for = controlFor) %>%
  mutate(unique_id = paste0(control_for, "_control"))

write_sf(pa_controls_final, "data/spatialData/protectedAreas/controls_for_strict_pas.gpkg")
