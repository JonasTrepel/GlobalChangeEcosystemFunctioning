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

pasCovsDTRaw <- fread("data/processedData/cleanData/pasWithCovs.csv") %>% 
  mutate(STATUS_YR = ifelse(STATUS_YR == 0, NA, STATUS_YR), 
         PaAge = 2023-STATUS_YR) %>% 
  filter(!LandCover %in% c("Water", "Ocean", "Urban", "Cultivated"))


paShapes <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg") %>% 
  filter(unique_id %in% unique(pasCovsDTRaw$unique_id)) %>%
  left_join(pasCovsDTRaw)

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


paControls <- NULL

paControls <- foreach(i = 1:nrow(paShapes),
               .packages = c('sf', 'ggplot2', 'data.table', 'tidyverse'),
               .options.snow = opts,
               .inorder = FALSE,
               .combine = rbind) %dopar% {

  
  pa <- paShapes[i,]
  
  mat <- pa$MAT
  map <- pa$MAP
  biome <- pa$FunctionalBiome
  cont <- pa$Continent 
  
  
  potentialSpace <- gridNotProt %>% 
   # filter(!MAT > (mat+2) & !MAT < (mat-2)) %>% 
   # filter(!MAP > (map+250) & !MAP < (map-250)) %>% 
    filter(FunctionalBiomeShort == biome & continent == cont) %>% 
    summarize()
  
  ## skip if there is no equivalent at all
  bBox <- st_bbox(potentialSpace)
  if(is.na(bBox$xmin)){return(NULL)}
 
  # mapview(potentialSpace) + mapview(pa, col.regions = "red")
  sf_use_s2(FALSE)
  newPoly <- movePolygon(ogPoly = pa, potSpace = potentialSpace, maxAttempts = 1000)
  sf_use_s2(TRUE)
  
  # mapview(potentialSpace) + mapview(newPoly, col.regions = "red") + mapview(pa, col.regions = "green")
  
  if(!is.null(newPoly)){
   
     newPoly <- newPoly %>% 
       mutate(controlFor = paShapes[i,]$unique_id, 
              Continent = cont, 
              FunctionalBiome = biome)
  }
  
  return(newPoly)
  
}



print(paste0("Loop done! Found controls for ", 
             round((nrow(paControls)/nrow(paShapes)*100), 1), "% of the PAs (", nrow(paControls), " in total)"))
#Found controls for 75.9% of the PAs (11184 in total)"
toc()
stopCluster(clust)
mapview(paControls)

write_sf(paControls, "data/spatialData/protectedAreas/controlsForPas.gpkg")
