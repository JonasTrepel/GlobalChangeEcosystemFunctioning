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
  rename(PaAreaKm2 = GIS_AREA, 
         PaYear = STATUS_YR) %>% 
  mutate(PaYear = ifelse(PaYear == 0, NA, PaYear), 
         PaAge = 2023-PaYear) %>% 
  filter(!LandCover == "Water")


paShapes <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg") %>% 
  filter(unique_id %in% unique(pasCovsDTRaw$unique_id)) %>%
  left_join(pasCovsDTRaw)


gridNotProt <- read_sf("data/spatialData/gridWithCovs.gpkg") %>% 
  filter(iucnCat == "Not Protected")

paShapes <- paShapes %>% sample_n(25)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

nCores <- parallel::detectCores()-3
# Create and register a cluster
clust <- makeCluster(nCores)
registerDoSNOW(clust)

paControls <- NULL

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
  mmp <- pa$MMP
  biome <- pa$Biome
  cont <- pa$continent 
  
  
  potentialSpace <- gridNotProt %>% 
    filter(!MAT > (mat+2) & !MAT < (mat-2)) %>% 
    filter(!MMP > (mmp+50) & !MMP < (mmp-50)) %>% 
    filter(Biome == biome & continent == cont) %>% 
    summarize()
  
  ## skip if there is equivalent at all
  bBox <- st_bbox(potentialSpace)
  if(is.na(bBox$xmin)){return(NULL)}
 
  # mapview(potentialSpace) + mapview(pa, col.regions = "red")
  
  newPoly <- movePolygon(ogPoly = pa, potSpace = potentialSpace, maxAttempts = 5000)
 
  # mapview(potentialSpace) + mapview(newPoly, col.regions = "red")
  
  if(!is.null(newPoly)){
   
     newPoly <- newPoly %>% 
       mutate(controlFor = paShapes[i,]$unique_id, 
              continent = cont, 
              Biome = biome)
  }
  
  return(newPoly)
  
  print(paste0(i, "/", nrow(paShapes), " done, ", nrow(paControls), " of them successful (", round((nrow(paControls)/i*100), 1), "%)"))
  
}



print(paste0("Loop done! Found controls for ", 
             round((nrow(paControls)/nrow(paShapes)*100), 1), "% of the PAs (", nrow(paControls), " in total)"))

toc()
stopCluster(clust)
mapview(paControls)

write_sf(paControls, "data/spatialData/protectedAreas/controlsForPas.gpkg")