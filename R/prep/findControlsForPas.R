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

paShapes <- paShapes %>% sample_n(100)

paControls <- NULL

for(i in 1:nrow(paShapes)){
  
  
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
  if(is.na(bBox$xmin)){next}
 
  # mapview(potentialSpace) + mapview(pa, col.regions = "red")
  
  newPoly <- movePolygon(ogPoly = pa, potSpace = potentialSpace, maxAttempts = 1000)
 
  # mapview(potentialSpace) + mapview(newPoly, col.regions = "red")
  
  
  if(is.null(paControls)){
   
     newPoly <- newPoly %>% 
      mutate(controlFor = paShapes[i,]$unique_id, 
             continent = cont, 
             Biome = biome)
    
    paControls <- newPoly
  }
  if(!is.null(newPoly)){
   
     newPoly <- newPoly %>% 
       mutate(controlFor = paShapes[i,]$unique_id, 
              continent = cont, 
              Biome = biome)
     
     paControls <- rbind(paControls, newPoly)
  }
  
  print(paste0(i, "/", nrow(paShapes), " done"))
}



mapview(paControls)
