Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}




library(MuMIn)
library(glmmTMB)
library(terra)
library(tidyverse)
library(data.table)
library(sf)
library(tidylog)
library(GGally)
library(scales)
library(rnaturalearth)
library(gridExtra)
library("qgam")
library(sf)
library(terra)


pas <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg")


## Covariates 

#### Environmental ####

## Elevation

## MAP

## MAT

## MaxTemp

#### Trends ####

## EVI trend

## Fire frequency trend

## SOS Trend 

## EOS Trend

## LOS Trend 

#### Global Change ####

## Degree increase 

## Percentage temp increase 

## Prec mm increase 

## Percentage mm prec increase 

## Temperature slope last 20 years

## Precipitation slope last 20 years 

## Nitrogen depo 

## Megafauna species loss 

## Body mass loss

#### Others ####

## Biome

## Landcover 


#### EXTRACT CONTINUOUS COVS #### ----------------------------------

colNames <- c(
  #### Environmental ####
  "Elevation", ## Elevation
  "MAP", ## MAP
  "MAT", ## MAT
  "MaxTemp", ## MaxTemp
  #### Trends ####
  "EVITrend", ## EVI trend
  "FireFreqTrend", ## Fire frequency trend
  "SOSTrend", ## SOS Trend 
  "EOSTrend", ## EOS Trend
  "LOSTrend", ## LOS Trend 
  #### Global Change ####
  "TempDiffDeg", ## Degree increase 
  "PrecDiffMM", ## Precipitation increase
  "RelTempDiff", ## Percentage temp increase 
  "RelPrecDiff", ## Percentage precipitation increase 
  "TempSlope", ## Temperature slope last 20 years
  "PrecSlope", ## Precipitation slope last 20 years 
  "NitrogenDepo", ## Nitrogen depo
  
)

covPaths <- c(
  #### Environmental ####
  "../../../../resources/spatial/Elevation/Elevation_Global_930m.tif", ## Elevation
  "../../../../resources/spatial/Chelsa_Climate/CHELSA_bio12_1981-2010_V.2.1.tif", ## MAP
  "../../../../resources/spatial/Chelsa_Climate/CHELSA_bio1_1981-2010_V.2.1.tif", ## MAT
  NA, ## MaxTemp
  #### Trends ####
  NA, ## EVI trend
  NA, ## Fire frequency trend
  NA, ## SOS Trend 
  NA, ## EOS Trend
  NA, ## LOS Trend 
  #### Global Change ####
  NA, ## Degree increase 
  NA, ## Precipitation increase
  NA, ## Percentage temp increase 
  NA, ## Percentage precipitation increase 
  NA, ## Temperature slope last 20 years
  NA, ## Precipitation slope last 20 years 
  "../../../../resources/spatial/N_deposition_Rubin_etal/total_N_dep.tif"## Nitrogen depo 
)

covs <- data.table(
  colName = colNames, 
  covPath = covPaths
) %>% filter(!is.na(covPaths))


paCovs <- pas %>% as.data.table() %>% mutate(geom = NULL)

for(i in 1:nrow(covs)){
  
  covR <- rast(covs[i, ]$covPath)
  
  paTrans <- st_transform(pas, crs = st_crs(covR))
  
  extr <- exactextractr::exact_extract(covR, 
                                       paTrans, 
                                       fun = "mean")
  extrDT <- data.table(
    extrCol = extr
  )
  setnames(extrDT, "extrCol", covs[i, ]$colName)
  
  paCovs <- cbind(paCovs, extrDT)
  
  
}


#### EXTRACT CATEGORICAl COVS #### ----------------------------------

## Biome

wwfBiome <- read_sf("../../../../resources/spatial/WWF_Olson_Ecoregions/WWF_BIOMES.gpkg")

biomeLeg <- wwfBiome %>% as.data.table() %>% mutate(geom = NULL) %>% unique()

rastTmp <- terra::rast(res = 0.045) #roughly 5 km at equator 
biomeR <- terra::rasterize(wwfBiome, rastTmp, field = "BIOME")
plot(biomeR)

paTransBiome <- pas %>%
  st_transform(crs(biomeR))

biomeExtr <- terra::extract(biomeR,
                             paTransBiome,
                             fun = Mode, na.rm = T)# 


biomeExtrFin <- biomeExtr %>% 
  as.data.table() %>% 
  cbind(pas[, "unique_id"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(biomeLeg, by = "BIOME") %>% 
  rename(Biome = BIOME_Name) %>% 
  dplyr::select(unique_id, Biome) 

## Landcover 

