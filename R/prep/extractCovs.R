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
  "MMP", ## MAP
  "MAT", ## MAT
  "MaxTemp", ## MaxTemp
  "MinTemp", 
  #### Trends ####
  "EVITrend", ## EVI trend
  "FireFreqTrend", ## Fire frequency trend
  "SOSTrend", ## SOS Trend 
  "EOSTrend", ## EOS Trend
  "LOSTrend", ## LOS Trend 
  #### Global Change ####
  "AbsMeanTempDiff", ## Mean temp change 
  "RelMeanTempDiff", ## Relative mean temp change
  "AbsMaxTempDiff", ## Max temp change
  "RelMaxTempDiff", ## Relative Max temp change
  "AbsMinTempDiff", ## Min temp change
  "RelMinTempDiff", ## Relative min temp change
  "AbsPrecDiff", ## Prec change
  "RelprecDiff", ## Relative prec change
  "SlopeMeanTemp", ## Mean temp slope
  "SlopeMeanTemp", ## Max temp slope
  "SlopeMeanTemp", ## Min temp slope
  "SlopePrec", ## Prec slope
  "NitrogenDepo", ## Nitrogen depo
  "SpeciesLoss", ## Species loss
  "BodyMassLoss" ## Body mass loss
  
)

covPaths <- c(
  #### Environmental ####
  "data/spatialData/otherCovariates/Elevation_Global_930m.tif", ## Elevation
  "data/spatialData/climateData/currentMeanMonthlyPrec20092019.tif", ## MMP
  "data/spatialData/climateData/currentMeanTemp20092019.tif", ##MAT
  "data/spatialData/climateData/currentMaxTemp20092019.tif", ##Max Temp
  "data/spatialData/climateData/currentMinTemp20092019.tif", ## MinTemp
  #### Trends ####
  "data/spatialData/trendData/eviTrend20032023.tif", ## EVI trend
  "data/spatialData/trendData/fireFreqTrend20032023.tif", ## Fire frequency trend
  "data/spatialData/trendData/sosTrend20032023.tif", ## SOS Trend 
  "data/spatialData/trendData/eosTrend20032023.tif", ## EOS Trend
  "data/spatialData/trendData/losTrend20032023.tif", ## LOS Trend 
  #### Global Change ####
  "data/spatialData/climateData/absoluteChangeMeanTemp.tif", ## Mean temp change 
  "data/spatialData/climateData/relativeChangeMeanTemp.tif", ## Relative mean temp change
  "data/spatialData/climateData/absoluteChangeMaxTemp.tif", ## Max temp change
  "data/spatialData/climateData/relativeChangeMaxTemp.tif", ## Relative Max temp change
  "data/spatialData/climateData/absoluteChangeMinTemp.tif", ## Min temp change
  "data/spatialData/climateData/relativeChangeMinTemp.tif", ## Relative min temp change
  "data/spatialData/climateData/absoluteChangeMonthlyPrec.tif", ## Prec change
  "data/spatialData/climateData/relativeChangeMonthlyPrec.tif", ## Relative prec change
  "data/spatialData/climateData/slopeMeanTempChelsa20002019.tif", ## Mean temp slope
  NA, ## Max temp slope
  NA, ## Min temp slope
  NA, ## Prec slope
  "data/spatialData/otherCovariates/total_N_dep.tif",## Nitrogen depo
  "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Processed Data Layers/DEFAUNATION_MASS.tif", ## Species loss
  "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Processed Data Layers/DEFAUNATION_RICHNESS.tif" ## Body mass loss
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

