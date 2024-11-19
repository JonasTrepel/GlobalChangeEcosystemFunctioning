Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


library(tidyverse)
library(sf)
library(rnaturalearth)
library(mapview)
library(terra)
library(data.table)


pas <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg")


## Load world and make Grid
world <- rnaturalearth::ne_countries() %>% dplyr::select(continent, sovereignt, geometry)

## For now 5km
worldGridRaw <- st_make_grid(world, cellsize = c(0.045, 0.045)) 
worldGridRaw1 <- worldGridRaw %>% st_as_sf()

sf_use_s2(FALSE)
worldGridRaw2 <- worldGridRaw1 %>%
  st_join(world) %>%
  filter(!is.na(sovereignt)) %>% 
  mutate(gridID = paste0("grid", 1:nrow(.))) %>% 
  st_join(pas) %>% 
  dplyr::select(-c(WDPA_PID, NAME, DESIG_ENG, IUCN_CAT, Continent)) %>%
  group_by(gridID) %>% 
  slice_min(iucnCatOrd) %>% distinct()

table(worldGridRaw2$continent)

#mapview::mapview(worldGridRaw2, zcol = "continent")

# extract landcover and remove urban and cultivated pixels
lcNum <-  c(0, 20, 30, 40, 50, 60, 70, 80, 90, 
          100, 111, 112, 113, 114, 115, 116,
          121, 122, 123, 124, 125, 126, 
          200)
LandCover <- c("Unknown", 
          "Shrubs", 
          "HerbacousVeg", 
          "Cultivated",
          "Urban",
          "Bare",
          "SnowIce",
          "Water",
          "HerbaceousWetland",
          "MossAndLichen",
          "ClosedForestEvergreenNeedle",
          "ClosedForestEvergreenBroad",
          "ClosedForestDeciduousNeedle",
          "ClosedForestDeciduousBroad",
          "ClosedForestMixed",
          "ClosedForestOther",
          "OpenForestEvergreenNeedle",
          "OpenForestEvergreenBroad",
          "OpenForestDeciduousNeedle",
          "OpenForestDeciduousBroad",
          "OpenForestMixed",
          "OpenForestOther",
          "Ocean")

lcLeg <- data.frame(LandCover=LandCover, lcNum=lcNum)

lc <- rast("data/spatialData/otherCovariates/GlobalLandCovercopernicus2019.tif")
plot(lc)

gridTransLC <- worldGridRaw2 %>%
  st_transform(crs(lc))

lcExtr <- exactextractr::exact_extract(lc,
                                       gridTransLC,
                                       summarize_df = TRUE,
                                       fun = function(df){
                                         dat <- df[!is.na(df$value) & df$coverage_fraction > 0.75, ]
                                         uniqueValues <- unique(dat$value)
                                         mode <- uniqueValues[which.max(tabulate(match(dat$value, uniqueValues)))]
                                         return(mode)
                                       } #https://rdrr.io/cran/exactextractr/man/exact_extract.html, see under User-defined summary functions
                                       )


lcExtrFin <- data.table(lcMode = lcExtr) %>% 
  cbind(gridTransLC[, "gridID"]) %>% 
  rename(lcNum = `lcMode`) %>%
  mutate(x = NULL) %>% 
  left_join(lcLeg, by = "lcNum") %>% 
  dplyr::select(gridID, LandCover) 

worldGridRaw3 <- worldGridRaw2 %>%
  left_join(lcExtrFin) %>% 
  filter(!LandCover %in% c("Unknown", "Cultivated", "Urban", "Ocean", "Water")) %>% 
  mutate(
    iucnCat = case_when(
      iucnCatOrd == 1 ~ "Ia",  
      iucnCatOrd == 2 ~ "Ib",  
      iucnCatOrd == 3 ~ "II",  
      iucnCatOrd == 4 ~ "III",  
      iucnCatOrd == 5 ~ "IV",  
      iucnCatOrd == 6 ~ "V",  
      iucnCatOrd == 7 ~ "VI",  
      iucnCatOrd == 8 ~ "UnknownOrNA",  
    )
  )

unique(worldGridRaw3$LandCover)

## Biome ---------------

wwfBiome <- read_sf("data/spatialData/otherCovariates/WWF_BIOMES.gpkg")

biomeLeg <- wwfBiome %>% as.data.table() %>% mutate(geom = NULL) %>% unique()

rastTmp <- terra::rast(res = 0.01) #roughly 1 km at equator 
biomeR <- terra::rasterize(wwfBiome, rastTmp, field = "BIOME")
plot(biomeR)

worldGridBiome <- worldGridRaw3 %>% st_as_sf() %>%
  st_transform(crs(biomeR))

biomeExtr <- exactextractr::exact_extract(biomeR,
                                          worldGridBiome,
                                          summarize_df = TRUE,
                                          fun = function(df){
                                            dat <- df[!is.na(df$value) & df$coverage_fraction > 0.25, ]
                                            uniqueValues <- unique(dat$value)
                                            mode <- uniqueValues[which.max(tabulate(match(dat$value, uniqueValues)))]
                                            return(mode)
                                          } #https://rdrr.io/cran/exactextractr/man/exact_extract.html, see under User-defined summary functions
)


biomeExtrFin <- data.table(BIOME = biomeExtr) %>% 
  cbind(worldGridBiome[, "gridID"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(biomeLeg, by = "BIOME") %>% 
  rename(Biome = BIOME_Name) %>% 
  dplyr::select(gridID, Biome) 

table(biomeExtr)
sum(is.na(biomeExtr))


## Functional Biome ------------

funBiome <- rast("data/spatialData/otherCovariates/higginsFunctionalBiomes.tif")
plot(funBiome)
funBiomeLeg <- data.table(
  num = c(1:24), 
  FunctionalBiome = c("SLC","SMC","SHC","TLC","TMC","THC",
                      "SLD","SMD","SHD","TLD","TMD","THD",
                      "SLB","SMB","SHB","TLB","TMB","THB","SLN","SMN","SHN",
                      "TLN","TMN","THN")
)

worldGridFunBiome <- worldGridRaw3 %>% st_as_sf() %>%
  st_transform(crs(funBiome))

funBiomeExtr <- exactextractr::exact_extract(funBiome,
                                             worldGridFunBiome,
                                             summarize_df = TRUE,
                                             fun = function(df){
                                               dat <- df[!is.na(df$value) & df$coverage_fraction > 0.00000001, ]
                                               uniqueValues <- unique(dat$value)
                                               mode <- uniqueValues[which.max(tabulate(match(dat$value, uniqueValues)))]
                                               return(mode)
                                             } #https://rdrr.io/cran/exactextractr/man/exact_extract.html, see under User-defined summary functions
)

sum(is.na(funBiomeExtr))
funBiomeExtrFin <- data.table(num = funBiomeExtr) %>% 
  cbind(worldGridFunBiome[, "gridID"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(funBiomeLeg, by = "num") %>% 
  dplyr::select(gridID, FunctionalBiome) 

table(biomeExtr)
sum(is.na(biomeExtr))

## Koppen Geiger Climatic Regions ------------------------------

climaticRegionsR <- rast("data/spatialData/otherCovariates/KoppenGeigerClimaticRegions1km.tif")

climRegLeg <- data.frame(ClimRegNum = 1:30,
                         ClimaticRegion = c(
                           rep("Tropical", 3),
                           rep("Arid", 4),
                           rep("Temperate", 9),
                           rep("Cold", 12),
                           rep("Polar", 2)
                         ))

sfClimRegTrans <- worldGridRaw3 %>%
  st_transform(crs(climaticRegionsR))

climRegExtr <- exactextractr::exact_extract(climaticRegionsR,
                                            sfClimRegTrans,
                                            summarize_df = TRUE,
                                            fun = function(df){
                                              dat <- df[!is.na(df$value) & df$coverage_fraction > 0.25, ]
                                              uniqueValues <- unique(dat$value)
                                              mode <- uniqueValues[which.max(tabulate(match(dat$value, uniqueValues)))]
                                              return(mode)
                                            } #https://rdrr.io/cran/exactextractr/man/exact_extract.html, see under User-defined summary functions
)


climRegExtrFin <- data.table(ClimRegNum = climRegExtr) %>% 
  cbind(sfClimRegTrans[, "gridID"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(climRegLeg) %>% 
  dplyr::select(gridID, ClimaticRegion) 

worldGridRaw4 <- worldGridRaw3 %>% 
  left_join(funBiomeExtrFin %>% unique()) %>% 
  left_join(biomeExtrFin %>% unique()) %>% 
  left_join(climRegExtrFin %>% unique())

write_sf(worldGridRaw4, "data/spatialData/protectedAreas/paGrid.gpkg", append = FALSE)



#### EXTRACT CONTINUOUS COVS #### ----------------------------------

colNames <- c(
  #### Environmental ####
  "Elevation", ## Elevation
  "MAP", ## MAP
  "MAT", ## MAT
  "MaxTemp", ## MaxTemp
  "BurnedAreaMean", ## Burned Area
  "EVI", ## EVI
  
  #### Trends ####
  "EviTrend", ## EVI trend
  "BurnedAreaTrend", ## Burned area trend
  "SOSTrend", ## SOS Trend 
  "EviTrendMasked", ## EVI trend masked
  "BurnedAreaTrendMasked", ## Burned area trend masked
  "SOSTrendMasked", ## SOS Trend masked
  
  #### Global Change ####
  "SlopeMeanTemp", ## Mean temp slope
  "SlopeMaxTemp", ## Max temp slope
  "SlopePrec", ## Prec slope
  "SlopeMeanTempMasked", ## Mean temp slope masked
  "SlopeMaxTempMasked", ## Max temp slope masked
  "SlopePrecMasked", ## Prec slope masked
  "NitrogenDepo", ## Nitrogen depo
  "HumanModification", #Human modification index
  
  #### R-squared ####
  "EviTrendR2", ## EVI trend R2
  "BurnedAreaTrendR2", ## Burned Area trend R2
  "SOSTrendR2", ## SOS trend R2 
  "SlopeMeanTempR2", ## MAT trend R2
  "SlopeMaxTempR2", ## MaxTemp trend R2
  "SlopePrecR2" ## Prec trend R2
)

covPaths <- c(
  #### Environmental ####
  "data/spatialData/otherCovariates/Elevation_Global_930m.tif", ## Elevation
  "data/spatialData/climateData/currentMeanMonthlyPrec20092019.tif", ## MAP
  "data/spatialData/climateData/currentMeanTemp20092019.tif", ##MAT
  "data/spatialData/climateData/currentMaxTemp20092019.tif", ##Max Temp
  "data/spatialData/otherCovariates/BurnedAreaMean20012023.tif", ## Burned Area
  "data/spatialData/otherCovariates/EviMean20012023.tif", ## EVI
  
  #### Trends ####
  "data/spatialData/trendData/EviTrend20012023.tif", ## EVI trend
  "data/spatialData/trendData/BurnedAreaTrend20012023.tif", ## Burned area trend
  "data/spatialData/trendData/SosTrend20012023.tif", ## SOS Trend 
  "data/spatialData/trendData/EviTrend20012023Masked.tif", ## EVI trend Masked
  "data/spatialData/trendData/BurnedAreaTrend20012023Masked.tif", ## Burned area trend Masked
  "data/spatialData/trendData/SosTrend20012023Masked.tif", ## SOS Trend Masked
  
  #### Global Change ####
  "data/spatialData/trendData/MatTrend19502023.tif", ## Mean temp slope
  "data/spatialData/trendData/MaxTempTrend19502023.tif", ## Max temp slope
  "data/spatialData/trendData/MapTrend19502023.tif", ## Prec slope
  "data/spatialData/trendData/MatTrend19502023Masked.tif", ## Mean temp slope masked
  "data/spatialData/trendData/MaxTempTrend19502023Masked.tif", ## Max temp slope masked
  "data/spatialData/trendData/MapTrend19502023Masked.tif", ## Prec slope masked
  "data/spatialData/otherCovariates/total_N_dep.tif",## Nitrogen depo
  "data/spatialData/otherCovariates/lulc-human-modification-terrestrial-systems_geographic.tif", #Human modification index
  
  #### R-squared ####
  "data/spatialData/trendData/EviTRsq20012023.tif", ## EVI trend R2
  "data/spatialData/trendData/BurnedAreaTRsq20012023.tif", ## Burned Area trend R2
  "data/spatialData/trendData/SosTRsq20012023.tif", ## SOS trend R2 
  "data/spatialData/trendData/MatRsq19502023.tif", ## MAT trend R2
  "data/spatialData/trendData/MaxTempRsq19502023.tif", ## MaxTemp trend R2
  "data/spatialData/trendData/MapRsq19502023.tif" ## Prec trend R2
)

covs <- data.table(
  colName = colNames, 
  covPath = covPaths
) %>% filter(!is.na(covPaths))


gridCovsRaw <- worldGridRaw4 %>% as.data.table() %>% mutate(geom = NULL, x = NULL, geometry = NULL)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

# Create and register a cluster
clust <- makeCluster(30)
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

dtCovs <- foreach(i = 1:nrow(covs),
                  .packages = c('tidyverse', 'exactextractr', 'data.table', 'terra', 'sf'),
                  .options.snow = opts,
                  .inorder = FALSE,
                  .combine = left_join) %dopar% {
                    
                    #for(i in 1:nrow(covs)){
                    
                    covR <- rast(covs[i, ]$covPath)
                    
                    worldGridTrans <- st_transform(gridCovsRaw, crs = st_crs(covR))
                    
                    extr <- exactextractr::exact_extract(covR, 
                                                         worldGridTrans, 
                                                         fun = "mean")
                    extrDT <- data.table(
                      extrCol = extr
                    )
                    setnames(extrDT, "extrCol", covs[i, ]$colName)
                    
                    extraDTFin <- cbind(worldGridTrans[, "gridID"], extrDT) %>%
                      as.data.table() %>%
                      mutate(geom = NULL, x = NULL, geometry = NULL) %>% 
                      unique()
                    return(extraDTFin)
                    
                    #print(paste0(i, "/", nrow(covs)))
                    
                  }




gridCovsDT <- gridCovsRaw %>% 
  left_join(dtCovs %>% unique()) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL)

fwrite(gridCovsDT, "data/processedData/cleanData/gridWithCovs.csv")



#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################


#mapview(worldGridRaw3, zcol = "iucnCatOrd")

#### EXTRACT CONTINUOUS COVS #### ----------------------------------

colNames <- c(
  #### Environmental ####
  "Elevation", ## Elevation
  "MMP", ## MAP
  "MAT", ## MAT
  "MaxTemp", ## MaxTemp
  "MinTemp", ##MinTemp
  "BurnedAreaMean", ## Burned Area
  "FireFrequencyMean", ## Fire Frequency
  "EVI", ## EVI
  "NPP", ##
  "EviSD", ##
 
   #### Trends ####
  "EviTrend", ## EVI trend
  "NppTrend", ## NPP Trend
  "FireFreqTrend", ## Fire frequency trend
  "BurnedAreaTrend", ## Burned area trend
  "SOSTrend", ## SOS Trend 
  "EviSDTrend", ##EVI SD Trend

  #### Global Change ####
  "AbsMeanTempDiff", ## Mean temp change 
  "RelMeanTempDiff", ## Relative mean temp change
  "AbsMaxTempDiff", ## Max temp change
  "RelMaxTempDiff", ## Relative Max temp change
  "AbsMinTempDiff", ## Min temp change
  "RelMinTempDiff", ## Relative min temp change
  "AbsPrecDiff", ## Prec change
  "RelPrecDiff", ## Relative prec change
  "SlopeMeanTemp", ## Mean temp slope
  "SlopeMaxTemp", ## Max temp slope
  "SlopeMinTemp", ## Min temp slope
  "SlopePrec", ## Prec slope
  "NitrogenDepo", ## Nitrogen depo
  "HumanModification", #Human moddification index
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
  "data/spatialData/otherCovariates/BurnedAreaMean20012023.tif", ## Burned Area
  "data/spatialData/otherCovariates/FireMeanFreqTrend20012023.tif", ## Fire Frequency
  "data/spatialData/otherCovariates/EviMean20012023.tif", ## EVI
  "data/spatialData/otherCovariates/NppMean20012023.tif", ##NPP
  "data/spatialData/otherCovariates/EviSdMean20012023.tif", ## EVI Sd trend
  
  #### Trends ####
  "data/spatialData/trendData/EviTrend20012023.tif", ## EVI trend
  "data/spatialData/trendData/NppTrend20012023.tif", ## NPP Trend
  "data/spatialData/trendData/FireFreqTrend20012023.tif", ## Fire frequency trend
  "data/spatialData/trendData/BurnedAreaTrend20012023.tif", ## Burned area trend
  "data/spatialData/trendData/SosTrend20012023.tif", ## SOS Trend 
  "data/spatialData/trendData/EviSdTrend20012023.tif", ##EVI SD Trend
  
  #### Global Change ####
  "data/spatialData/climateData/absoluteChangeMeanTemp.tif", ## Mean temp change 
  "data/spatialData/climateData/relativeChangeMeanTemp.tif", ## Relative mean temp change
  "data/spatialData/climateData/absoluteChangeMaxTemp.tif", ## Max temp change
  "data/spatialData/climateData/relativeChangeMaxTemp.tif", ## Relative Max temp change
  "data/spatialData/climateData/absoluteChangeMinTemp.tif", ## Min temp change
  "data/spatialData/climateData/relativeChangeMinTemp.tif", ## Relative min temp change
  "data/spatialData/climateData/absoluteChangeMonthlyPrec.tif", ## Prec change
  "data/spatialData/climateData/relativeChangeMonthlyPrec.tif", ## Relative prec change
  "data/spatialData/climateData/slopeMeanTempChelsa20012019.tif", ## Mean temp slope
  "data/spatialData/climateData/slopeMaxTempChelsa20012019.tif", ## Max temp slope
  "data/spatialData/climateData/slopeMinTempChelsa20012019.tif", ## Min temp slope
  "data/spatialData/climateData/slopePrecChelsa20012019.tif", ## Prec slope
  "data/spatialData/otherCovariates/total_N_dep.tif",## Nitrogen depo
  "data/spatialData/otherCovariates/lulc-human-modification-terrestrial-systems_geographic.tif", #Human moddification index
  "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Processed Data Layers/DEFAUNATION_MASS.tif", ## Species loss
  "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Processed Data Layers/DEFAUNATION_RICHNESS.tif" ## Body mass loss
)

covs <- data.table(
  colName = colNames, 
  covPath = covPaths
) %>% filter(!is.na(covPaths))


worldGridRawCovs <- worldGridRaw3 %>% as.data.table() %>% mutate(geom = NULL)

for(i in 1:nrow(covs)){
  
  covR <- rast(covs[i, ]$covPath)
  
  worldGridTrans <- st_transform(worldGridRaw3, crs = st_crs(covR))
  
  extr <- exactextractr::exact_extract(covR, 
                                       worldGridTrans, 
                                       fun = "mean")
  extrDT <- data.table(
    extrCol = extr
  )
  setnames(extrDT, "extrCol", covs[i, ]$colName)
  
  worldGridRawCovs <- cbind(worldGridRawCovs, extrDT)
  
  print(paste0(i, "/", nrow(covs)))
  
}

worldGridCovsDT <- worldGridRawCovs %>% 
  as.data.table() %>% 
  mutate(x = NULL)

fwrite(worldGridCovsDT, "data/processedData/cleanData/gridWithCovs.csv")

worldGridCovs <- worldGridRawCovs

write_sf(worldGridCovs, "data/spatialData/gridWithCovs.gpkg")

#### EXTRACT CATEGORICAl COVS #### ----------------------------------

## Biome

wwfBiome <- read_sf("data/spatialData/otherCovariates/WWF_BIOMES.gpkg")

biomeLeg <- wwfBiome %>% as.data.table() %>% mutate(geom = NULL) %>% unique()

rastTmp <- terra::rast(res = 0.045) #roughly 5 km at equator 
biomeR <- terra::rasterize(wwfBiome, rastTmp, field = "BIOME")
plot(biomeR)

worldGridBiome <- worldGridRawCovs %>% st_as_sf() %>%
  st_transform(crs(biomeR))

biomeExtr <- exactextractr::exact_extract(biomeR,
                                          worldGridBiome,
                                       summarize_df = TRUE,
                                       fun = function(df){
                                         dat <- df[!is.na(df$value) & df$coverage_fraction > 0.5, ]
                                         uniqueValues <- unique(dat$value)
                                         mode <- uniqueValues[which.max(tabulate(match(dat$value, uniqueValues)))]
                                         return(mode)
                                       } #https://rdrr.io/cran/exactextractr/man/exact_extract.html, see under User-defined summary functions
)


biomeExtrFin <- data.table(BIOME = biomeExtr) %>% 
  cbind(worldGridRawCovs[, "gridID"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(biomeLeg, by = "BIOME") %>% 
  rename(Biome = BIOME_Name) %>% 
  dplyr::select(gridID, Biome) 


######## Summarize and Write Out
worldGridCovsDT <- worldGridRawCovs %>% 
  left_join(biomeExtrFin) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL, 
         iucnCat = ifelse(is.na(iucnCat), "Not Protected", iucnCat), 
         FireFreqTrend = ifelse(is.na(FireFreqTrend), 0 , FireFreqTrend),
         BurnedAreaMean = ifelse(is.na(BurnedAreaMean), 0 , BurnedAreaMean),
         FireFrequencyMean = ifelse(is.na(FireFrequencyMean), 0 , FireFrequencyMean)
         )

fwrite(worldGridCovsDT, "data/processedData/cleanData/gridWithCovs.csv")

worldGridCovs <- worldGridRawCovs %>% 
  left_join(biomeExtrFin) %>% 
  mutate(iucnCat = ifelse(is.na(iucnCat), "Not Protected", iucnCat),
         FireFreqTrend = ifelse(is.na(FireFreqTrend), 0 , FireFreqTrend),
         BurnedAreaMean = ifelse(is.na(BurnedAreaMean), 0 , BurnedAreaMean),
         FireFrequencyMean = ifelse(is.na(FireFrequencyMean), 0 , FireFrequencyMean))

write_sf(worldGridCovs, "data/spatialData/gridWithCovs.gpkg", append = FALSE)
