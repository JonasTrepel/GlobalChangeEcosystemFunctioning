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

shp1.raw <- read_sf("data/spatialData/protectedAreas/rawPAs/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_0/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

table(shp1.raw$IUCN_CAT)
shp1 <- shp1.raw %>% 
  filter(!grepl("arine", DESIG_ENG)) 


shp2.raw <- read_sf("data/spatialData/protectedAreas/rawPAs/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_1/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

shp2 <- shp2.raw %>% 
  filter(!grepl("arine", DESIG_ENG))


shp3.raw <- read_sf("data/spatialData/protectedAreas/rawPAs/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_2/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

shp3 <- shp3.raw %>% 
  filter(!grepl("arine", DESIG_ENG)) 

pas <- rbind(shp1, shp2, shp3)
names(pas)
pas <- pas[!duplicated(data.frame(pas)),] %>% 
  mutate(iucnCatOrd = case_when(
    IUCN_CAT == "Ia" ~  1,  
    IUCN_CAT == "Ib" ~  2,  
    IUCN_CAT == "II" ~  3,  
    IUCN_CAT == "III" ~  4,  
    IUCN_CAT == "IV" ~  5,  
    IUCN_CAT == "V" ~  6,  
    IUCN_CAT == "VI" ~  7,  
    IUCN_CAT %in% c("Not Applicable", "Not Reported","Not Assigned") ~  8,  
  )) %>% 
  filter(MARINE == 0) %>% 
  dplyr::select(iucnCatOrd, IUCN_CAT, DESIG_ENG, WDPA_PID, NAME)

## Load world and make Grid
world <- rnaturalearth::ne_countries() %>% dplyr::select(continent, sovereignt, geometry)

#world <- world %>% filter(sovereignt == "Germany")

## For now 5km
worldGridRaw <- st_make_grid(world, cellsize = c(0.045, 0.045)) 
worldGridRaw1 <- worldGridRaw %>% st_as_sf()

mapview::mapview(worldGridRaw1)
sf_use_s2(FALSE)
worldGridRaw2 <- worldGridRaw1 %>%
  st_join(world) %>%
  filter(!is.na(sovereignt)) %>% 
  mutate(gridID = paste0("grid", 1:nrow(.))) %>% 
  st_join(pas) %>% 
  dplyr::select(-c(WDPA_PID, NAME,DESIG_ENG, IUCN_CAT)) %>%
  group_by(gridID) %>% 
  slice_min(iucnCatOrd) %>% distinct()

mapview::mapview(worldGridRaw2, zcol = "iucnCatOrd")

#pasGE <- pas %>% st_join(world) %>% filter(sovereignt == "Germany")
#mapview(pasGE)

# get area of intesection (only partially covered: remove)

write_sf(worldGridRaw2, "data/spatialData/protectedAreas/paGrid.gpkg", append = FALSE)

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
  cbind(worldGridRaw2[, "gridID"]) %>% 
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
  "RelprecDiff", ## Relative prec change
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
                                         dat <- df[!is.na(df$value) & df$coverage_fraction > 0.75, ]
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
