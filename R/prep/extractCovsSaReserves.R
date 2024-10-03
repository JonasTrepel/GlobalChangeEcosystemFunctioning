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



reserveMeta <- fread("data/processedData/dataFragments/ReserveDataSouthAfricaFinal.csv")
pas <- read_sf("data/spatialData/protectedAreas/reserve_shapes.gpkg") %>% 
  left_join(reserveMeta) %>% 
  filter(!is.na(reserve_name)) %>% 
  rename(MeanBodyMassCwm = CW_mean_species_body_mass, 
         MaxBodyMass = max_species_body_mass,
         HerbivoreSpeciesRichness = n_herbi_sp_reserve, 
         HerbivoreFunDiv = herbi_fun_div_distq1, 
         HerbivoreBiomassKgKm2 = herbi_biomass_kgkm2,
         GrazerBiomassKgKm2 = grazer_biomass_kgkm2,
         BrowserBiomassKgKm2 = browser_biomass_kgkm2,
         MixedFeederBiomassKgKm2 = mixed_feeder_biomass_kgkm2,
         ElephantBiomassKgKm2 = elephant_biomass_kgkm2 
         ) %>% 
  mutate(PaAge = 2023-establishment_year) %>%
  dplyr::select(
           reserve_name, area_ha, geom, PaAge,
           MeanBodyMassCwm, MaxBodyMass, HerbivoreSpeciesRichness, 
           HerbivoreFunDiv, HerbivoreBiomassKgKm2, GrazerBiomassKgKm2,
           BrowserBiomassKgKm2, MixedFeederBiomassKgKm2, ElephantBiomassKgKm2
         )


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


pasRawCovs <- pas %>% as.data.table() %>% mutate(geom = NULL)

for(i in 1:nrow(covs)){
  
  covR <- rast(covs[i, ]$covPath)
  
  worldGridTrans <- st_transform(pas, crs = st_crs(covR))
  
  extr <- exactextractr::exact_extract(covR, 
                                       worldGridTrans, 
                                       fun = "mean")
  extrDT <- data.table(
    extrCol = extr
  )
  setnames(extrDT, "extrCol", covs[i, ]$colName)
  
  pasRawCovs <- cbind(pasRawCovs, extrDT)
  
  print(paste0(i, "/", nrow(covs)))
  
}

pasCovsDT <- pasRawCovs %>% 
  as.data.table() %>% 
  mutate(x = NULL)


#### EXTRACT CATEGORICAl COVS #### ----------------------------------

## Biome

wwfBiome <- read_sf("data/spatialData/otherCovariates/WWF_BIOMES.gpkg")

biomeLeg <- wwfBiome %>% as.data.table() %>% mutate(geom = NULL) %>% unique()

rastTmp <- terra::rast(res = 0.01) #roughly 5 km at equator 
biomeR <- terra::rasterize(wwfBiome, rastTmp, field = "BIOME")
plot(biomeR)

worldGridBiome <- pas %>% st_as_sf() %>%
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
  cbind(pas[, "reserve_name"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(biomeLeg, by = "BIOME") %>% 
  rename(Biome = BIOME_Name) %>% 
  dplyr::select(reserve_name, Biome) 

table(biomeExtr)
sum(is.na(biomeExtr))

## Koppen Geiger Climatic Regions

climaticRegionsR <- rast("data/spatialData/otherCovariates/KoppenGeigerClimaticRegions1km.tif")

climRegLeg <- data.frame(ClimRegNum = 1:30,
                         ClimaticRegion = c(
                           rep("Tropical", 3),
                           rep("Arid", 4),
                           rep("Temperate", 9),
                           rep("Cold", 12),
                           rep("Polar", 2)
                         ))

sfClimRegTrans <- pas %>%
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
  cbind(sfClimRegTrans[, "reserve_name"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(climRegLeg) %>% 
  dplyr::select(reserve_name, ClimaticRegion) 

table(ClimaticRegion)
sum(is.na(ClimaticRegion))


## Land cover 

# extract land cover and remove urban and cultivated pixels
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

gridTransLC <- pas %>%
  st_transform(crs(lc))


lcExtr <- exactextractr::exact_extract(lc,
                                       gridTransLC,
                                       summarize_df = TRUE,
                                       fun = function(df){
                                         dat <- df[!is.na(df$value) & df$coverage_fraction > 0.5, ]
                                         uniqueValues <- unique(dat$value)
                                         mode <- uniqueValues[which.max(tabulate(match(dat$value, uniqueValues)))]
                                         return(mode)
                                       } #https://rdrr.io/cran/exactextractr/man/exact_extract.html, see under User-defined summary functions
)

sum(is.na(lcExtr))


lcExtrFin <- data.table(lcMode = lcExtr) %>% 
  cbind(pas[, "reserve_name"]) %>% 
  rename(lcNum = `lcMode`) %>%
  mutate(x = NULL) %>% 
  left_join(lcLeg, by = "lcNum") %>% 
  dplyr::select(reserve_name, LandCover) 

### get coordinates  
sf_use_s2(FALSE)
coords <- st_coordinates(st_centroid(pas)) %>% as.data.table()
pas$Longitude <- coords$X
pas$Latitude <- coords$Y

dtCoords <- pas %>% as.data.table() %>% 
  mutate(geom = NULL) %>% 
  dplyr::select(Longitude, Latitude, reserve_name)


######## Summarize and Write Out
pasCovsDT <- pasRawCovs %>% 
  left_join(biomeExtrFin) %>% 
  left_join(lcExtrFin) %>% 
  left_join(dtCoords) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL, 
         FireFreqTrend = ifelse(is.na(FireFreqTrend), 0 , FireFreqTrend),
         BurnedAreaMean = ifelse(is.na(BurnedAreaMean), 0 , BurnedAreaMean),
         FireFrequencyMean = ifelse(is.na(FireFrequencyMean), 0 , FireFrequencyMean)
  ) %>% filter(area_ha > 100)


fwrite(pasCovsDT, "data/processedData/cleanData/saReservesWithCovs.csv")
table(pasCovsDT$Biome)

pa.shapes <- pasCovsDT %>% left_join(pas) %>% st_as_sf

mapview::mapview(pa.shapes[is.na(pa.shapes$Biome), ])

round(quantile(pa.shapes[is.na(pa.shapes$Biome), ]$GIS_AREA, c(1, .99, .5, .01, 0), na.rm = T), 2)
round(quantile(pa.shapes[is.na(pa.shapes$Biome), ]$GIS_AREA, c(1, .99, .5, .01, 0), na.rm = T), 2)
