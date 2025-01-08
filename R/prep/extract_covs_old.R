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
#  "EviTrend", ## EVI trend
#  "BurnedAreaTrend", ## Burned area trend
#  "SOSTrend", ## SOS Trend 
#  "EviTrendMasked", ## EVI trend masked
#  "BurnedAreaTrendMasked", ## Burned area trend masked
#  "SOSTrendMasked", ## SOS Trend masked

  #### Global Change ####
#  "SlopeMeanTemp", ## Mean temp slope
#  "SlopeMaxTemp", ## Max temp slope
#  "SlopePrec", ## Prec slope
#  "SlopeMeanTempMasked", ## Mean temp slope masked
#  "SlopeMaxTempMasked", ## Max temp slope masked
#  "SlopePrecMasked", ## Prec slope masked
  "NitrogenDepo", ## Nitrogen depo
  "HumanModification" #Human modification index
  
  #### R-squared ####
#  "EviTrendR2", ## EVI trend R2
#  "BurnedAreaTrendR2", ## Burned Area trend R2
#  "SOSTrendR2", ## SOS trend R2 
#  "SlopeMeanTempR2", ## MAT trend R2
#  "SlopeMaxTempR2", ## MaxTemp trend R2
#  "SlopePrecR2" ## Prec trend R2
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
#  "data/spatialData/trendData/EviTrend20012023.tif", ## EVI trend
#  "data/spatialData/trendData/BurnedAreaTrend20012023.tif", ## Burned area trend
#  "data/spatialData/trendData/SosTrend20012023.tif", ## SOS Trend 
#  "data/spatialData/trendData/EviTrend20012023Masked.tif", ## EVI trend Masked
#  "data/spatialData/trendData/BurnedAreaTrend20012023Masked.tif", ## Burned area trend Masked
#  "data/spatialData/trendData/SosTrend20012023Masked.tif", ## SOS Trend Masked

  #### Global Change ####
#  "data/spatialData/trendData/MatTrend19502023.tif", ## Mean temp slope
#  "data/spatialData/trendData/MaxTempTrend19502023.tif", ## Max temp slope
#  "data/spatialData/trendData/MapTrend19502023.tif", ## Prec slope
#  "data/spatialData/trendData/MatTrend19502023Masked.tif", ## Mean temp slope masked
#  "data/spatialData/trendData/MaxTempTrend19502023Masked.tif", ## Max temp slope masked
#  "data/spatialData/trendData/MapTrend19502023Masked.tif", ## Prec slope masked
  "data/spatialData/otherCovariates/total_N_dep.tif",## Nitrogen depo
  "data/spatialData/otherCovariates/lulc-human-modification-terrestrial-systems_geographic.tif" #Human modification index
  
  #### R-squared ####
#  "data/spatialData/trendData/EviTRsq20012023.tif", ## EVI trend R2
#  "data/spatialData/trendData/BurnedAreaTRsq20012023.tif", ## Burned Area trend R2
#  "data/spatialData/trendData/SosTRsq20012023.tif", ## SOS trend R2 
#  "data/spatialData/trendData/MatRsq19502023.tif", ## MAT trend R2
#  "data/spatialData/trendData/MaxTempRsq19502023.tif", ## MaxTemp trend R2
#  "data/spatialData/trendData/MapRsq19502023.tif" ## Prec trend R2
)

covs <- data.table(
  colName = colNames, 
  covPath = covPaths
) %>% filter(!is.na(covPaths))


pasRawCovs <- pas %>% as.data.table() %>% mutate(geom = NULL)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

# Create and register a cluster
clust <- makeCluster(9)
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
  
  worldGridTrans <- st_transform(pas, crs = st_crs(covR))
  
  extr <- exactextractr::exact_extract(covR, 
                                       worldGridTrans, 
                                       fun = "mean")
  extrDT <- data.table(
    extrCol = extr
  )
  setnames(extrDT, "extrCol", covs[i, ]$colName)
  
  extraDTFin <- cbind(worldGridTrans[, "unique_id"], extrDT) %>%
    as.data.table() %>%
    mutate(geom = NULL) %>% 
    unique()
  return(extraDTFin)
  
  #print(paste0(i, "/", nrow(covs)))
  
}




pasCovsDT <- pasRawCovs %>% 
  left_join(dtCovs) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL)


#### EXTRACT CATEGORICAl COVS #### ----------------------------------

## Biome ---------------

wwfBiome <- read_sf("data/spatialData/otherCovariates/WWF_BIOMES.gpkg")

biomeLeg <- wwfBiome %>% as.data.table() %>% mutate(geom = NULL) %>% unique()

rastTmp <- terra::rast(res = 0.01) #roughly 1 km at equator 
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
  cbind(pas[, "unique_id"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(biomeLeg, by = "BIOME") %>% 
  rename(Biome = BIOME_Name) %>% 
  dplyr::select(unique_id, Biome) 

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


worldGridFunBiome <- pas %>% st_as_sf() %>%
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
  cbind(pas[, "unique_id"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(funBiomeLeg, by = "num") %>% 
  dplyr::select(unique_id, FunctionalBiome) 

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
  cbind(sfClimRegTrans[, "unique_id"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(climRegLeg) %>% 
  dplyr::select(unique_id, ClimaticRegion) 


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
  cbind(pas[, "unique_id"]) %>% 
  rename(lcNum = `lcMode`) %>%
  mutate(x = NULL) %>% 
  left_join(lcLeg, by = "lcNum") %>% 
  dplyr::select(unique_id, LandCover) 


### load more covariates 

#megafaunaCovs <- fread("data/processedData/cleanData/pasWithMegafauna.csv")

### get coordinates  
sf_use_s2(FALSE)
coords <- st_coordinates(st_centroid(pas)) %>% as.data.table()
pas$Longitude <- coords$X
pas$Latitude <- coords$Y

dtCoords <- pas %>% as.data.table() %>% 
  mutate(geom = NULL) %>% 
  dplyr::select(Longitude, Latitude, unique_id)



######## Summarize and Write Out
pasCovsDT <- pasCovsDT %>% 
  left_join(biomeExtrFin %>% unique()) %>% 
  left_join(lcExtrFin %>% unique()) %>% 
 # left_join(megafaunaCovs) %>% 
  left_join(climRegExtrFin %>% unique() ) %>% 
  left_join(funBiomeExtrFin %>% unique()) %>% 
  left_join(dtCoords %>% unique()) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL, 
         BurnedAreaMean = ifelse(is.na(BurnedAreaMean), 0 , BurnedAreaMean),
  ) %>% filter(GIS_AREA > 1)


fwrite(pasCovsDT, "data/processedData/cleanData/pasWithCovs.csv")
table(pasCovsDT$Biome)
table(pasCovsDT$FunctionalBiome)


pa.shapes <- pasCovsDT %>% left_join(pas) %>% st_as_sf

mapview::mapview(pa.shapes[is.na(pa.shapes$Biome), ])

round(quantile(pa.shapes[is.na(pa.shapes$Biome), ]$GIS_AREA, c(1, .99, .5, .01, 0), na.rm = T), 2)
round(quantile(pa.shapes[is.na(pa.shapes$Biome), ]$GIS_AREA, c(1, .99, .5, .01, 0), na.rm = T), 2)
