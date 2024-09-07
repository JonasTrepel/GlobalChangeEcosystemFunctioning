library(terra)
library(tidyverse)

############ BUILD CLIMATE CHANGE METRICS #############
########## Using Chelsa and Paleoclim Data ########### 



############ BASELINE: 1600-1900 ############ 

#mean Temperature (C) of 1600-1900 -------------
paleoMeanTFiles <- c("O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/BIOCLIM/CHELSA_TraCE21k_bio01_16_V1.0.tif",
                     "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/BIOCLIM/CHELSA_TraCE21k_bio01_17_V1.0.tif",
                     "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/BIOCLIM/CHELSA_TraCE21k_bio01_18_V1.0.tif")

paleoMeanT <- rast(paleoMeanTFiles)

preIndustMeanT <- mean(paleoMeanT, na.rm = T)
plot(preIndustMeanT)

#max Temperature (C) of 1600-1900 ------------
paleoMaxTFiles16 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/tasmax/", pattern = "_16_V1.0", full.names = T)
paleoMaxTFiles17 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/tasmax/", pattern = "_17_V1.0", full.names = T)
paleoMaxTFiles18 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/tasmax/", pattern = "_18_V1.0", full.names = T)

paleoMaxT <- rast(c(paleoMaxTFiles16, paleoMaxTFiles17, paleoMaxTFiles18))

preIndustMaxTRaw <- max(paleoMaxT, na.rm = T)
plot(preIndustMaxTRaw)
#for some reason the unit is K/10
preIndustMaxT <- ((preIndustMaxTRaw/10) -273.15)
plot(preIndustMaxT)

#min Temperature (C) of 1600-1900 ---------------
paleoMinTFiles16 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/tasmin/", pattern = "_16_V1.0", full.names = T)
paleoMinTFiles17 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/tasmin/", pattern = "_17_V1.0", full.names = T)
paleoMinTFiles18 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/tasmin/", pattern = "_18_V1.0", full.names = T)

paleoMinT <- rast(c(paleoMinTFiles16, paleoMinTFiles17, paleoMinTFiles18))

preIndustMinTRaw <- min(paleoMinT, na.rm = T)
plot(preIndustMinTRaw)
#for some reason the unit is K/10
preIndustMinT <- ((preIndustMinTRaw/10) -273.15)
plot(preIndustMinT)

#mean Prec of 1600-1900 ------------------------
paleoPrecTFiles16 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/pr/", pattern = "_16_V1.0", full.names = T)
paleoPrecTFiles17 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/pr/", pattern = "_17_V1.0", full.names = T)
paleoPrecTFiles18 <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/pr/", pattern = "_18_V1.0", full.names = T)

paleoPrecT <- rast(c(paleoPrecTFiles16, paleoPrecTFiles17, paleoPrecTFiles18))

preIndustPrec <- mean(paleoPrecT, na.rm = T)
plot(preIndustPrec)

############ Current: 2009-2019 ############ 

## Current mean
file.exists("O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tas")
currentMeanTPathsRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tas/",  full.names = T)
currentMeanTFilesRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tas/",  full.names = F)
currentMeanPaths <- data.frame(paths = currentMeanTPathsRaw, 
                               files = currentMeanTFilesRaw) %>% 
  mutate(year = gsub("CHELSA_tas_", "", files),
         year = gsub("_V.2.1.tif", "", year), 
         year = str_split_i(year, "_", 2) ) %>%
  filter(year %in% c(2009:2019)) %>% 
  dplyr::select(paths) %>%
  pull()

currentMeanTStack <- rast(currentMeanPaths)

currentMeanTRaw <- mean(currentMeanTStack, na.rm = T)
plot(currentMeanTRaw)

currentMeanT <- (currentMeanTRaw/10)
plot(currentMeanT)

## Current max
currentMaxTPathsRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tasmax/",  full.names = T)
currentMaxTFilesRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tasmax/",  full.names = F)
currentMaxPaths <- data.frame(paths = currentMaxTPathsRaw, 
                              files = currentMaxTFilesRaw) %>% 
  mutate(year = gsub("CHELSA_tasmax_", "", files),
         year = gsub("_V.2.1.tif", "", year), 
         year = str_split_i(year, "_", 2) ) %>%
  filter(year %in% c(2009:2019)) %>% 
  dplyr::select(paths) %>%
  pull()

currentMaxTStack <- rast(currentMaxPaths)

currentMaxTRaw <- max(currentMaxTStack, na.rm = T)
plot(currentMaxTRaw)

currentMaxT <- (currentMaxTRaw/10)
plot(currentMaxT)

## Current min
currentMinTPathsRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tasmin/",  full.names = T)
currentMinTFilesRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tasmin/",  full.names = F)
currentMinPaths <- data.frame(paths = currentMinTPathsRaw, 
                              files = currentMinTFilesRaw) %>% 
  mutate(year = gsub("CHELSA_tasmin_", "", files),
         year = gsub("_V.2.1.tif", "", year), 
         year = str_split_i(year, "_", 2) ) %>%
  filter(year %in% c(2009:2019)) %>% 
  dplyr::select(paths) %>%
  pull()

currentMinTStack <- rast(currentMinPaths)

currentMinTRaw <- min(currentMinTStack, na.rm = T)
plot(currentMinTRaw)

currentMinT <- (currentMinT/10)
plot(currentMinT)

## Current mean monthly precipitation 
currentPrecPathsRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_pr/",  full.names = T)
currentPrecFilesRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_pr/",  full.names = F)
currentPrecPaths <- data.frame(paths = currentPrecPathsRaw, 
                               files = currentPrecFilesRaw) %>% 
  mutate(year = gsub("CHELSA_pr_", "", files),
         year = gsub("_V.2.1.tif", "", year), 
         year = str_split_i(year, "_", 2) ) %>%
  filter(year %in% c(2009:2019)) %>% 
  dplyr::select(paths) %>%
  pull()

currentPrecStack <- rast(currentPrecPaths)

currentPrec <- min(currentPrecStack, na.rm = T)
plot(currentPrec)

#write data to avoid redoing all these time-consuming calculations 
writeRaster(preIndustMeanT, "data/spatialData/climateData/preIndustrialMeanTemp16001900.tif")
writeRaster(preIndustMaxT , "data/spatialData/climateData/preIndustrialMaxTemp16001900.tif")
writeRaster(preIndustMinT , "data/spatialData/climateData/preIndustrialMinTemp16001900.tif")
writeRaster(preIndustPrec, "data/spatialData/climateData/preIndustrialMeanMonthlyPrec16001900.tif")
writeRaster(currentMeanT , "data/spatialData/climateData/currentMeanTemp20092019.tif")
writeRaster(currentMaxT , "data/spatialData/climateData/currentMaxTemp20092019.tif")
writeRaster(currentMinT, "data/spatialData/climateData/currentMinTemp20092019.tif")
writeRaster(currentPrec , "data/spatialData/climateData/currentMeanMonthlyPrec20092019.tif")


NA, ## Degree increase 
NA, ## Precipitation increase
NA, ## Percentage temp increase 
NA, ## Percentage precipitation increase 
NA, ## Temperature slope last 20 years
NA, ## Precipitation slope last 20 years 
NA, ## temp max slope last 20 years
NA, ## temp min slope last 20 years
NA, ## temp max diff last 20 years
NA, ## temp min diff last 20 years