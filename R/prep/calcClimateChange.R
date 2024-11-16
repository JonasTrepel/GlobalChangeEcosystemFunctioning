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

currentMeanT <- ((currentMeanTRaw/10) -273.15)
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

currentMaxT <- ((currentMaxTRaw/10) -273.15)
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

currentMinT <- ((currentMinTRaw/10) -273.15)
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

currentPrecRaw <- mean(currentPrecStack, na.rm = T)
currentPrec <- (currentPrecRaw/100)
plot(currentPrec)

#write data to avoid redoing all these time-consuming calculations 
writeRaster(preIndustMeanT, "data/spatialData/climateData/preIndustrialMeanTemp16001900.tif", overwrite = T)
writeRaster(preIndustMaxT , "data/spatialData/climateData/preIndustrialMaxTemp16001900.tif", overwrite = T)
writeRaster(preIndustMinT , "data/spatialData/climateData/preIndustrialMinTemp16001900.tif", overwrite = T)
writeRaster(preIndustPrec, "data/spatialData/climateData/preIndustrialMeanMonthlyPrec16001900.tif", overwrite = T)
writeRaster(currentMeanT , "data/spatialData/climateData/currentMeanTemp20092019.tif", overwrite = T)
writeRaster(currentMaxT , "data/spatialData/climateData/currentMaxTemp20092019.tif", overwrite = T)
writeRaster(currentMinT, "data/spatialData/climateData/currentMinTemp20092019.tif", overwrite = T)
writeRaster(currentPrec , "data/spatialData/climateData/currentMeanMonthlyPrec20092019.tif", overwrite = T)

preIndustMeanT <-  rast("data/spatialData/climateData/preIndustrialMeanTemp16001900.tif")
preIndustMaxT <-  rast("data/spatialData/climateData/preIndustrialMaxTemp16001900.tif")
preIndustMinT <-  rast("data/spatialData/climateData/preIndustrialMinTemp16001900.tif")
preIndustPrec <-  rast("data/spatialData/climateData/preIndustrialMeanMonthlyPrec16001900.tif")
currentMeanT <- rast("data/spatialData/climateData/currentMeanTemp20092019.tif")
currentMaxT <- rast("data/spatialData/climateData/currentMaxTemp20092019.tif")
currentMinT <- rast("data/spatialData/climateData/currentMinTemp20092019.tif")
currentPrec <- rast("data/spatialData/climateData/currentMeanMonthlyPrec20092019.tif")



## Mean Degree Change 
absChangeMeanT <- (currentMeanT-preIndustMeanT)
absChangeMeanT[is.infinite(absChangeMeanT)] <- NA
absChangeMeanT <- clamp(absChangeMeanT,
                        lower=quantile(values(absChangeMeanT), c(.005), na.rm = T),
                        upper=quantile(values(absChangeMeanT), c(.995), na.rm = T),
                        values=TRUE)
plot(absChangeMeanT)
writeRaster(absChangeMeanT, "data/spatialData/climateData/absoluteChangeMeanTemp.tif", overwrite = T)

relChangeMeanT <- (currentMeanT/preIndustMeanT)
relChangeMeanT[is.infinite(relChangeMeanT)] <- NA
relChangeMeanT <- clamp(relChangeMeanT,
                       lower=quantile(values(relChangeMeanT), c(.005), na.rm = T),
                       upper=quantile(values(relChangeMeanT), c(.995), na.rm = T),
                       values=TRUE)
plot(relChangeMeanT)
writeRaster(relChangeMeanT, "data/spatialData/climateData/relativeChangeMeanTemp.tif", overwrite = T)

# Max Degree Change
absChangeMaxT <- (currentMaxT-preIndustMaxT)
absChangeMaxT[is.infinite(absChangeMaxT)] <- NA
absChangeMaxT <- clamp(absChangeMaxT,
                        lower=quantile(values(absChangeMaxT), c(.005), na.rm = T),
                        upper=quantile(values(absChangeMaxT), c(.995), na.rm = T),
                        values=TRUE)
plot(absChangeMaxT)
writeRaster(absChangeMaxT, "data/spatialData/climateData/absoluteChangeMaxTemp.tif", overwrite = T)

relChangeMaxT <- (currentMaxT/preIndustMaxT)
quantile(values(relChangeMaxT))
relChangeMaxT[is.infinite(relChangeMaxT)] <- NA
relChangeMaxT <- clamp(relChangeMaxT,
                       lower=quantile(values(relChangeMaxT), c(.005), na.rm = T),
                       upper=quantile(values(relChangeMaxT), c(.995), na.rm = T),
                       values=TRUE)
plot(relChangeMaxT)
writeRaster(relChangeMaxT, "data/spatialData/climateData/relativeChangeMaxTemp.tif", overwrite = T)

# Min Degree Change
absChangeMinT <- (currentMinT-preIndustMinT)
absChangeMinT[is.infinite(absChangeMinT)] <- NA
absChangeMinT <- clamp(absChangeMinT,
                       lower=quantile(values(absChangeMinT), c(.005), na.rm = T),
                       upper=quantile(values(absChangeMinT), c(.995), na.rm = T),
                       values=TRUE)
plot(absChangeMinT)
writeRaster(absChangeMinT, "data/spatialData/climateData/absoluteChangeMinTemp.tif", overwrite = T)

relChangeMinT <- (currentMinT/preIndustMinT)
relChangeMinT[is.infinite(relChangeMinT)] <- NA
relChangeMinT <- clamp(relChangeMinT,
                       lower=quantile(values(relChangeMinT), c(.005), na.rm = T),
                       upper=quantile(values(relChangeMinT), c(.995), na.rm = T),
                       values=TRUE)
plot(relChangeMinT)
writeRaster(relChangeMinT, "data/spatialData/climateData/relativeChangeMinTemp.tif", overwrite = T)


# Prec. change
plot(currentPrec)
plot(preIndustPrec)
absChangePrec <- (currentPrec-preIndustPrec)
absChangePrec[is.infinite(absChangePrec)] <- NA
absChangePrec <- clamp(absChangePrec,
                       lower=quantile(values(absChangePrec), c(.005), na.rm = T),
                       upper=quantile(values(absChangePrec), c(.995), na.rm = T),
                       values=TRUE)
plot(absChangePrec)
writeRaster(absChangePrec, "data/spatialData/climateData/absoluteChangeMonthlyPrec.tif", overwrite = T)

relChangePrec <- (currentPrec/preIndustPrec)
quantile(values(relChangePrec), c(1, .99, .5, .01, 0), na.rm  = T)
relChangePrec[is.infinite(relChangePrec)] <- NA
relChangePrec <- clamp(relChangePrec,
                       lower=quantile(values(relChangePrec), c(.005), na.rm = T),
                       upper=quantile(values(relChangePrec), c(.995), na.rm = T),
                       values=TRUE)
plot(relChangePrec)
writeRaster(relChangePrec, "data/spatialData/climateData/relativeChangeMonthlyPrec.tif", overwrite = T)

## SLOPES ----------------------

#Define a custom function to fit a linear model and extract slope coefficient, R squared and p value
fitRasterLm <- function(x, allStats = FALSE) {
  # Create an index of the layers e.g. representing months or years (assumes the layers are in order)
  yearI <- 1:length(names(x))
  
  if(sum(is.na(x)) >= (length(names(x))/2)){
    slope = NA
    return((slope))
  }else{
    
    # Fit a linear model (Y = a*X + b) using lm()
    lmFit <- lm(x ~ yearI)
    
    # Extract stats from the linear model
    if(allStats == TRUE){ # calculate all summary stats
      slope <- summary(lmFit)$coefficients['yearI', 'Estimate']
      pValue <- summary(lmFit)$coefficients['yearI', 'Pr(>|t|)']
      rSquared <- summary(lmFit)$r.squared
      
      # Return a vector containing slope, p-value, and R-squared
      return(c(slope = slope, pValue = pValue, rSquared = rSquared))
    }
    
    else{ # get only slope
      slope <- summary(lmFit)$coefficients['yearI', 'Estimate']
      return(slope)
    }
  }
}

terraOptions(memfrac = .9)
terraOptions()

## Temperature slope last 20 years

currentMeanTPathsRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tas/",  full.names = T)
currentMeanTFilesRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tas/",  full.names = F)
currentMeanPaths <- data.frame(paths = currentMeanTPathsRaw, 
                               files = currentMeanTFilesRaw) %>% 
  mutate(year = gsub("CHELSA_tas_", "", files),
         year = gsub("_V.2.1.tif", "", year), 
         year = str_split_i(year, "_", 2) ) %>%
  filter(year %in% c(1980:2019))

listMeanT <- list()

for(Year in unique(currentMeanPaths$year)){
  
  yearPaths <- currentMeanPaths %>%
    filter(year %in% c(Year)) %>% 
    dplyr::select(paths) %>%
    pull()
  
  annualStack <- rast(yearPaths)
  
  annualMeanRaw <- mean(annualStack, na.rm = T)
  annualMean <- ((annualMeanRaw/10)-273.15)
  names(annualMean) <- paste0("year", Year)
  
  listMeanT[[length(listMeanT) + 1]] <- annualMean
  
  print(paste0(Year, " done"))
  
}


#stackMeanT <- rast(listMeanT)
stackMeanT <- do.call(c, listMeanT)


slopeMeanT <- app(stackMeanT, fun=fitRasterLm, allStats = FALSE, cores = 100, overwrite = TRUE,
                  filename = "data/spatialData/climateData/slopeMeanTempChelsa19802019.tif")

plot(slopeMeanT)


## max Temp
currentMaxTPathsRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tasmax/",  full.names = T)
currentMaxTFilesRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tasmax/",  full.names = F)
currentMaxPaths <- data.frame(paths = currentMaxTPathsRaw, 
                              files = currentMaxTFilesRaw) %>% 
  mutate(year = gsub("CHELSA_tasmax_", "", files),
         year = gsub("_V.2.1.tif", "", year), 
         year = str_split_i(year, "_", 2) ) %>%
  filter(year %in% c(1980:2019))

listMaxT <- list()

for(Year in unique(currentMaxPaths$year)){
  
  yearPaths <- currentMaxPaths %>%
    filter(year %in% c(Year)) %>% 
    dplyr::select(paths) %>%
    pull()
  
  annualStack <- rast(yearPaths)
  
  annualMaxRaw <- max(annualStack, na.rm = T)
  
  names(annualMaxRaw) <- paste0("year", Year)
  annualMax <- ((annualMaxRaw/10)-273.15)
  
  
  listMaxT[[length(listMaxT) + 1]] <- annualMax
  
  
  print(paste0(Year, " done"))
  
}

#stackMaxT <- rast(listMaxT)
stackMaxT <- do.call(c, listMaxT)


slopeMaxT <- app(stackMaxT, fun=fitRasterLm, allStats = FALSE, cores = 100, overwrite = TRUE,
                  filename = "data/spatialData/climateData/slopeMaxTempChelsa19802019.tif")

plot(slopeMaxT)

## temp min slope last 20 years
currentMinTPathsRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tasmin/",  full.names = T)
currentMinTFilesRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_tasmin/",  full.names = F)
currentMinPaths <- data.frame(paths = currentMinTPathsRaw, 
                              files = currentMinTFilesRaw) %>% 
  mutate(year = gsub("CHELSA_tasmin_", "", files),
         year = gsub("_V.2.1.tif", "", year), 
         year = str_split_i(year, "_", 2) ) %>%
  filter(year %in% c(1980:2019))

listMinT <- list()

for(Year in unique(currentMinPaths$year)){
  
  yearPaths <- currentMinPaths %>%
    filter(year %in% c(Year)) %>% 
    dplyr::select(paths) %>%
    pull()
  
  annualStack <- rast(yearPaths)
  
  annualMinRaw <- min(annualStack, na.rm = T)
  annualMin <- ((annualMinRaw/10)-273.15)
  
  names(annualMin) <- paste0("year", Year)
  
  listMinT[[length(listMinT) + 1]] <- annualMin
  
  print(paste0(Year, " done"))
  
  
}

stackMinT <- do.call(c, listMinT)


slopeMinT <- app(stackMinT, fun=fitRasterLm, allStats = FALSE, cores = 100, overwrite = TRUE,
                 filename = "data/spatialData/climateData/slopeMinTempChelsa19802019.tif")

plot(slopeMinT)

## prec slope
currentPrecPathsRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_pr/",  full.names = T)
currentPrecFilesRaw <- list.files(path = "O:/Nat_Ecoinformatics/C_Write/_User/MatthewKerr_au738027/Novelty Exposure Mapping/Data/Climate/Modern/CHELSA_pr/",  full.names = F)
currentPrecPaths <- data.frame(paths = currentPrecPathsRaw, 
                               files = currentPrecFilesRaw) %>% 
  mutate(year = gsub("CHELSA_pr_", "", files),
         year = gsub("_V.2.1.tif", "", year), 
         year = str_split_i(year, "_", 2) ) %>%
  filter(year %in% c(1980:2019))

listPrec <- list()

for(Year in unique(currentPrecPaths$year)){
  
  yearPaths <- currentPrecPaths %>%
    filter(year %in% c(Year)) %>% 
    dplyr::select(paths) %>%
    pull()
  
  annualStack <- rast(yearPaths)
  
  annualPrecRaw <- mean(annualStack, na.rm = T)
  annualMPrec <- (annualPrecRaw/100)
  names(annualMPrec) <- paste0("year", Year)
  
  listPrec[[length(listPrec) + 1]] <- annualMPrec
  
  print(paste0(Year, " done"))
  
  
}

stackPrec <- do.call(c, listPrec)


slopePrec <- app(stackPrec, fun=fitRasterLm, allStats = FALSE, cores = 100, overwrite = TRUE,
                 filename = "data/spatialData/climateData/slopePrecChelsa19802019.tif")

plot(slopePrec)

quantile(values(slopePrec), na.rm = T)
