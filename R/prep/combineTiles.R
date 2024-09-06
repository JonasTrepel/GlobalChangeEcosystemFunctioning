library(terra)


### Combine the different tiles from the Google Earth Engine Output ###


#### EVI TREND ----------------------------------------
eviFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "EviTrend", full.names = T)

eviR1 <- rast(eviFiles[1])
eviR2 <- rast(eviFiles[2])
eviR3 <- rast(eviFiles[3])
eviR4 <- rast(eviFiles[4])
eviR5 <- rast(eviFiles[5])
eviR6 <- rast(eviFiles[6])
eviR7 <- rast(eviFiles[7])
eviR8 <- rast(eviFiles[8])



globalEVI <- merge(eviR1, eviR2, eviR3, eviR4, eviR5, eviR6, eviR7, eviR8, 
                filename = "data/spatialData/trendData/eviTrend20032023.tif", 
                overwrite = TRUE)
globalEVI[globalEVI > 100] <- NA
globalEVI[globalEVI < -100] <- NA

plot(globalEVI)

writeRaster(globalEVI,  "data/spatialData/trendData/eviTrend20032023.tif", 
            overwrite = TRUE)

#### Fire Frequency TREND ----------------------------------------
fireFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "FireFreq", full.names = T)

fireR1 <- rast(fireFiles[1])
fireR2 <- rast(fireFiles[2])
fireR3 <- rast(fireFiles[3])
fireR4 <- rast(fireFiles[4])
fireR5 <- rast(fireFiles[5])
fireR6 <- rast(fireFiles[6])
fireR7 <- rast(fireFiles[7])
fireR8 <- rast(fireFiles[8])



globalFireTrend <- merge(fireR1, fireR2, fireR3, fireR4, fireR5, fireR6, fireR7, fireR8, 
                   filename = "data/spatialData/trendData/fireFreqTrend20032023.tif", 
                   overwrite = TRUE)

plot(globalFireTrend)

#### Start of Season TREND ----------------------------------------
sosFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "sosTrend", full.names = T)

sosR1 <- rast(sosFiles[1])
sosR2 <- rast(sosFiles[2])
sosR3 <- rast(sosFiles[3])
sosR4 <- rast(sosFiles[4])
sosR5 <- rast(sosFiles[5])
sosR6 <- rast(sosFiles[6])
sosR7 <- rast(sosFiles[7])
sosR8 <- rast(sosFiles[8])



globalSosTrend <- merge(sosR1, sosR2, sosR3, sosR4, sosR5, sosR6, sosR7, sosR8, 
                         filename = "data/spatialData/trendData/sosTrend20032023.tif", 
                         overwrite = TRUE)

plot(globalSosTrend)

#### End of Season TREND ----------------------------------------
eosFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "eosTrend", full.names = T)

eosR1 <- rast(eosFiles[1])
eosR2 <- rast(eosFiles[2])
eosR3 <- rast(eosFiles[3])
eosR4 <- rast(eosFiles[4])
eosR5 <- rast(eosFiles[5])
eosR6 <- rast(eosFiles[6])
eosR7 <- rast(eosFiles[7])
eosR8 <- rast(eosFiles[8])



globalEosTrend <- merge(eosR1, eosR2, eosR3, eosR4, eosR5, eosR6, eosR7, eosR8, 
                        filename = "data/spatialData/trendData/eosTrend20032023.tif", 
                        overwrite = TRUE)

plot(globalEosTrend)


#### Length of Season TREND ----------------------------------------
losFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "losTrend", full.names = T)

losR1 <- rast(losFiles[1])
losR2 <- rast(losFiles[2])
losR3 <- rast(losFiles[3])
losR4 <- rast(losFiles[4])
losR5 <- rast(losFiles[5])
losR6 <- rast(losFiles[6])
losR7 <- rast(losFiles[7])
losR8 <- rast(losFiles[8])



globalLosTrend <- merge(losR1, losR2, losR3, losR4, losR5, losR6, losR7, losR8, 
                        filename = "data/spatialData/trendData/losTrend20032023.tif", 
                        overwrite = TRUE)

plot(globalLosTrend)
#quantile(values(globalLosTrend))

