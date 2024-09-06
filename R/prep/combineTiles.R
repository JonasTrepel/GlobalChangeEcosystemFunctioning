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
quantile(values(globalEVI), na.rm = T, c(.99, 0.95, 0.75, 0.5, 0.25, 0.05, 0.01))
summary(values(globalEVI))
plot(globalEVI)

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
quantile(values(globalFireTrend))

