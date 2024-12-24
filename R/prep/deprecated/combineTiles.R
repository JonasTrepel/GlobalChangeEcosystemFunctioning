library(terra)


## Mask according to R2

### Mean Annual Temperature 

mat <- rast("data/spatialData/rawTiles/MatTrend19502023.tif")
plot(mat)

matR2 <-  rast("data/spatialData/rawTiles/MatRsq19502023.tif")
plot(matR2)

MatMask <- (matR2 > 0.2)*1
plot(MatMask)

MatMasked <- mask(mat, MatMask, maskvalues = 0, updatevalue = 0)

plot(MatMasked)

writeRaster(MatMasked,  "data/spatialData/trendData/MatTrend19502023Masked.tif", 
            overwrite = TRUE)

### Max Annual Temperature 

MaxTemp <- rast("data/spatialData/rawTiles/MaxTempTrend19502023.tif")
plot(MaxTemp)

MaxTempR2 <-  rast("data/spatialData/rawTiles/MaxTempRsq19502023.tif")
plot(MaxTempR2)

MaxTempMask <- (MaxTempR2 > 0.2)*1
plot(MaxTempMask)

MaxTempMasked <- mask(MaxTemp, MaxTempMask, maskvalues = 0, updatevalue = 0)

plot(MaxTempMasked)

writeRaster(MaxTempMasked,  "data/spatialData/trendData/MaxTempTrend19502023Masked.tif", 
            overwrite = TRUE)


### Mean Annual Precipitation 


Map <- rast("data/spatialData/rawTiles/MapTrend19502023.tif")
plot(Map)

MapR2 <-  rast("data/spatialData/rawTiles/MapRsq19502023.tif")
plot(MapR2)

MapMask <- (MapR2 > 0.2)*1
plot(MapMask)

MapMasked <- mask(Map, MapMask, maskvalues = 0, updatevalue = 0)

plot(MapMasked)

writeRaster(MapMasked,  "data/spatialData/trendData/MapTrend19502023Masked.tif", 
            overwrite = TRUE)


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
                filename = "data/spatialData/trendData/EviTrend20012023.tif", 
                overwrite = TRUE)
plot(globalEVI)


# globalEVI <- clamp(globalEVI,
#                        lower=-200,
#                        upper=200,
#                        values=FALSE)

# plot(globalEVI)

# writeRaster(globalEVI,  "data/spatialData/trendData/eviTrend20012023.tif", 
#             overwrite = TRUE)

#### EVI Trend R2
eviTRsqFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "EviRsq", full.names = T)

eviTRsqR1 <- rast(eviTRsqFiles[1])
eviTRsqR2 <- rast(eviTRsqFiles[2])
eviTRsqR3 <- rast(eviTRsqFiles[3])
eviTRsqR4 <- rast(eviTRsqFiles[4])
eviTRsqR5 <- rast(eviTRsqFiles[5])
eviTRsqR6 <- rast(eviTRsqFiles[6])
eviTRsqR7 <- rast(eviTRsqFiles[7])
eviTRsqR8 <- rast(eviTRsqFiles[8])




globalEviTRsq <- merge(eviTRsqR1, eviTRsqR2, eviTRsqR3, eviTRsqR4, eviTRsqR5, eviTRsqR6, eviTRsqR7, eviTRsqR8,
                   filename = "data/spatialData/trendData/EviTRsq20012023.tif", 
                   overwrite = TRUE)

#globalEviTRsq[globalEviTRsq < 0] <- NA
plot(globalEviTRsq)

# Mask EVO Trend 
EviMask <- (globalEviTRsq > 0.2)*1
plot(EviMask)

EviMasked <- mask(globalEVI, EviMask, maskvalues = 0, updatevalue = 0)

plot(EviMasked)
writeRaster(EviMasked,  "data/spatialData/trendData/EviTrend20012023Masked.tif", 
            overwrite = TRUE)



#### EVI Mean ----------------------------------------
eviMeanFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "EviMean", full.names = T)

eviMeanR1 <- rast(eviMeanFiles[1])
eviMeanR2 <- rast(eviMeanFiles[2])
eviMeanR3 <- rast(eviMeanFiles[3])
eviMeanR4 <- rast(eviMeanFiles[4])
eviMeanR5 <- rast(eviMeanFiles[5])
eviMeanR6 <- rast(eviMeanFiles[6])
eviMeanR7 <- rast(eviMeanFiles[7])
eviMeanR8 <- rast(eviMeanFiles[8])



globaleviMean <- merge(eviMeanR1, eviMeanR2, eviMeanR3, eviMeanR4, eviMeanR5, eviMeanR6, eviMeanR7, eviMeanR8, 
                    filename = "data/spatialData/otherCovariates/EviMean20012023.tif", 
                   overwrite = TRUE)
plot(globaleviMean)


#### Start of Season TREND ----------------------------------------
sosFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "SosTrend", full.names = T)

sosR1 <- rast(sosFiles[1])
sosR2 <- rast(sosFiles[2])
sosR3 <- rast(sosFiles[3])
sosR4 <- rast(sosFiles[4])
sosR5 <- rast(sosFiles[5])
sosR6 <- rast(sosFiles[6])
sosR7 <- rast(sosFiles[7])
sosR8 <- rast(sosFiles[8])

globalSosTrend <- merge(sosR1, sosR2, sosR3, sosR4, sosR5, sosR6, sosR7, sosR8, 
                         filename = "data/spatialData/trendData/SosTrend20012023.tif", 
                         overwrite = TRUE)

#### Sos Trend R2
SosTRsqFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "SosRsq", full.names = T)

SosTRsqR1 <- rast(SosTRsqFiles[1])
SosTRsqR2 <- rast(SosTRsqFiles[2])
SosTRsqR3 <- rast(SosTRsqFiles[3])
SosTRsqR4 <- rast(SosTRsqFiles[4])
SosTRsqR5 <- rast(SosTRsqFiles[5])
SosTRsqR6 <- rast(SosTRsqFiles[6])
SosTRsqR7 <- rast(SosTRsqFiles[7])
SosTRsqR8 <- rast(SosTRsqFiles[8])




globalSosTRsq <- merge(SosTRsqR1, SosTRsqR2, SosTRsqR3, SosTRsqR4, SosTRsqR5, SosTRsqR6, SosTRsqR7, SosTRsqR8,
                       filename = "data/spatialData/trendData/SosTRsq20012023.tif", 
                       overwrite = TRUE)

plot(globalSosTRsq)

# Mask Sos Trend 
SosMask <- (globalSosTRsq > 0.2)*1
plot(SosMask)

SosMasked <- mask(globalSosTrend, SosMask, maskvalues = 0, updatevalue = 0)

plot(SosMasked)
writeRaster(SosMasked,  "data/spatialData/trendData/SosTrend20012023Masked.tif", 
            overwrite = TRUE)


#### Sos Mean ----------------------------------------
SosMeanFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "SosMean", full.names = T)

SosMeanR1 <- rast(SosMeanFiles[1])
SosMeanR2 <- rast(SosMeanFiles[2])
SosMeanR3 <- rast(SosMeanFiles[3])
SosMeanR4 <- rast(SosMeanFiles[4])
SosMeanR5 <- rast(SosMeanFiles[5])
SosMeanR6 <- rast(SosMeanFiles[6])
SosMeanR7 <- rast(SosMeanFiles[7])
SosMeanR8 <- rast(SosMeanFiles[8])



globalSosMean <- merge(SosMeanR1, SosMeanR2, SosMeanR3, SosMeanR4, SosMeanR5, SosMeanR6, SosMeanR7, SosMeanR8, 
                       filename = "data/spatialData/otherCovariates/SosMean20012023.tif", 
                       overwrite = TRUE)
plot(globalSosMean)



#### Burned Area TREND ----------------------------------------
BurnedAreaFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "BurnedAreaTrend", full.names = T)

BurnedAreaR1 <- rast(BurnedAreaFiles[1])
BurnedAreaR2 <- rast(BurnedAreaFiles[2])
BurnedAreaR3 <- rast(BurnedAreaFiles[3])
BurnedAreaR4 <- rast(BurnedAreaFiles[4])
BurnedAreaR5 <- rast(BurnedAreaFiles[5])
BurnedAreaR6 <- rast(BurnedAreaFiles[6])
BurnedAreaR7 <- rast(BurnedAreaFiles[7])
BurnedAreaR8 <- rast(BurnedAreaFiles[8])

globalBurnedAreaTrend <- merge(BurnedAreaR1, BurnedAreaR2, BurnedAreaR3, BurnedAreaR4, BurnedAreaR5, BurnedAreaR6, BurnedAreaR7, BurnedAreaR8, 
                        filename = "data/spatialData/trendData/BurnedAreaTrend20012023.tif", 
                        overwrite = TRUE)

#### BurnedArea Trend R2
BurnedAreaTRsqFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "BurnedAreaRsq", full.names = T)

BurnedAreaTRsqR1 <- rast(BurnedAreaTRsqFiles[1])
BurnedAreaTRsqR2 <- rast(BurnedAreaTRsqFiles[2])
BurnedAreaTRsqR3 <- rast(BurnedAreaTRsqFiles[3])
BurnedAreaTRsqR4 <- rast(BurnedAreaTRsqFiles[4])
BurnedAreaTRsqR5 <- rast(BurnedAreaTRsqFiles[5])
BurnedAreaTRsqR6 <- rast(BurnedAreaTRsqFiles[6])
BurnedAreaTRsqR7 <- rast(BurnedAreaTRsqFiles[7])
BurnedAreaTRsqR8 <- rast(BurnedAreaTRsqFiles[8])




globalBurnedAreaTRsq <- merge(BurnedAreaTRsqR1, BurnedAreaTRsqR2, BurnedAreaTRsqR3, BurnedAreaTRsqR4, BurnedAreaTRsqR5, BurnedAreaTRsqR6, BurnedAreaTRsqR7, BurnedAreaTRsqR8,
                       filename = "data/spatialData/trendData/BurnedAreaTRsq20012023.tif", 
                       overwrite = TRUE)

plot(globalBurnedAreaTRsq)

# Mask BurnedArea Trend 
BurnedAreaMask <- (globalBurnedAreaTRsq > 0.2)*1
plot(BurnedAreaMask)

BurnedAreaMasked <- mask(globalBurnedAreaTrend, BurnedAreaMask, maskvalues = 0, updatevalue = 0)

plot(BurnedAreaMasked)
writeRaster(BurnedAreaMasked,  "data/spatialData/trendData/BurnedAreaTrend20012023Masked.tif", 
            overwrite = TRUE)


#### BurnedArea Mean ----------------------------------------
BurnedAreaMeanFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "BurnedAreaMean", full.names = T)

BurnedAreaMeanR1 <- rast(BurnedAreaMeanFiles[1])
BurnedAreaMeanR2 <- rast(BurnedAreaMeanFiles[2])
BurnedAreaMeanR3 <- rast(BurnedAreaMeanFiles[3])
BurnedAreaMeanR4 <- rast(BurnedAreaMeanFiles[4])
BurnedAreaMeanR5 <- rast(BurnedAreaMeanFiles[5])
BurnedAreaMeanR6 <- rast(BurnedAreaMeanFiles[6])
BurnedAreaMeanR7 <- rast(BurnedAreaMeanFiles[7])
BurnedAreaMeanR8 <- rast(BurnedAreaMeanFiles[8])



globalBurnedAreaMean <- merge(BurnedAreaMeanR1, BurnedAreaMeanR2, BurnedAreaMeanR3, BurnedAreaMeanR4, BurnedAreaMeanR5, BurnedAreaMeanR6, BurnedAreaMeanR7, BurnedAreaMeanR8, 
                       filename = "data/spatialData/otherCovariates/BurnedAreaMean20012023.tif", 
                       overwrite = TRUE)
plot(globalBurnedAreaMean)



#### Landcover ----------------------------------------
lcFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "GlobalLand", full.names = T)

lcR1 <- rast(lcFiles[1])
lcR2 <- rast(lcFiles[2])


globalLC <- merge(lcR1, lcR2, 
                  filename = "data/spatialData/otherCovariates/GlobalLandCoverCopernicus2019.tif", 
                  overwrite = TRUE)


plot(globalLC)
