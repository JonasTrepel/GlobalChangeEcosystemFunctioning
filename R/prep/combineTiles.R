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
               # filename = "data/spatialData/trendData/eviTrend20032023.tif", 
                overwrite = TRUE)
plot(globalEVI)
globalEVI <- clamp(globalEVI,
                       lower=-100,
                       upper=100,
                       values=FALSE)

plot(globalEVI)

writeRaster(globalEVI,  "data/spatialData/trendData/eviTrend20012023.tif", 
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

#### NPP Trend ----------------------------------------
nppTrendFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "NppModisSimple", full.names = T)

nppTrendR1 <- rast(nppTrendFiles[1])
nppTrendR2 <- rast(nppTrendFiles[2])
nppTrendR3 <- rast(nppTrendFiles[3])
nppTrendR4 <- rast(nppTrendFiles[4])
nppTrendR5 <- rast(nppTrendFiles[5])
nppTrendR6 <- rast(nppTrendFiles[6])
nppTrendR7 <- rast(nppTrendFiles[7])
nppTrendR8 <- rast(nppTrendFiles[8])



globalNppTrend <- merge(nppTrendR1, nppTrendR2, nppTrendR3, nppTrendR4, nppTrendR5, nppTrendR6, nppTrendR7, nppTrendR8, 
                       filename = "data/spatialData/trendData/NppTrend20012023.tif", 
                       overwrite = TRUE)
plot(globalNppTrend)

#### NPP Mean ----------------------------------------
nppMeanFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "NppModisMean", full.names = T)

nppMeanR1 <- rast(nppMeanFiles[1])
nppMeanR2 <- rast(nppMeanFiles[2])
nppMeanR3 <- rast(nppMeanFiles[3])
nppMeanR4 <- rast(nppMeanFiles[4])
nppMeanR5 <- rast(nppMeanFiles[5])
nppMeanR6 <- rast(nppMeanFiles[6])
nppMeanR7 <- rast(nppMeanFiles[7])
nppMeanR8 <- rast(nppMeanFiles[8])



globalNppMean <- merge(nppMeanR1, nppMeanR2, nppMeanR3, nppMeanR4, nppMeanR5, nppMeanR6, nppMeanR7, nppMeanR8, 
                        filename = "data/spatialData/otherCovariates/NppMean20012023.tif", 
                        overwrite = TRUE)
plot(globalNppMean)

#### Landcover ----------------------------------------
lcFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "GlobalLand", full.names = T)

lcR1 <- rast(lcFiles[1])
lcR2 <- rast(lcFiles[2])


globalLC <- merge(lcR1, lcR2, 
                   filename = "data/spatialData/otherCovariates/GlobalLandCoverCopernicus2019.tif", 
                   overwrite = TRUE)


plot(globalLC)


#### Fire Frequency TREND ----------------------------------------
fireFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "FireFreqTrend", full.names = T)

fireR1 <- rast(fireFiles[1])
fireR2 <- rast(fireFiles[2])
fireR3 <- rast(fireFiles[3])
fireR4 <- rast(fireFiles[4])
fireR5 <- rast(fireFiles[5])
fireR6 <- rast(fireFiles[6])
fireR7 <- rast(fireFiles[7])
fireR8 <- rast(fireFiles[8])



globalFireTrend <- merge(fireR1, fireR2, fireR3, fireR4, fireR5, fireR6, fireR7, fireR8, 
                   filename = "data/spatialData/trendData/FireFreqTrend20012023.tif", 
                   overwrite = TRUE)

plot(globalFireTrend)

#### Fire Frequency MEAN ----------------------------------------
fireMeanFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "FireFreqMean", full.names = T)

fireMeanR1 <- rast(fireMeanFiles[1])
fireMeanR2 <- rast(fireMeanFiles[2])
fireMeanR3 <- rast(fireMeanFiles[3])
fireMeanR4 <- rast(fireMeanFiles[4])
fireMeanR5 <- rast(fireMeanFiles[5])
fireMeanR6 <- rast(fireMeanFiles[6])




globalfireMeanTrend <- merge(fireMeanR1, fireMeanR2, fireMeanR3, fireMeanR4, fireMeanR5, fireMeanR6,  
                         filename = "data/spatialData/otherCovariates/FireMeanFreqTrend20012023.tif", 
                         overwrite = TRUE)

plot(globalfireMeanTrend)

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

globalSosTrend <- clamp(globalSosTrend,
                   lower=-50,
                   upper=50,
                   values=FALSE)

plot(globalSosTrend)


writeRaster(globalSosTrend,  "data/spatialData/trendData/SosTrend20012023.tif", 
            overwrite = TRUE)

#### SD (evi; 0.25 at 1km) TREND ----------------------------------------
sdTrendFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "SDTrend", full.names = T)

sdTrendR1 <- rast(sdTrendFiles[1])
sdTrendR2 <- rast(sdTrendFiles[2])
sdTrendR3 <- rast(sdTrendFiles[3])
sdTrendR4 <- rast(sdTrendFiles[4])
sdTrendR5 <- rast(sdTrendFiles[5])
sdTrendR6 <- rast(sdTrendFiles[6])




globalsdTrendTrend <- merge(sdTrendR1, sdTrendR2, sdTrendR3, sdTrendR4, sdTrendR5, sdTrendR6,  
                             filename = "data/spatialData/otherCovariates/sdTrendFreqTrend20012023.tif", 
                             overwrite = TRUE)

plot(globalsdTrendTrend)

#### SD (evi; 0.25 at 1km) TREND ----------------------------------------
sdTrendFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "SDTrend", full.names = T)

sdTrendR1 <- rast(sdTrendFiles[1])
sdTrendR2 <- rast(sdTrendFiles[2])

globalEviSdTrend <- merge(sdTrendR1, sdTrendR2, 
                            filename = "data/spatialData/trendData/EviSdTrend20012023.tif", 
                            overwrite = TRUE)

plot(globalEviSdTrend)

globalEviSdTrend <- clamp(globalEviSdTrend,
                        lower=quantile(values(globalEviSdTrend), .01, na.rm = T),
                        upper=quantile(values(globalEviSdTrend), .99, na.rm = T),
                        values=FALSE)

plot(globalEviSdTrend)


writeRaster(globalEviSdTrend,  "data/spatialData/trendData/EviSdTrend20012023.tif", 
            overwrite = TRUE)

#### SD (evi; 0.25 at 1km) MEAN ----------------------------------------
sdMeanFiles <- list.files(path = "data/spatialData/rawTiles/", pattern = "SDMean", full.names = T)

sdMeanR1 <- rast(sdMeanFiles[1])
sdMeanR2 <- rast(sdMeanFiles[2])

globalEviSdMean <- merge(sdMeanR1, sdMeanR2,   
                            filename = "data/spatialData/otherCovariates/EviSdMean20012023.tif", 
                            overwrite = TRUE)

plot(globalEviSdMean)
