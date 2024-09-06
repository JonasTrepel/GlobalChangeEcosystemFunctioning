
############ BUILD CLIMATE CHANGE METRICS ###############
## Using Chelsa and Paleoclim Data 
## 
library(terra)

paleoMeanTFiles <- c("O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/BIOCLIM/CHELSA_TraCE21k_bio01_16_V1.0.tif",
                     "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/BIOCLIM/CHELSA_TraCE21k_bio01_17_V1.0.tif",
                     "O:/Nat_Ecoinformatics/C_Write/_User/AlejoOrdonez_au467796/GLOBAL/PALEOCLIMATE/CHELSA-TraCE21k/Data/BIOCLIM/CHELSA_TraCE21k_bio01_18_V1.0.tif")

paleoMeanT <- rast(paleoMeanTFiles)

preIndustMeanT <- mean(paleoMeanT, na.rm = T)
plot(preIndustMeanT)

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