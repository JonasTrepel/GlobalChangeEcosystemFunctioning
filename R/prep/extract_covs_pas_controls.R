
library(terra)
library(tidyverse)
library(data.table)
library(sf)
library(tidylog)
library(sf)


### define if we want to run it for control or PA 

param <- "pa_grid"
#param <- "pa"
#param <- "control"

if(param == "pa"){
  pas <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg")
} else if(param == "control"){
  pas <- read_sf("data/spatialData/protectedAreas/controls_for_strict_pas.gpkg")
} else if(param == "pa_grid") {
  pas <- read_sf("data/spatialData/protectedAreas/pa_and_control_grid_1km.gpkg")
}

# get legends -----------

## Land cover 

# extract land cover and remove urban and cultivated pixels
land_cover_num <-  c(0, 20, 30, 40, 50, 60, 70, 80, 90, 
            100, 111, 112, 113, 114, 115, 116,
            121, 122, 123, 124, 125, 126, 
            200)
land_cover <- c("Unknown", 
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

land_cover_leg <- data.frame(land_cover_num = land_cover_num, land_cover = land_cover)

wwfBiome <- read_sf("data/spatialData/otherCovariates/WWF_BIOMES.gpkg")

biome_leg <- wwfBiome %>%
  as.data.table() %>%
  mutate(geom = NULL) %>%
  unique() %>% 
  rename(olson_biome_num = BIOME, 
         olson_biome = BIOME_Name)

clim_reg_leg <- data.frame(climatic_region_num = 1:30,
                           climatic_region = c(
                             rep("Tropical", 3),
                             rep("Arid", 4),
                             rep("Temperate", 9),
                             rep("Cold", 12),
                             rep("Polar", 2) ))

fun_biome_leg <- data.table(
  functional_biome_num = c(1:24), 
  functional_biome = c("SLC","SMC","SHC","TLC","TMC","THC",
                       "SLD","SMD","SHD","TLD","TMD","THD",
                       "SLB","SMB","SHB","TLB","TMB","THB","SLN","SMN","SHN",
                       "TLN","TMN","THN"))

#### EXTRACT CONTINUOUS COVS #### ----------------------------------
colNames <- c(
  #### Environmental ####
  "elevation", ## Elevation
  "MAP", ## MAP
  "MAT", ## MAT

  #### Global Change ####
  "nitrogen_depo", ## Nitrogen depo
  "human_modification", #Human modification index
  
  #### categorical ####
  "land_cover_num", # land cover
  "functional_biome_num", # functional biome 
  "olson_biome_num", # olson biome 
  "climatic_region_num" # climatic region 
)

covPaths <- c(
  #### Environmental ####
  "data/spatialData/otherCovariates/Elevation_Global_930m.tif", ## Elevation
  "data/spatialData/climateData/currentMeanMonthlyPrec20092019.tif", ## MAP
  "data/spatialData/climateData/currentMeanTemp20092019.tif", ##MAT

  #### Global Change ####
  "data/spatialData/otherCovariates/total_N_dep.tif",## Nitrogen depo
  "data/spatialData/otherCovariates/lulc-human-modification-terrestrial-systems_geographic.tif", #Human modification index
  
  #### categorical ####
  "data/spatialData/otherCovariates/GlobalLandCovercopernicus2019.tif", #land cover
  "data/spatialData/otherCovariates/higginsFunctionalBiomes.tif", # functional biome 
  "data/spatialData/otherCovariates/wwf_olson_biome.tif", # olson biome 
  "data/spatialData/otherCovariates/KoppenGeigerClimaticRegions1km.tif"# climatic region 
)


funcs <- c(rep("mean", 5), rep("mode", 4))

covs <- data.table(
  colName = colNames, 
  covPath = covPaths, 
  func = funcs
) %>% filter(!is.na(covPaths))


pas_covs_raw <- pas %>% as.data.table() %>% mutate(geom = NULL, x = NULL, geometry = NULL)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

# Create and register a cluster
clust <- makeCluster(10)
registerDoSNOW(clust)

## progress bar 
iterations <- nrow(covs)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

##############################################################################            
################################## LOOOOOOOOOOOOP ############################            
##############################################################################    

tic()

dtCovs <- foreach(i = 1:nrow(covs),
                  .packages = c('tidyverse', 'exactextractr', 'data.table', 'terra', 'sf'),
                  .options.snow = opts,
                  .inorder = FALSE,
                  .verbose = TRUE, 
                  .combine = left_join) %dopar% {
                    
                    #for(i in 1:nrow(covs)){
                    
                    covR <- rast(covs[i, ]$covPath)
                    
                    world_grid_t <- st_transform(pas, crs = st_crs(covR))
                    
                    all_polys <- pas %>% dplyr::select(unique_id)
                    
                    func <- covs[i, ]$func
                    
                    
                    extr <- exactextractr::exact_extract(covR, 
                                                         world_grid_t, 
                                                         append_cols = c("unique_id"),
                                                         fun = func)
                    
                    setnames(extr, func, covs[i, ]$colName)
                    
                    extraDTFin <- all_polys %>%
                      left_join(extr) %>%
                      as.data.table() %>%
                      mutate(geom = NULL, x = NULL, geometry = NULL) %>% 
                      unique()
                    
                    return(extraDTFin)
                    
                  }

stopCluster(clust)
toc()


dt_pas_covs <- pas_covs_raw %>% 
  left_join(dtCovs) %>% 
  left_join(clim_reg_leg) %>% 
  left_join(biome_leg) %>%
  left_join(fun_biome_leg) %>% 
  left_join(land_cover_leg) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL) %>% 
  dplyr::select(-land_cover_num, -functional_biome_num, -olson_biome_num, -climatic_region_num)
summary(dt_pas_covs)
pa.shapes <- dt_pas_covs %>% left_join(pas) %>% st_as_sf


if(param == "pa"){
  fwrite(dt_pas_covs, "data/processedData/cleanData/pas_with_covs.csv")
} else if(param == "control"){
  fwrite(dt_pas_covs, "data/processedData/cleanData/controls_with_covs.csv")
} else if(param == "pa_grid") {
  fwrite(dt_pas_covs, "data/processedData/cleanData/pa_and_control_grid_with_covs.csv")
}
