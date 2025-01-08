Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


library(tidyverse)
library(sf)
library(rnaturalearth)
library(mapview)
library(terra)
library(data.table)

be_crs <- st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

pas <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg")

pas_t <- pas %>% st_transform(crs = be_crs) %>% 
  mutate(
    iucn_cat = case_when(
      iucnCatOrd == 1 ~ "Ia",  
      iucnCatOrd == 2 ~ "Ib",  
      iucnCatOrd == 3 ~ "II",  
      iucnCatOrd == 4 ~ "III",  
      iucnCatOrd == 5 ~ "IV",  
      iucnCatOrd == 6 ~ "V",  
      iucnCatOrd == 7 ~ "VI",  
      iucnCatOrd == 8 ~ "unknown_or_NA",  
    ))

############# skip this if first grid is already done #########

## Load world and make Grid
world <- rnaturalearth::ne_countries() %>%
  dplyr::select(continent, sovereignt, geometry) %>% 
  st_transform(crs = st_crs(be_crs))

## For now 5km
world_grid <- st_make_grid(world, cellsize = c(5000, 5000)) %>% st_as_sf()

world_grid2 <- world_grid %>%
  st_join(world) %>%
  filter(!is.na(sovereignt)) %>% 
  mutate(gridID = paste0("grid", 1:nrow(.)))

# extract landcover ----------------------------
# and remove urban and cultivated pixels
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
plot(lc)

grid_trans_lc <- world_grid2 %>%
  st_transform(crs(lc))

lcExtr <- exactextractr::exact_extract(lc,
                                       grid_trans_lc,
                                       summarize_df = TRUE,
                                       fun = function(df){
                                         dat <- df[!is.na(df$value) & df$coverage_fraction > 0.75, ]
                                         uniqueValues <- unique(dat$value)
                                         mode <- uniqueValues[which.max(tabulate(match(dat$value, uniqueValues)))]
                                         return(mode)
                                       } #https://rdrr.io/cran/exactextractr/man/exact_extract.html, see under User-defined summary functions
)


lcExtrFin <- data.table(lcMode = lcExtr) %>% 
  cbind(grid_trans_lc[, "gridID"]) %>% 
  rename(lcNum = `lcMode`) %>%
  mutate(x = NULL) %>% 
  left_join(lcLeg, by = "lcNum") %>% 
  dplyr::select(gridID, LandCover) 

world_grid3 <- world_grid2 %>%
  left_join(lcExtrFin) %>% 
  filter(!LandCover %in% c("Unknown", "Cultivated", "Urban", "Ocean", "Water")) 

unique(world_grid3$LandCover)

#get grid coordinates
coords <- world_grid3 %>% st_centroid() %>% st_coordinates()

coords_lat_lon <- world_grid3 %>% st_transform(4326) %>% st_centroid() %>% st_coordinates()


world_grid3$X <- coords[,1]
world_grid3$Y <- coords[,2]

world_grid3$lon <- coords_lat_lon[,1]
world_grid3$lat <- coords_lat_lon[,2]

write_sf(world_grid3, "data/spatialData/protectedAreas/paGrid.gpkg", append = FALSE)


## Biome raster ---------------

wwfBiome <- read_sf("data/spatialData/otherCovariates/WWF_BIOMES.gpkg")


rastTmp <- terra::rast(res = 0.01) #roughly 1 km at equator 
biomeR <- terra::rasterize(wwfBiome, rastTmp, field = "BIOME")
plot(biomeR)


writeRaster(biomeR, "data/spatialData/otherCovariates/wwf_olson_biome.tif")


####### continue here #####
# get legends -----------
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

pa_meta <- pas_t %>%
  as.data.table() %>%
  mutate(geom = NULL) %>%
  dplyr::select(own_pa_id, WDPA_PID, NAME, STATUS_YR, area_km2, GIS_AREA)


#### EXTRACT CONTINUOUS COVS #### ----------------------------------
world_grid3 <- read_sf("data/spatialData/protectedAreas/paGrid.gpkg")

#grid_coords <- world_grid3 %>% as.data.table() %>% mutate(geom = NULL) %>% dplyr::select(gridID, X, Y, lon, lat)


colNames <- c(
  #### Environmental ####
  "elevation", ## Elevation
  "MAP", ## MAP
  "MAT", ## MAT
  "MaxTemp", ## MaxTemp
  "BurnedAreaMean", ## Burned Area
  "EVI", ## EVI
  
  #### Global Change ####
  "nitrogen_depo", ## Nitrogen depo
  "human_modification", #Human modification index
  
  ### PA fract
  "fract_prot", # fraction protected
  
  #### categorical ####
  "iucn_cat_ord",# PA 
  "own_pa_id", #own pa id
  "functional_biome_num", # functional biome 
  "olson_biome_num", # olson biome 
  "climatic_region_num" # climatic region 
  
)

covPaths <- c(
  #### Environmental ####
  "data/spatialData/otherCovariates/Elevation_Global_930m.tif", ## Elevation
  "data/spatialData/climateData/currentMeanMonthlyPrec20092019.tif", ## MAP
  "data/spatialData/climateData/currentMeanTemp20092019.tif", ##MAT
  "data/spatialData/climateData/currentMaxTemp20092019.tif", ##Max Temp
  "data/spatialData/otherCovariates/BurnedAreaMean20012023.tif", ## Burned Area
  "data/spatialData/otherCovariates/EviMean20012023.tif", ## EVI
  
  #### Global Change ####
#  "data/spatialData/trendData/MapTrend19502023Masked.tif", ## Prec slope masked
  "data/spatialData/otherCovariates/total_N_dep.tif",## Nitrogen depo
  "data/spatialData/otherCovariates/lulc-human-modification-terrestrial-systems_geographic.tif", #Human modification index
 
  #### fraction protected ####
  "data/spatialData/protectedAreas/pa_raster_binary.tif",


  #### categorical ####
 "data/spatialData/protectedAreas/pa_raster.tif", # PA 
 "data/spatialData/protectedAreas/pa_id_raster.tif", #pwn pa id
 "data/spatialData/otherCovariates/higginsFunctionalBiomes.tif", # functional biome 
 "data/spatialData/otherCovariates/wwf_olson_biome.tif", # olson biome 
 "data/spatialData/otherCovariates/KoppenGeigerClimaticRegions1km.tif"# climatic region 

)


funcs <- c(rep("mean", 9), rep("mode", 5))

covs <- data.table(
  colName = colNames, 
  covPath = covPaths, 
  func = funcs
) %>% filter(!is.na(covPaths))


grid_covs_raw <- world_grid3 %>% as.data.table() %>% mutate(geom = NULL, x = NULL, geometry = NULL)

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

# Create and register a cluster
clust <- makeCluster(14)
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
                    
                    world_grid_t <- st_transform(world_grid3, crs = st_crs(covR))
                    
                    all_grids <- grid_covs_raw %>% dplyr::select(gridID)
                    
                    func <- covs[i, ]$func
                    
                    
                    extr <- exactextractr::exact_extract(covR, 
                                                         world_grid_t, 
                                                         append_cols = c("gridID"),
                                                         fun = func)
                   
                    setnames(extr, func, covs[i, ]$colName)
                    
                    extraDTFin <- all_grids %>%
                      left_join(extr) %>%
                      as.data.table() %>%
                      mutate(geom = NULL, x = NULL, geometry = NULL) %>% 
                      unique()
                    
                    return(extraDTFin)
                    
  }

stopCluster(clust)
toc()


grid_covs_raw2 <- grid_covs_raw %>% 
  left_join(dtCovs %>% unique()) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL) %>% 
  left_join(clim_reg_leg) %>% 
  left_join(biome_leg) %>%
  left_join(fun_biome_leg) %>% 
  left_join(pa_meta) %>% 
  #left_join(grid_coords) %>% 
  mutate(
    iucn_cat = case_when(
      iucn_cat_ord == 1 ~ "Ia",  
      iucn_cat_ord == 2 ~ "Ib",  
      iucn_cat_ord == 3 ~ "II",  
      iucn_cat_ord == 4 ~ "III",  
      iucn_cat_ord == 5 ~ "IV",  
      iucn_cat_ord == 6 ~ "V",  
      iucn_cat_ord == 7 ~ "VI",  
      iucn_cat_ord == 8 ~ "unknown_or_NA", 
      is.na(iucn_cat_ord) ~ "unprotected"
    )) %>%
  dplyr::select(-climatic_region_num, -olson_biome_num, -functional_biome_num)

summary(grid_covs_raw2)

fwrite(grid_covs_raw2, "data/processedData/cleanData/gridWithCovs.csv")

gridWithCovs <- world_grid3 %>% 
  left_join(grid_covs_raw2 %>% unique())

write_sf(gridWithCovs, "data/spatialData/protectedAreas/gridWithCovs.gpkg", append = FALSE)

#### get polygon based coverage fraction  


# temp_grid <- world_grid3 %>%
#  mutate(WDPA_PID = "NA", 
#         NAME = "NA", 
#         GIS_AREA = 0, 
#         DESIG_ENG = "no_designation", 
#         fract_prot = 0, 
#         area_km2 = 0,
#         n_pas_overlap = 0) %>% 
#  left_join(grid_covs_raw2) 


# grid_pa <- temp_grid %>% 
#   filter(!is.na(iucn_cat_ord))
# 
# grid_no_pa <- temp_grid %>% 
#   filter(is.na(iucn_cat_ord))
# 
# 
# ############### create cluster ####################
# library(doSNOW)
# library(foreach)
# library(tictoc)
# 
# # Create and register a cluster
# clust <- makeCluster(100)
# registerDoSNOW(clust)
# 
# ## progress bar 
# iterations <- nrow(grid_pa)
# pb <- txtProgressBar(max = iterations, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# 
# ##############################################################################            
# ####################### LOOOOOOOOOOOOP 2          ############################            
# ##############################################################################    
# #dt.tier <- dt.tier[3,]
# tic()
# 
# grid_pa2 <- foreach(i = 1:nrow(grid_pa),
#                     .packages = c('tidyverse', 'exactextractr', 'data.table', 'terra', 'sf'),
#                     .options.snow = opts,
#                     .inorder = FALSE,
#                     .combine = rbind) %dopar% {
#                       
#                       #for(i in 1:nrow(grid_pa)){
#                       
#                       iucn.cat <- grid_pa[i, ]$iucn_cat
#                       
#                       pa_sub <- pas_t %>%
#                         rename(pa_id = WDPA_PID) %>% 
#                         filter(iucn_cat %in% iucn.cat) %>% 
#                         dplyr::select(pa_id) %>% 
#                         unique()
#                       
#                       int <- st_intersection(grid_pa[i, ], pa_sub) 
#                       
#                       if(nrow(int) == 0){
#                         
#                         grid_pa2_row <- grid_pa[i, ]   
#                         
#                       }else{
#                         
#                         int$area <- st_area(int)
#                         int_area <- sum(int$area)
#                         grid_area <- st_area(grid_pa[i, ])
#                         
#                         larg_int_id <- int %>% 
#                           slice_max(area) %>% #select PA with largest overlap
#                           as.data.frame() %>%
#                           dplyr::select(pa_id) %>%
#                           pull()
#                         
#                         pa_larg_int <- pas_t %>% 
#                           filter(WDPA_PID == larg_int_id) %>%
#                           dplyr::select(-unique_id, -Country) %>%
#                           unique()
#                         
#                         grid_pa[i, ]$fract_prot <- as.numeric(int_area/grid_area)
#                         
#                         grid_pa[i, ]$WDPA_PID <- pa_larg_int$WDPA_PID
#                         grid_pa[i, ]$NAME <- pa_larg_int$NAME
#                         grid_pa[i, ]$DESIG_ENG <- pa_larg_int$DESIG_ENG
#                         grid_pa[i, ]$GIS_AREA <- pa_larg_int$GIS_AREA
#                         grid_pa[i, ]$area_km2 <- pa_larg_int$area_km2
#                         grid_pa[i, ]$n_pas_overlap <- nrow(int)
#                         
#                         grid_pa2_row <- grid_pa[i, ]
#                         
#                       }
#                       
#                       
#                       # if(i == 1){
#                       #    grid_pa2 <- grid_pa2_row
#                       #  }else{
#                       #    grid_pa2 <- rbind(grid_pa2_row, grid_pa2)
#                       #  }
#                       
#                       return(grid_pa2_row)
#                       
#                       # print(paste(round((i/nrow(grid_pa))*100, 3)))
#                       
#                     }
# 
# toc()
# 
# grid_covs <- rbind(grid_pa2, grid_no_pa) %>% 
#   mutate()
