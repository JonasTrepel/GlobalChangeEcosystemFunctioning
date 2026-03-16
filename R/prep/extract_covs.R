### extract covariates 
#library(terra)
library(tidyverse)
library(data.table)
library(sf)
library(tidylog)
library(sf)
library(tictoc)
library(furrr)
library(terra)
library(exactextractr)
### define if we want to run it for control or PA 

#param <- "world_grid"
param = "pas"
#param = "pa_controls"

if(param == "world_grid"){
 
   vect <- read_sf("data/spatial_data/grids/world_grid.gpkg")
  
} else if(param == "pas"){
  
  vect <- read_sf("data/spatial_data/protected_areas/pa_shapes.gpkg")
  
} else if(param == "pa_controls"){
  
  vect <- read_sf("data/spatial_data/protected_areas/controls_for_pas.gpkg")
  
} 

# get legends -----------
wwf_biome <- read_sf("data/spatial_data/covariates/WWF_BIOMES.gpkg")

biome_leg <- wwf_biome %>%
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


lc_num <-  c(0, 20, 30, 40, 50, 60, 70, 80, 90, 
             100, 111, 112, 113, 114, 115, 116,
             121, 122, 123, 124, 125, 126, 
             200)
land_cover <- c("unknown", 
                "shrubs", 
                "herbacous_veg", 
                "cultivated",
                "urban",
                "bare",
                "snow_ice",
                "water",
                "herbaceous_wetland",
                "moss_and_lichen",
                "closed_forest_evergreen_needle",
                "closed_forest_evergreen_broad",
                "closed_forest_deciduous_needle",
                "closed_forest_deciduous_broad",
                "closed_forest_mixed",
                "closed_forest_other",
                "open_forest_evergreen_needle",
                "open_forest_evergreen_broad",
                "open_forest_deciduous_needle",
                "open_forest_deciduous_broad",
                "open_forest_mixed",
                "open_forest_other",
                "ocean")

lc_leg <- data.frame(land_cover=land_cover, lc_num=lc_num)

#get other file paths sorted -------------
col_names <- c(
  
  #### Environmental ####
  "evi", 
  "max_temp", 
  "mat", 
  "map", 
  "burned_area", 
  "fire_frequency", 
  "elevation", ## Elevation
  "nitrogen_depo", ## Nitrogen depo
  "hmi_2020", #Human modification index
  "hmi_change", #Human modification index
#  "land_cover_exclusion", #if values > 0, it should be excluded

  ### PA fract
  "fract_prot", # fraction protected
  
  #### categorical ####
  "iucn_cat_ord",# PA 
  "functional_biome_num", # functional biome 
  "olson_biome_num", # olson biome 
  "climatic_region_num", # climatic region 
  "lc_num"

  )



cov_paths <- c(
  
  "data/spatial_data/covariates/modis_evi.tif", 
  "data/spatial_data/covariates/era5_max_temp.tif", 
  "data/spatial_data/covariates/era5_mat.tif", 
  "data/spatial_data/covariates/era5_prec.tif", 
  "data/spatial_data/covariates/modis_burned_area.tif", 
  "data/spatial_data/covariates/n_fires_500m_2001_2024.tif", 
  "data/spatial_data/covariates/Elevation_Global_930m.tif", ## Elevation
  "data/spatial_data/covariates/total_N_dep.tif", ## Nitrogen depo
  "data/spatial_data/covariates/HMv20240801_2020c_AA.tif", #Human modification index
  "data/spatial_data/covariates/hmi_change.tif", #Human modification index
#  "data/spatial_data/covariates/modis_unsuitable_landcover_500m_2001_2024.tif",
  
  ### PA fract
  "data/spatial_data/covariates/pa_raster_binary.tif", # fraction protected
  
  #### categorical ####
  "data/spatial_data/covariates/pa_raster.tif",# PA 
  "data/spatial_data/covariates/higginsFunctionalBiomes.tif", # functional biome 
  "data/spatial_data/covariates/wwf_olson_biome.tif", # olson biome 
  "data/spatial_data/covariates/KoppenGeigerClimaticRegions1km.tif", # climatic region 
  "data/spatial_data/covariates/GlobalLandCovercopernicus2019.tif"

)


funcs <- c(rep("mean", 11), rep("mode", 5))

covs <- data.table(
  col_name = col_names, 
  cov_path = cov_paths, 
  func = funcs
)

for(i in 1:nrow(covs)) {
  cov_r <- rast(covs[i, ]$cov_path)
  plot(cov_r, main = paste0(covs[i, ]$col_name))
  Sys.sleep(2)
  
}


################################## LOOOOOOOOOOOOP ############################            
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
plan(multisession, workers = 16)
tic()



dt_covs_list <- future_map(1:nrow(covs),
                           .progress = TRUE,
                           .options = furrr_options(seed = TRUE),
                           function(i) {
                             
                             #for(i in 1:nrow(covs)){
                             
                             func <- covs[i, ]$func
                             
                             
                             cov_r <- rast(covs[i, ]$cov_path)
                           #  cov_r <- terra::project(cov_r, st_crs(vect)$wkt)

                             vect_trans <- st_transform(vect, crs = st_crs(cov_r))
                            # vect_trans <- vect
                                                         
                             if(param %in% c("pa_controls", "world_grid")){
                             vect_trans <- vect_trans %>% st_make_valid()
                             vect_trans <- vect_trans[!st_is_empty(vect_trans) & !is.na(st_dimension(vect_trans)), ]
                             }
                           
                             extr <- exactextractr::exact_extract(cov_r, 
                                                                  append_cols = c("unique_id"),
                                                                  vect_trans, 
                                                                  fun = func)
                             
                             setnames(extr, func, covs[i, ]$col_name)
                             
                             dt_extr_fin <- extr %>% 
                               as.data.table() %>%
                               mutate(geom = NULL, x = NULL, geometry = NULL) %>% 
                               unique()
                             
                             print(paste0(covs[i, ]$col_name,"; i = ", i))
                             
                             return(dt_extr_fin)
                             
                           }
)


toc()
plan(sequential)
Sys.time()
dt_covs <- dt_covs_list %>%
  reduce(~ left_join(.x, .y, by = "unique_id"))

#combine
vect_covs <- vect %>% 
  as.data.table() %>% 
  mutate(x = NULL, geom = NULL, geometry = NULL) %>% 
  left_join(dt_covs) %>% 
  as.data.table() %>% 
  mutate(x = NULL, 
         geom = NULL, 
         geometry = NULL) %>% 
  left_join(clim_reg_leg) %>% 
  left_join(biome_leg) %>%
  left_join(fun_biome_leg) %>% 
  left_join(lc_leg) %>% 
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
  dplyr::select(-climatic_region_num, -olson_biome_num, -functional_biome_num, -lc_num)

summary(vect_covs)
glimpse(vect_covs)

if(param == "world_grid"){
  
  fwrite(vect_covs, "data/processed_data/grid_with_covariates.csv")
  
} else if(param == "pas"){
  
  fwrite(vect_covs, "data/processed_data/pas_with_covariates.csv")
  
} else if(param == "pa_controls"){
  
  fwrite(vect_covs, "data/processed_data/pa_controls_with_covariates.csv")
  
} 
gc()