
source("scripts/functions/get.heterogeneity.R")
source("scripts/functions/ask.log.R")

library(MuMIn)
library(glmmTMB)
library(terra)
library(tidyverse)
library(data.table)
library(sf)
library(tidylog)
library(GGally)
library(scales)
library(rnaturalearth)
library(gridExtra)
library("qgam")
library(sf)
library(terra)

dt.raw <- fread("data/clean_data/dt_with_covs.csv") 

dt <- dt.raw %>%
  unique() %>% 
  mutate(n_herbi_sp_reserve = ifelse(is.na(n_herbi_sp_reserve), 0, n_herbi_sp_reserve)) %>% 
  mutate(herbi_biomass_ha = ifelse(is.na(herbi_biomass_ha), 0, herbi_biomass_ha), 
         days_since_last_fire = ifelse(days_since_last_fire == 0, 7000, days_since_last_fire), 
         fire_events_since_2001 = ifelse(is.na(fire_events_since_2001), 0, fire_events_since_2001))


## load EVI from 2019 which corresponds to the year of the tree cover layer 
dt.shape <-  st_read("data/clean_data/prelim_dt_with_covs.shp")  %>%
  rename(reserve_name = rsrv_nm) %>% left_join(dt.raw)

### build grids 

sa.sh <- ne_countries(country = "south africa", type = "map_units")

sa.sh.utm <- st_transform(sa.sh, crs = 22235)

grid5k <- st_make_grid(sa.sh.utm, cellsize = c(5000, 5000), what = "polygons", square = F) %>%
  st_as_sf() %>% st_intersection(sa.sh.utm) %>% 
  mutate(grid_ID = paste0("grid_", 1:nrow(.)), 
         area_km2 = as.numeric(st_area(.)/1000000)) %>% 
  dplyr::select(c("grid_ID", "area_km2")) %>% 
  filter(area_km2 > quantile(area_km2, 0.02)) %>% st_transform(crs = 4326) #to keep only complete cells 

quantile(grid5k$area_km2, 0.5)

mapview::mapview(grid5k)

sf_use_s2(FALSE)
grid1k <- st_make_grid(sa.sh.utm, cellsize = c(1000, 1000), what = "polygons", square = T) %>%
  st_as_sf() %>% st_intersection(sa.sh.utm) %>% 
  mutate(grid1k_ID = paste0("grid1k_", 1:nrow(.)), 
         area_km2 = as.numeric(st_area(.)/1000000)) %>% 
  dplyr::select(c("grid1k_ID", "area_km2")) %>% 
  filter(area_km2 >= 1) %>% st_transform(crs = 4326) #to keep only complete cells 


## Extract covariates ------------------------------------------

### Landcover ---------------
lc <- rast("../../../../resources/spatial/Landcover_ZAF/lc_sa_discrete.tif") 

plot(lc)


discrete_classification <- c(0, 20, 30, 40, 50, 60, 70, 80, 90, 
                       100, 111, 112, 113, 114, 115, 116, 121, 
                       122, 123, 124, 125, 126, 200)

landcover <- c("unknown", "shrubs", "herbaceous_vegetation",
               "agriculture", "urban", "bare_or_sparse_vegetation", "snow", 
               "permanent_water", "herbaceous_wetland", 
               "moss", "evergreen_needle_forest", "evergreen_broad_leaf_forest", 
               "deciduous_needle_leaf_forest", "deciduous_broad_leaf_forest", 
               "mixed_forest", "other_forest", "open_evergreen_needle_forest", 
               "open_evergreen_broad_leaf_forest", "open_deciduous_needle_forest", 
               "open_deciduous_brod_leaf_forest", "open_mixed_forest", 
               "open_other_forest", "large_water_bodies")

lc.leg <- data.table(landcover = landcover, discrete_classification = discrete_classification)


Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

#### grid--------
grid5k.trans <- grid5k %>%
  st_transform(crs(lc))

grid.lc.extraction <- terra::extract(lc,
                                    grid5k.trans,
                                    fun = Mode, na.rm = T)# 



setDT(grid.lc.extraction)

grid.lc.extraction <- grid.lc.extraction %>% rename(discrete_classification = discrete_classification_mode) %>% 
  left_join(lc.leg)

grid5k.1 <- grid5k %>% 
  cbind(grid.lc.extraction) %>% 
  filter(!landcover %in% c("unknown", "agriculture", "urban", "large_water_bodies", "permanent_water")) ## remove unwanted landuse categories 
mapview::mapview(grid5k.1, zcol = "landcover")

table(grid5k.1$landcover)
#### reserves ---------

dt.shape.trans <- dt.shape %>%
  st_transform(crs(lc))

pa.lc.extraction <- terra::extract(lc,
                                   dt.shape.trans,
                                   fun = Mode, na.rm = T)# 



setDT(pa.lc.extraction)

pa.lc.extraction <- pa.lc.extraction %>% rename(discrete_classification = discrete_classification_mode) %>% 
  left_join(lc.leg)

dt.shape.1 <- dt.shape %>% 
  cbind(pa.lc.extraction) %>% 
  filter(!landcover %in% c("unknown", "agriculture", "urban", "large_water_bodies", "permanent_water"))

table(pa.lc.extraction$landcover)

### Mean IAP cover -------------------------
iap.mean <- rast("../../../../resources/spatial/WfW_SA_invasive_plant_cover/iap_cover_mean_100m.tif")
plot(iap.mean)

#### grids ------------------------

iap.mean.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = iap.mean)
setnames(iap.mean.extr.grid, c("mean"), 
         c("iap_cover_mean"))

iap.mean.extr.grid$x <- NULL

grid5k.2 <- grid5k.1 %>% 
  left_join(iap.mean.extr.grid, by = "grid_ID")  %>%
  unique()

#### reserves -------------------------

iap.mean.extr <- get.heterogeneity(vector = dt.shape.1, grid = NULL, id.col = "reserve_name", raster = iap.mean)
setnames(iap.mean.extr, c("mean"), 
         c("iap_cover_mean"))

dt.shape.2 <- dt.shape.1 %>% 
  as.data.table() %>% 
  mutate(geometry = NULL) %>%
  left_join(iap.mean.extr, by = "reserve_name")  %>%
  unique()

### Max IAP cover -------------------------
iap.max <- rast("../../../../resources/spatial/WfW_SA_invasive_plant_cover/iap_cover_max_100m.tif")
plot(iap.max)


#### grids -------------------------------

grid5k.trans <- grid5k.1 %>%
  st_transform(crs(iap.max))

iap.max.extr.grid <- terra::extract(iap.max, # the rast layers
                               vect(grid5k.trans), # the spatial polygons, which you have to convert to a terra format
                               max, na.rm = T) # since we're dealing with polygons we need to summarize per overlapping pixel somehow

setDT(iap.max.extr.grid)
names(iap.max.extr.grid) <- c("ID", "iap_cover_max")
iap.max.extr.grid$ID <- NULL

grid5k.3 <- grid5k.2 %>% 
  cbind(iap.max.extr.grid)  %>%
  unique()


#### reserves -------------------------
v.trans <- dt.shape.1 %>%
  st_transform(crs(iap.max))

iap.max.extr <- terra::extract(iap.max, # the rast layers
                                  vect(v.trans), # the spatial polygons, which you have to convert to a terra format
                                  max, na.rm = T) # since we're dealing with polygons we need to summarize per overlapping pixel somehow

setDT(iap.max.extr)
names(iap.max.extr) <- c("ID", "iap_cover_max")
iap.max.extr$ID <- NULL

dt.shape.3 <- dt.shape.2 %>% 
  cbind(iap.max.extr)  %>%
  unique()


### temperature novelty ---------------------------------------------

temp.nov <- rast("../../../../resources/spatial/Kerr et al - Global Novelty Layers/CLIMATE_MEANTEMP_300.tif")
plot(temp.nov)

#to make the scale more easy to interpret. The larger the number, the longer ago the last equivalent, the more novel. 
values(temp.nov) <- values(temp.nov) + 21800

#### grids ------------------------

temp.nov.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = temp.nov)
setnames(temp.nov.extr.grid, c("mean"), 
         c("temperature_novelty"))

temp.nov.extr.grid$x <- NULL

grid5k.4 <- grid5k.3 %>% 
  left_join(temp.nov.extr.grid, by = "grid_ID")  %>%
  unique()

#### reserves -------------------------

temp.nov.extr <- get.heterogeneity(vector = dt.shape.1, grid = NULL, id.col = "reserve_name", raster = temp.nov)
setnames(temp.nov.extr, c("mean"), 
         c("temperature_novelty"))


dt.shape.4 <- dt.shape.3 %>% 
  left_join(temp.nov.extr, by = "reserve_name")  %>%
  unique()

### precipitation novelty ---------------------------------------------

prec.nov <- rast("../../../../resources/spatial/Kerr et al - Global Novelty Layers/CLIMATE_MEANPREC_300.tif")
plot(prec.nov)
values(prec.nov) <- values(prec.nov) + 21800


#### grids ------------------------

prec.nov.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = prec.nov)
setnames(prec.nov.extr.grid, c("mean"), 
         c("precipitation_novelty"))

prec.nov.extr.grid$x <- NULL

grid5k.5 <- grid5k.4 %>% 
  left_join(prec.nov.extr.grid, by = "grid_ID")  %>%
  unique()

#### reserves ------------------

prec.nov.extr <- get.heterogeneity(vector = dt.shape.1, grid = NULL, id.col = "reserve_name", raster = prec.nov)
setnames(prec.nov.extr, c("mean"), 
         c("precipitation_novelty"))


dt.shape.5 <- dt.shape.4 %>% 
  left_join(prec.nov.extr, by = "reserve_name")  %>%
  unique()


### biotic novelty ---------------------------------------------

bio.nov.raw <- rast("../../../../resources/spatial/Kerr et al - Global Novelty Layers/GLOBALISATION 1950/GLOBALISATION_BIOTIC.tif")
plot(bio.nov.raw)

gbif.sample <- rast("../../../../resources/spatial/Kerr et al - Global Novelty Layers/GLOBALISATION 1950/GBIF_SAMPLE.tif")
plot(gbif.sample)
quantile(values(gbif.sample, na.rm = T))


gbif.richness <- rast("../../../../resources/spatial/Kerr et al - Global Novelty Layers/GLOBALISATION 1950/GBIF_RICHNESS.tif")
plot(gbif.richness)
quantile(values(gbif.richness, na.rm = T))


#clean gbif mess
bio.nov <- bio.nov.raw

#set cells with less than 10 samples to NA 
values(bio.nov) <- ifelse(values(gbif.sample) < 10, NA, values(bio.nov))

#set cells with less than 5 species to NA
values(bio.nov) <- ifelse(values(gbif.richness) < 5, NA, values(bio.nov))

#### grids ------------------------

bio.nov.raw.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = bio.nov.raw)
setnames(bio.nov.raw.extr.grid, c("mean"), 
         c("biotic_novelty_raw"))
bio.nov.raw.extr.grid$x <- NULL

bio.nov.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = bio.nov)
setnames(bio.nov.extr.grid, c("mean"), 
         c("biotic_novelty"))

bio.nov.extr.grid$x <- NULL

grid5k.6 <- grid5k.5 %>% 
  left_join(bio.nov.raw.extr.grid, by = "grid_ID")  %>%
  unique() %>% 
  left_join(bio.nov.extr.grid, by = "grid_ID")  %>%
  unique()



#### reserves ------------------------

bio.nov.raw.extr <- get.heterogeneity(vector = dt.shape.1, grid = NULL, id.col = "reserve_name", raster = bio.nov.raw)
setnames(bio.nov.raw.extr, c("mean"), 
         c("biotic_novelty_raw"))
bio.nov.raw.extr$x <- NULL

bio.nov.extr <- get.heterogeneity(vector = dt.shape.1, grid = NULL, id.col = "reserve_name", raster = bio.nov)
setnames(bio.nov.extr, c("mean"), 
         c("biotic_novelty"))

bio.nov.extr$x <- NULL

dt.shape.6 <- dt.shape.5 %>% 
  left_join(bio.nov.raw.extr, by = "reserve_name")  %>%
  unique() %>% 
  left_join(bio.nov.extr, by = "reserve_name")  %>%
  unique()

## extract remaining stuff for the grids ------------

### Canopy height --------
ch <- rast("../../../../resources/spatial/Canopy_height_Tolan_etal_SA/sa_canopy_height_100m_tolan_et_al_2024.tif")
res(ch)
plot(ch)


ch.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = grid1k, id.col = "grid_ID", raster = ch)
setnames(ch.extr.grid, c("id.col", "grid.cv", "grid.sd", "grid.mean", "mean"), 
         c("grid_ID", "canopy_height_cv", "canopy_height_sd", "canopy_height_mean_grid", "canopy_height_mean"))
ch.extr.grid$geometry <- NULL
ch.extr.grid$x <- NULL
ch.extr.grid$canopy_height_mean_grid <- NULL


### elevation --------

ele <- rast("../../../../resources/spatial/Elevation_ZAF/Elevation_SA_90m.tif")
crs(ele) # 
res(ele)

ele.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = grid1k, id.col = "grid_ID", raster = ele)

setnames(ele.extr.grid, c("id.col", "grid.cv", "grid.sd", "grid.mean", "mean"), 
         c("grid_ID", "elevation_cv", "elevation_sd", "elevation_mean_grid", "elevation_mean"))
ele.extr.grid$x <- NULL
ele.extr.grid$elevation_mean_grid <- NULL

### MAT---------------

mat <- rast("../../../../resources/spatial/Chelsa_Climate/CHELSA_bio1_1981-2010_V.2.1.tif") 

mat.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = mat)
setnames(mat.extr.grid, c("mean"), 
         c("MAT"))
mat.extr.grid$x <- NULL


### MAP---------------
map <- rast("../../../../resources/spatial/Chelsa_Climate/CHELSA_bio12_1981-2010_V.2.1.tif") 
plot(map)

map.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = map)
setnames(map.extr.grid, c("mean"), 
         c("MAP"))
map.extr.grid$x <- NULL

### Tree cover ---------------
tc <- rast("../../../../resources/spatial/Africa_Tree_Cover_Reiner_et_al_2023/ps_africa_treecover_2019_100m_v1.1.tif") 

tc.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = grid1k, id.col = "grid_ID", raster = tc)
setnames(tc.extr.grid, c("id.col", "grid.cv", "grid.sd", "grid.mean", "mean"), 
         c("grid_ID", "tree_cover_cv", "tree_cover_sd", "tree_cover_mean_grid", "tree_cover_mean"))
tc.extr.grid$x <- NULL
tc.extr.grid$tree_cover_mean_grid <- NULL

### evi mean and heterogeneity -------------------------------

evi <- rast("../../../../resources/spatial/SA_Terra_Vegetation_Indices/mean_EVI_for_All_2017.tif")


evi.grid <- get.heterogeneity(vector = grid5k.1, grid = grid1k, id.col = "grid_ID", raster = evi)
setnames(evi.grid, c("id.col", "grid.cv", "grid.sd", "grid.mean", "mean"), 
         c("grid_ID", "evi_cv", "evi_sd", "evi_mean_grid", "evi_mean"))
evi.grid$x <- NULL
evi.grid$evi_mean_grid <- NULL


### combine
grid5k.7 <- grid5k.6 %>% 
  left_join(ch.extr.grid, by = "grid_ID")  %>%
  unique() %>% 
  left_join(ele.extr.grid, by = "grid_ID")  %>%
  unique() %>% 
  left_join(mat.extr.grid, by = "grid_ID")  %>%
  unique() %>% 
  left_join(map.extr.grid, by = "grid_ID")  %>%
  unique() %>% 
  left_join(tc.extr.grid, by = "grid_ID")  %>%
  unique() %>% 
  left_join(evi.grid, by = "grid_ID")  %>%
  unique()


### Biome -------------------------


veg <- sf::read_sf("../../../../resources/spatial/Vegetation_Map_SA_2018/NVM2018_AEA_V22_7_16082019_final.shp")


veg <- st_transform(veg, crs = 4326)
veg <- st_sf(veg) %>% 
  mutate(BIOMEID_18 = ifelse(veg$BIOME_18 == "Albany Thicket", 5, veg$BIOMEID_18))
# Perform a spatial join to associate "pa" polygons with "veg" polygons

table(veg$BIOMEID_18)

# Get the extent of the "veg" data
veg_bbox <- st_bbox(veg)

# Calculate the desired resolution in degrees
desired_resolution <- c(0.01, 0.01)  # Approximately 1 km at the equator

# Calculate the number of rows and columns based on the desired resolution
num_rows <- as.integer((32.9-16.4) / desired_resolution[2])
num_cols <- as.integer((34.8-22.1) / desired_resolution[1])

# Create a raster with the matching extent and resolution
matching_raster <- rast(
  extent = veg_bbox,
  nrows = num_rows,
  ncols = num_cols,
  crs = crs(veg)
)

# Rasterize the "veg" polygons into the matching raster


table(veg$BIOME_18)
veg_raster <- rasterize(veg, matching_raster, field = "BIOMEID_18")
#plot(veg_raster)


leg <- unique(veg %>% dplyr::select(c("BIOME_18", "BIOMEID_18")))
leg$geometry <- NULL
setDT(leg)
leg <- unique(leg)
leg <- leg[!BIOMEID_18 == 0]
leg

grid5k.1.trans <- grid5k.1 %>%
  st_transform(crs(veg_raster))

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


veg.extraction.grid <- terra::extract(veg_raster,
                                 grid5k.1.trans,
                                 fun = Mode, na.rm = T)# 


table(veg.extraction.grid$BIOMEID_18)

#veg.df <- as.data.frame(veg_raster, xy = T)

veg.final <- veg.extraction.grid %>% 
  # mutate(BIOMEID_18 = ifelse(BIOMEID_18 %in% c(5.5, 7.5, 8.5), NA, BIOMEID_18)) %>% 
  as.data.table() %>% 
  cbind(grid5k.1[, "grid_ID"]) %>% 
  mutate(x = NULL) %>% 
  left_join(leg, by = "BIOMEID_18") %>% 
  rename(Biome = BIOME_18) %>% 
  dplyr::select(grid_ID, Biome) 
table(veg.final$Biome)


grid5k.8 <- grid5k.7 %>% 
  left_join(veg.final, by = "grid_ID")  %>% 
  unique() 


## Fire frequency ---------------
fir <- rast("../../../../resources/spatial/Fire/FireEventsBetween2001And2024_SA.tif")
#plot(fir)
range(values(fir), na.rm = T)


values(fir) <- ifelse(is.na(values(fir)), 0, values(fir))


fir.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = fir)
setnames(fir.extr.grid, c("mean"), 
         c("fire_events_since_2001"))
fir.extr.grid$x <- NULL


## Proportion burned area ---------------
dsf <- rast("../../../../resources/spatial/Fire/DaysSinceLastFireJune2019_SA.tif")
#plot(dsf)

dsf.burned <- dsf
values(dsf.burned) <- ifelse(values(dsf) < 3650 & values(dsf) > 1, 1, 0)

plot(dsf.burned)


dsf.burned.extr.grid <- get.heterogeneity(vector = grid5k.1, grid = NULL, id.col = "grid_ID", raster = dsf.burned)
setnames(dsf.burned.extr.grid, c("mean"), 
         c("prop_burned_area"))
dsf.burned.extr.grid$x <- NULL

grid5k.9 <- grid5k.8 %>% 
  left_join(dsf.burned.extr.grid, by = "grid_ID")  %>% 
  left_join(fir.extr.grid, by = "grid_ID")  %>% 
  unique() 


############ SAVE DATA 

#grids ------------
names(grid5k.9)
dt.grid5k <- grid5k.9 %>% 
  as.data.table() %>% 
  select(-c(discrete_classification, ID)) %>%
  mutate(x = NULL)

fwrite(dt.grid5k, "data/clean_data/sa_novelty_grids.csv")

shp.grid5k <- grid5k.9 %>% 
  select(c(grid_ID, x)) 
st_write(shp.grid5k, "data/clean_data/sa_novelty_grids.gpkg", append = TRUE)



### reserves 
names(dt.shape.6)
dt.reserves <- dt.shape.6 %>% 
  select(-c(discrete_classification, ID)) %>%
  as.data.table() %>% 
  mutate(geometry = NULL)

fwrite(dt.reserves, "data/clean_data/sa_novelty_reserves.csv")

shp.reserves <- dt.shape.1 %>% 
  select(c(geometry, reserve_name)) 

st_write(shp.reserves, "data/clean_data/sa_novelty_reserves.gpkg", append = FALSE)

ggplot() + 
  geom_point(data = dt.reserves, aes(y = iap_cover_mean, x = herbi_biomass_ha))

ggplot() + 
  geom_point(data = dt.reserves, aes(y = iap_cover_mean, x = herbivore_fg))

ggplot() + 
  geom_point(data = dt.reserves, aes(y = iap_cover_mean, x = n_herbi_sp_reserve))

ggplot() + 
  geom_point(data = dt.reserves, aes(y = iap_cover_mean, x = MAP))



####### Corr plot #######
library(ggcorrplot)
cor.grid <- grid5k.9  %>% as.data.table() %>% select(-c(x, area_km2, grid_ID, landcover, discrete_classification, ID, Biome)) %>% filter(complete.cases(.))
corr <- round(cor(cor.grid), 1)
head(corr[, 1:6])

ggcorrplot(corr, hc.order = F, type = "lower",
           lab = TRUE)

# log variables that need to get pushed towards a normal distribution 
dt.mod <- dt.reserves %>% 
  mutate(
    MAT_scaled = as.numeric(scale(MAT)),
    MAP_scaled = as.numeric(scale(MAP)),
    log_days_since_last_fire_scaled = as.numeric(scale(log1p(days_since_last_fire))),
    log_fire_events_since_2001_scaled = as.numeric(scale(log1p(fire_events_since_2001))), 
    prop_burned_area_scaled = as.numeric(scale(prop_burned_area)), 
    
    
    temperature_novelty_trans = 1/(temperature_novelty*-1 + 0.001), 
    precipitation_novelty_trans = 1/(precipitation_novelty*-1 + 0.001), 

    
    biotic_novelty_raw_scaled = scale(biotic_novelty_raw), 
    biotic_novelty_scaled = scale(biotic_novelty), 
    
    precipitation_novelty_scaled = scale(precipitation_novelty_trans), 
    temperature_novelty_scaled = scale(temperature_novelty_trans), 
    iap_cover_mean_scaled = scale(iap_cover_mean), 
    iap_cover_max_scaled = scale(iap_cover_max), 
    
    
    
    log_tree_cover_mean_scaled = as.numeric(scale(log1p(tree_cover_mean))), 
    log_canopy_height_mean_scaled = as.numeric(scale(log1p(canopy_height_mean))), 
    
    herbivore_fg_scaled = as.numeric(scale(herbivore_fg)), 
    log_herbi_biomass_ha_scaled = as.numeric(scale(log1p(herbi_biomass_ha))), 
    log_grazer_biomass_ha_scaled = as.numeric(scale(log1p(grazer_biomass_ha))), 
    log_browser_biomass_ha_scaled = as.numeric(scale(log1p(browser_biomass_ha))), 
    log_mixed_feeder_biomass_ha_scaled = as.numeric(scale(log1p(mixed_feeder_biomass_ha))), 
    n_herbi_sp_reserve_scaled = as.numeric(scale(n_herbi_sp_reserve)), 
    log_CW_mean_species_body_mass_scaled = as.numeric(scale(log1p(CW_mean_species_body_mass))), 
    Biome = as.factor(Biome), 
    elephant_yn = ifelse(elephant_biomass_ha > 0, "elephants", "no_elephants")
  ) %>% filter(!is.na(CW_mean_species_body_mass))



## Dredge GLMM ----------------------------------
global.glmm.tc <- glmmTMB(tree_cover_mean ~
                            MAT_scaled +
                            MAP_scaled +
                            biotic_novelty_raw_scaled + 
                            biotic_novelty_scaled + 
                    
                            iap_cover_max_scaled + 
                            iap_cover_mean_scaled + 
                            temperature_novelty_scaled + 
                            precipitation_novelty_scaled + 
                            #log_days_since_last_fire_scaled +
                            log_fire_events_since_2001_scaled +
                            #prop_burned_area_scaled +
                            herbivore_fg_scaled +
                            log_herbi_biomass_ha_scaled +
                            log_grazer_biomass_ha_scaled +
                            log_browser_biomass_ha_scaled +
                            log_mixed_feeder_biomass_ha_scaled +
                            n_herbi_sp_reserve_scaled +
                            log_CW_mean_species_body_mass_scaled, 
                          data = dt.mod)
summary(global.glmm.tc)
r.squaredGLMM(global.glmm.tc)



options(na.action = "na.fail")
dd.tc <- MuMIn::dredge(global.glmm.tc, beta = "sd", trace = 2)
options(na.action = "na.omit")

#best model
summary(get.models(dd.tc, 1)[[1]])
b.glmm.tc <- get.models(dd.tc, 1)[[1]]
r.squaredGLMM(b.glmm.tc) 
