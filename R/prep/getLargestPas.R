## join WDPA 

rm(list = ls())

source("../../../../resources/R_functions/get.heterogeneity.R")


library(sf)
library(tidyverse)
library(data.table)
library(terra)

shp1.raw <- read_sf("../../../../resources/spatial/WDPA/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_0/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

table(shp1.raw$IUCN_CAT)
shp1 <- shp1.raw %>% 
  filter(IUCN_CAT %in% c("Ia", "Ib", "II") | grepl("ional park", DESIG_ENG) | grepl("ional Park", DESIG_ENG)) %>% 
  filter(!grepl("arine", DESIG_ENG)) %>% 
  filter(!grepl("Hunting Area", DESIG_ENG))
  

shp2.raw <- read_sf("../../../../resources/spatial/WDPA/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_1/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

shp2 <- shp2.raw %>% 
  filter(IUCN_CAT %in% c("Ia", "Ib", "II") | grepl("ional park", DESIG_ENG) | grepl("ional Park", DESIG_ENG)) %>% 
  filter(!grepl("arine", DESIG_ENG)) %>% 
  filter(!grepl("Hunting Area", DESIG_ENG))


shp3.raw <- read_sf("../../../../resources/spatial/WDPA/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_2/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

shp3 <- shp3.raw %>% 
  filter(IUCN_CAT %in% c("Ia", "Ib", "II") | grepl("ional park", DESIG_ENG) | grepl("ional Park", DESIG_ENG)) %>% 
  filter(!grepl("arine", DESIG_ENG)) %>% 
  filter(!grepl("Hunting Area", DESIG_ENG))

pas <- rbind(shp1, shp2, shp3)
pas <- pas[!duplicated(data.frame(pas)),]

names(pas)


world <- rnaturalearth::ne_countries() %>% dplyr::select(continent, sovereignt, geometry)

sf_use_s2(FALSE)
pas.cont <- st_join(pas, world)

## select largest X PA's for each country
table(pas.cont$MARINE)
pas.largest <- pas.cont %>% 
  filter(MARINE == 0) %>% 
  mutate(continent = ifelse(sovereignt == "Russia", "Russia", continent)) %>%
  arrange(desc(GIS_AREA)) %>% group_by(sovereignt) %>% slice(1:100000) ### adjust size as wanted 
mapview::mapview(pas.largest)


## check range maps 

library(sf)
library(tidyverse)
library(mapview)
library(terra)
## load data 

## terrestrial mammals 
mam.to.raw <- read_sf("../../../../resources/spatial/IUCN_range_maps/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

mam.to <- mam.to.raw %>% 
  mutate(nativeness = "IUCN_Native") %>% 
  dplyr::select(sci_name, category, nativeness, geometry)
mapview(mam.to[mam.to$sci_name == "Lynx lynx",])

## erick's ranges ------------

mam.ejl <- read_sf("../../../../resources/spatial/Lundgren_introduced_ranges/Native_Intro_Mammals_Updated_Taxonomy.gpkg") %>%
  filter(nativeness %in% c("EJL_Introduced", "IUCN_Introduced"))

unique(mam.ejl$data_quality)
unique(mam.ejl$poly_quality)
names(mam.ejl)

mam.intro <- mam.ejl %>% 
  filter(!data_quality %in% c("Medium","Poor", "Location qlty poor", "Low")) %>% 
  filter(!poly_quality %in% c("Range uncertain, political borders only")) %>% 
  filter(!grepl("Ovis_ammon", orig_Species)) %>%
  mutate(sci_name = gsub("_", " ", orig_Species)) %>% 
  rename(geometry = geom, 
         category = threat_status) %>% 
  dplyr::select(sci_name, category, nativeness, geometry)

mam.ejl[grepl("Ovis", mam.ejl$orig_Species), ]$orig_Species



mam.ejl[grepl("Camel", mam.ejl$orig_Species),]

mapview(mam.intro[grepl("Camel", mam.intro$sci_name),])


## phylacine 

mam.traits <- fread("../../../../resources/traits/Using_PHYLACINE/data/PHYLACINE_1.2.1/Data/Traits/Trait_data.csv") %>% 
  mutate(sci_name = gsub("_", " ", Binomial.1.2), 
         mass_kg = round(Mass.g/1000, 3)) %>% 
  filter(Terrestrial == 1 & !IUCN.Status.1.2 == "EP")

sf_use_s2(FALSE)

mams <- mam.to %>% 
  rbind(mam.intro) %>%
  mutate(range_name = ifelse(nativeness == "IUCN_Native", paste0(sci_name, " native range"), paste0(sci_name, " introduced range"))) %>% 
  left_join(mam.traits) %>% 
  filter(mass_kg >= 10) %>% 
  st_sf() %>% 
  group_by(range_name) %>%
  summarize()

mapview(mams[grepl("Loxodonta", mams$range_name),])
mapview(mams[grepl("Ovis ammon", mams$range_name),])

sf_use_s2(TRUE)

pas.loop <- pas.largest %>% 
  dplyr::select(NAME, WDPA_PID, DESIG_ENG, IUCN_CAT, GIS_AREA, STATUS_YR, geometry, sovereignt, continent) %>% 
  mutate(unique_id = paste0(sovereignt, "_", WDPA_PID))
pas.loop <- pas.loop[!duplicated(data.frame(pas.loop)),]



for(i in 1:nrow(mams)){
  
  range <- mams[i,]$range_name
  
  range.shp <- vect(mams[mams$range_name == range,])
  
  r0.1deg <- terra::rast(res=0.1)
  range.r <- terra::rasterize(range.shp, r0.1deg)
  values(range.r) <- ifelse(is.na(values(range.r)), 0, 1)
  

  extr <- terra::extract(range.r, #
                         vect(pas.loop), #
                         max, na.rm = T) #
  
  setnames(extr, "layer", (range))
  extr$ID <- NULL
  
  pas.loop <- cbind(pas.loop, extr)
  
  print(paste0(range, " done (", i, "/", nrow(mams), ")"  ))  
  
} 

glimpse(pas.loop)

pas.loop <- pas.loop[!duplicated(data.frame(pas.loop)),]

st_write(pas.loop, "data/clean_data/all_pas_with_megafauna_cols.gpkg", append = FALSE)

### extract MAP, MAT and elevation 

## MAT---------------

mat <- rast("../../../../resources/spatial/Chelsa_Climate/CHELSA_bio1_1981-2010_V.2.1.tif") 

mat.extr <- get.heterogeneity(vector = pas.loop, grid = NULL, id.col = "unique_id", raster = mat)
setnames(mat.extr, c("mean"), 
         c("MAT"))

pas.trans.mat <- st_transform(pas.loop, crs = terra::crs(mat))

mat.sd.extr <- terra::extract(mat,
                              pas.trans.mat,
                              fun = sd, na.rm = T)# 

setnames(mat.sd.extr, c("CHELSA_bio1_1981-2010_V.2.1"), 
         c("MAT_sd"))
mat.sd.extr$ID <- NULL

mat.extr.fin <- mat.sd.extr %>% 
  as.data.table() %>% 
  cbind(pas.loop[, "unique_id"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(mat.extr) 


## MAP---------------
map <- rast("../../../../resources/spatial/Chelsa_Climate/CHELSA_bio12_1981-2010_V.2.1.tif") 

map.extr <- get.heterogeneity(vector = pas.loop, grid = NULL, id.col = "unique_id", raster = map)
setnames(map.extr, c("mean"), 
         c("MAP"))

pas.trans.map <- st_transform(pas.loop, crs = terra::crs(map))

map.sd.extr <- terra::extract(map,
                              pas.trans.map,
                              fun = sd, na.rm = T)# 

setnames(map.sd.extr, c("CHELSA_bio12_1981-2010_V.2.1"), 
         c("MAP_sd"))
map.sd.extr$ID <- NULL

map.extr.fin <- map.sd.extr %>% 
  as.data.table() %>% 
  cbind(pas.loop[, "unique_id"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(map.extr) 


## elevation ---------------------------------
ele <- rast("../../../../resources/spatial/Elevation/Elevation_Global_930m.tif")
plot(ele)
ele.extr <- get.heterogeneity(vector = pas.loop, grid = NULL, id.col = "unique_id", raster = ele)
setnames(ele.extr, c("mean"), 
         c("elevation"))

pas.trans.ele <- st_transform(pas.loop, crs = terra::crs(ele))

ele.sd.extr <- terra::extract(ele,
                             pas.trans.ele,
                             fun = sd, na.rm = T)# 

setnames(ele.sd.extr, c("elevation"), 
         c("elevation_sd"))
ele.sd.extr$ID <- NULL

ele.extr.fin <- ele.sd.extr %>% 
  as.data.table() %>% 
  cbind(pas.loop[, "unique_id"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(ele.extr) 


## Biome ---------------------------

wwf.biome <- read_sf("../../../../resources/spatial/WWF_Olson_Ecoregions/WWF_BIOMES.gpkg")
mapview(wwf.biome, zcol = "BIOME_Name")

biome.leg <- wwf.biome %>% as.data.table() %>% mutate(geom = NULL) %>% unique()

rast.tmp <- terra::rast(res = 0.045) #roughly 5 km at equator 
biome.r <- terra::rasterize(wwf.biome, rast.tmp, field = "BIOME")
plot(biome.r)

pas.trans.biome <- pas.loop %>%
  st_transform(crs(biome.r))

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


biome.extr <- terra::extract(biome.r,
                                 pas.trans.biome,
                                 fun = Mode, na.rm = T)# 


biome.extr.fin <- biome.extr %>% 
  as.data.table() %>% 
  cbind(pas.loop[, "unique_id"]) %>% 
  mutate(geometry = NULL) %>% 
  left_join(biome.leg, by = "BIOME") %>% 
  rename(Biome = BIOME_Name) %>% 
  dplyr::select(unique_id, Biome) 


### temperature novelty ---------------------------------------------

temp.nov <- rast("../../../../resources/spatial/Kerr et al - Global Novelty Layers/CLIMATE_MEANTEMP_300.tif")
plot(temp.nov)

#to make the scale more easy to interpret. The larger the number, the longer ago the last equivalent, the more novel. 
values(temp.nov) <- values(temp.nov) + 21800

temp.nov.extr <- get.heterogeneity(vector = pas.loop, grid = NULL, id.col = "unique_id", raster = temp.nov)
setnames(temp.nov.extr, c("mean"), 
         c("temperature_novelty"))
temp.nov.extr$x <- NULL


### precipitation novelty ---------------------------------------------

prec.nov <- rast("../../../../resources/spatial/Kerr et al - Global Novelty Layers/CLIMATE_MEANPREC_300.tif")
plot(prec.nov)
values(prec.nov) <- values(prec.nov) + 21800

prec.nov.extr <- get.heterogeneity(vector = pas.loop, grid = NULL, id.col = "unique_id", raster = prec.nov)
setnames(prec.nov.extr, c("mean"), 
         c("precipitation_novelty"))

prec.nov.extr$x <- NULL



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

bio.nov.raw.extr<- get.heterogeneity(vector = pas.loop, grid = NULL, id.col = "unique_id", raster = bio.nov.raw)
setnames(bio.nov.raw.extr, c("mean"), 
         c("biotic_novelty_raw"))
bio.nov.raw.extr$x <- NULL

bio.nov.extr <- get.heterogeneity(vector = pas.loop, grid = NULL, id.col = "unique_id", raster = bio.nov)
setnames(bio.nov.extr, c("mean"), 
         c("biotic_novelty"))
bio.nov.extr$x <- NULL

#### get the species stuff sorted -------------------------


species.cols <- names(pas.loop %>% dplyr::select(-c(unique_id, NAME, WDPA_PID, DESIG_ENG, IUCN_CAT, STATUS_YR, GIS_AREA, geometry, sovereignt, continent)))
species.cols <- species.cols[!grepl("geometry", species.cols)]
dt.pas <- pas.loop %>% 
  left_join(map.extr.fin) %>% 
  left_join(mat.extr.fin) %>% 
  left_join(ele.extr.fin) %>% 
  left_join(temp.nov.extr) %>% 
  left_join(prec.nov.extr) %>% 
  left_join(bio.nov.raw.extr) %>% 
  left_join(bio.nov.extr) %>% 
  as.data.table() %>% 
  mutate(geometry = NULL) %>% 
  pivot_longer(cols = all_of(species.cols), 
               names_to = "range", 
               values_to = "presence") %>% 
  mutate(species = gsub(".native.range", "", range), 
         species = gsub(".introduced.range", "", species), 
         species = gsub("\\.", "_", species), 
         nativeness = ifelse(grepl("native", range), "native", "introduced"), 
         country = case_when(
           .default = sovereignt,
           NAME == "South Baranof" ~ "United States of America", 
           NAME == "Galápagos" ~ "Ecuador", 
           NAME == "Siberut" ~ "Indonesia", 
           NAME == "Rakiura National Park" ~ "New Zealand", 
         )
         ) %>% 
  filter(!presence == 0) %>% as.data.table() %>% unique()
dt.pas
fwrite(dt.pas, "data/clean_data/5_largest_PAs_per_country_megafauna_cols_long.csv")
unique(dt.pas[grepl("Lille Vildmose", NAME), ]$species)

unique(dt.pas[grepl("Ilulissat Isfjord", NAME), ]$species)
unique(dt.pas[country == "Denmark", ]$NAME)
unique(dt.pas$country)

## Rangifer tarandus range is kinda bullshit

### calculate functional groups etc... 


# Load megafauna traits and groups----------

mam.traits <- fread("../../Woody cover/R_stuff/data/traits/mammal_traits.csv")

mam.traits <- mam.traits %>% 
  dplyr::select(-species) %>% 
  rename(species = Binomial.1.2) %>% 
  mutate(
    mass_kg = ifelse(species == "Bos_primigenius", 400, mass_kg), 
    grazer = ifelse(species == "Bos_primigenius", 1, grazer), 
    mixed_feeder = ifelse(species == "Bos_primigenius", 1, mixed_feeder)) %>% filter(mass_kg >= 10)

mam.traits[mam.traits == "",] <- NA

## combine ----------

#mam.traits[grepl("Sus", Species_), .(Species_, mass_kg)]

pas.traits <- left_join(dt.pas, mam.traits, by = "species")
unique(pas.traits[is.na(pas.traits$mass_kg), species]) #good

## get area 

sf_use_s2(FALSE)
pas.loop$pa_size_km2 <- st_area(pas.loop) %>% units::set_units("km2")
pa_size_km2 <- pas.loop %>% 
  as.data.table() %>%
  mutate(geometry = NULL) %>% 
  units::drop_units() %>% 
  dplyr::select(c(WDPA_PID, pa_size_km2)) %>% as.data.table()



pas.traits2 <- pas.traits %>% 
  mutate(herbi_or_carni = ifelse(!is.na(trophic_level_carni), "Predator", "Herbivore"),
         one_per_row = 1,
         maximum_prey_size = case_when(
           .default = maximum_prey_size, 
           herbi_or_carni == "Predator" & mass_kg < 15 ~ mass_kg,
           herbi_or_carni != "Predator" ~ 0, 
           species == "Orycteropus_afer" ~ 0.1, 
           species == "Helarctos_malayanus" ~ 0.5, 
           species == "Smutsia_gigantea" ~ 0.1, 
           species == "Myrmecophaga_tridactyla" ~ 0.1, 
           species == "Priodontes_maximus" ~ 0.1,
           species == "Melursus_ursinus" ~ 0.1),
         largest_accessible_prey_size = 
           case_when(
             .default = largest_accessible_prey_size, 
             herbi_or_carni == "Predator" & mass_kg < 15 ~ mass_kg,
             herbi_or_carni != "Predator" ~ 0, 
             species == "Orycteropus_afer" ~ 0.1, 
             species == "Helarctos_malayanus" ~ 0.5, 
             species == "Smutsia_gigantea" ~ 0.1, 
             species == "Myrmecophaga_tridactyla" ~ 0.1, 
             species == "Priodontes_maximus" ~ 0.1,
             species == "Melursus_ursinus" ~ 0.1)
  ) %>% unique()

unique(pas.traits2[is.na(maximum_prey_size), ]$species)
pas.traits2[is.na(Herbi)]$WDPA_PID

pas.traits2[, n_species_pa := uniqueN(species), by = unique_id]

pas.traits2[, maximum_prey_size_pa := max(maximum_prey_size, na.rm = TRUE), by = unique_id]
quantile(pas.traits2$maximum_prey_size_pa)
hist(pas.traits2$n_species_pa)


quantile(pas.traits2$largest_accessible_prey_size, na.rm = T)
pas.traits2[, largest_accessible_prey_size_pa := max(largest_accessible_prey_size), by = unique_id]
quantile(pas.traits2$largest_accessible_prey_size_pa)


#### elephants -------

pas.ele <- pas.traits2[species %in% c("Loxodonta_africana", "Loxodonta_cyclotis", "Elephas_maximus"), ]
pas.ele[, elephant_yn := "yes", by = unique_id]
pas.ele2 <- pas.ele %>% dplyr::select(unique_id, elephant_yn) %>% unique()


#### carnivores -----
pas.pred <- pas.traits2[!herbi_or_carni == "Herbivore", ]
pas.pred[, n_carni_sp_pa := uniqueN(species, na.rm = TRUE), by = unique_id]
pas.pred2 <- pas.pred %>% dplyr::select(unique_id, n_carni_sp_pa) %>% unique()


quantile(pas.pred$n_carni_sp_pa)

#### herbivores ----
pas.herbi_or_carni <- pas.traits2[herbi_or_carni == "Herbivore", ]
pas.herbi[, n_herbi_sp_pa := uniqueN(species, na.rm = TRUE), by = unique_id]
quantile(pas.herbi$n_herbi_sp_pa)



pas.herbi2 <- pas.herbi_or_carni %>% mutate(
  regular_predator_present = ifelse(mass_kg < largest_accessible_prey_size_pa, 1, 0), 
  occasional_predator_present =  ifelse(mass_kg < maximum_prey_size_pa, 1, 0), 
) %>% filter(!n_herbi_sp_pa < 1) %>% unique()


pas.herbi2[, species_w_occasional_predator := sum(occasional_predator_present), by = unique_id]
pas.herbi2[, species_w_regular_predator := sum(regular_predator_present), by = unique_id]
quantile(pas.herbi2$species_w_regular_predator)
quantile(pas.herbi2$species_w_occasional_predator)

quantile(pas.herbi2$n_herbi_sp_pa)

unique(pas.herbi2$species_w_regular_predator)
unique(pas.herbi$WDPA_PID)
n_distinct(pas.traits2$NAME)
n_distinct(pas.traits2$WDPA_PID)
n_distinct(pas.traits2$unique_id)

pas.herbi3 <- pas.herbi2 %>% mutate(
  prop_species_w_regular_predator = species_w_regular_predator/n_herbi_sp_pa, 
  prop_species_w_occasional_predator = species_w_occasional_predator/n_herbi_sp_pa) %>%
  dplyr::select(
    c(prop_species_w_regular_predator, prop_species_w_occasional_predator, n_herbi_sp_pa,
      unique_id)
  ) %>% as.data.table() %>% 
  units::drop_units() %>% 
  mutate(prop_species_w_regular_predator = ifelse(prop_species_w_regular_predator > 1, 1, prop_species_w_regular_predator),
         prop_species_w_occasional_predator = ifelse(prop_species_w_occasional_predator > 1, 1, prop_species_w_occasional_predator)) %>% 
  unique()

summary(pas.herbi3)

## Browser ---------- 
pas.browsers <- pas.traits2[guild_only_herbivory == "Browser"]
pas.browsers[, n_browser_sp_pa := uniqueN(species, na.rm = TRUE), by = unique_id]
pas.browsers2 <- pas.browsers %>% 
  filter(!is.na(n_browser_sp_pa)) %>% 
  dplyr::select(n_browser_sp_pa, unique_id) %>% 
  unique()

## Grazer ----------------
pas.grazers <- pas.traits2[guild_only_herbivory == "Grazer"]
pas.grazers[, n_grazer_sp_pa := uniqueN(species, na.rm = TRUE), by = unique_id]
pas.grazers2 <- pas.grazers %>% 
  filter(!is.na(n_grazer_sp_pa)) %>% 
  dplyr::select(n_grazer_sp_pa, unique_id) %>% 
  unique()

## Mixed feeder ----------------
pas.mf <- pas.traits2[guild_only_herbivory == "Mixed Feeder"]
pas.mf[, n_mixed_feeder_sp_pa := uniqueN(species, na.rm = TRUE), by = unique_id]
pas.mf2 <- pas.mf %>% 
  filter(!is.na(n_mixed_feeder_sp_pa)) %>% 
  dplyr::select(n_mixed_feeder_sp_pa, unique_id) %>% 
  unique()

## Natives -------------------
pas.nat <- pas.traits2[nativeness == "native"]
pas.nat[, n_native_sp_pa := uniqueN(species, na.rm = TRUE), by = unique_id]
pas.nat2 <- pas.nat %>% 
  filter(!is.na(n_native_sp_pa)) %>% 
  dplyr::select(n_native_sp_pa, unique_id) %>% 
  unique()

## Intros --------
pas.intro <- pas.traits2[!nativeness == "native"]
pas.intro[, n_introduced_sp_pa := uniqueN(species, na.rm = TRUE), by = unique_id]
pas.intro2 <- pas.intro %>% 
  filter(!is.na(n_introduced_sp_pa)) %>% 
  dplyr::select(n_introduced_sp_pa, unique_id) %>% 
  unique()



pas.traits3 <- pas.traits2 %>% 
  left_join(pas.herbi3) %>% 
  left_join(pas.pred2) %>% 
  left_join(pas.ele2) %>% 
  left_join(pas.browsers2) %>% 
  left_join(pas.grazers2) %>% 
  left_join(pas.mf2) %>% 
  left_join(pas.nat2) %>% 
  left_join(pas.intro2) %>% 
  mutate(
    n_grazer_sp_pa = ifelse(!is.na(n_grazer_sp_pa), n_grazer_sp_pa, 0), 
    n_browser_sp_pa = ifelse(!is.na(n_browser_sp_pa), n_browser_sp_pa, 0), 
    n_mixed_feeder_sp_pa = ifelse(!is.na(n_mixed_feeder_sp_pa), n_mixed_feeder_sp_pa, 0), 
    n_introduced_sp_pa = ifelse(!is.na(n_introduced_sp_pa), n_introduced_sp_pa, 0), 
    elephant_yn = ifelse(!is.na(elephant_yn), elephant_yn, "no")
  ) %>%
  mutate(grazer_browser_ratio = ifelse(n_grazer_sp_pa == 0 | n_browser_sp_pa == 0, NA, n_grazer_sp_pa/n_browser_sp_pa),
         grazer_mixed_ratio = ifelse(n_mixed_feeder_sp_pa == 0 | n_browser_sp_pa == 0, NA, n_grazer_sp_pa/n_mixed_feeder_sp_pa),
         browser_mixed_ratio = ifelse(n_mixed_feeder_sp_pa == 0 & n_browser_sp_pa == 0, NA, n_browser_sp_pa/n_mixed_feeder_sp_pa), 
         percent_native_sp = n_native_sp_pa/n_species_pa, 
         percent_introduced_sp = n_introduced_sp_pa/n_species_pa, 
         intro_native_ratio = n_introduced_sp_pa/n_native_sp_pa
         )

  
n_distinct(pas.traits3$herbivore_fg_own)

pas.traits3[, `:=` (max_species_body_mass = max(mass_kg, na.rm = TRUE),
              mean_species_body_mass = mean(mass_kg, na.rm = TRUE),
              median_species_body_mass = median(mass_kg, na.rm = TRUE),
              herbivore_fg = n_distinct(herbivore_fg_own, na.rm = TRUE),
              predator_fg = n_distinct(carnivore_fg_own, na.rm = TRUE),
              trophic_groups = n_distinct(trophic_level, na.rm = TRUE),
              predator_fg_simple = n_distinct(trophic_level_carni, na.rm = TRUE),
              herbivore_fg_simple = n_distinct(trophic_level_herbivore, na.rm = TRUE)), by = unique_id]
unique(pas.traits3$herbivore_fg)


names(pas.traits3)
dt.final <- unique(pas.traits3[,.(unique_id, WDPA_PID, NAME,  STATUS_YR, country, continent, 
                                 MAT, MAT_sd, MAP, MAP_sd, elevation, elevation_sd, 
                                 max_species_body_mass, mean_species_body_mass, median_species_body_mass, 
                           herbivore_fg, predator_fg, trophic_groups, predator_fg_simple,  herbivore_fg_simple, 
                           n_species_pa, n_herbi_sp_pa, n_carni_sp_pa, 
                           n_mixed_feeder_sp_pa, n_grazer_sp_pa, n_browser_sp_pa,
                           grazer_browser_ratio, grazer_mixed_ratio, browser_mixed_ratio,
                           n_native_sp_pa, n_introduced_sp_pa, percent_native_sp, percent_introduced_sp, intro_native_ratio,
                           prop_species_w_occasional_predator, prop_species_w_regular_predator,
                           elephant_yn,
                           GIS_AREA)]) #%>% filter(!is.na())

unique(dt.final$percent_introduced_sp)

nrow(dt.final[STATUS_YR <= 2013])
nrow(dt.final[STATUS_YR <= 2003])



fwrite(x = dt.final, file = "data/clean_data/all_pas_w_megafauna_and_novelty.csv")

pa.shapes <- pas.loop %>% dplyr::select(unique_id, WDPA_PID, NAME)
st_write(pas.loop, "data/clean_data/all_pas_w_megafauna_and_novelty.gpkg", append = FALSE)

st_write(pas.loop, "data/clean_data/all_pas_w_megafauna_and_novelty.gpkg", append = FALSE)
