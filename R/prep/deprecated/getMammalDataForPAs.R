## join WDPA 

rm(list = ls())


library(sf)
library(tidyverse)
library(data.table)
library(terra)

shp1.raw <- read_sf("data/spatialData/protectedAreas/rawPAs/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_0/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

table(shp1.raw$IUCN_CAT)
shp1 <- shp1.raw %>% 
  filter(IUCN_CAT %in% c("Ia", "Ib", "II") | grepl("ional park", DESIG_ENG) | grepl("ional Park", DESIG_ENG)) %>% 
  filter(!grepl("arine", DESIG_ENG)) %>% 
  filter(!grepl("Hunting Area", DESIG_ENG))


shp2.raw <- read_sf("data/spatialData/protectedAreas/rawPAs/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_1/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

shp2 <- shp2.raw %>% 
  filter(IUCN_CAT %in% c("Ia", "Ib", "II") | grepl("ional park", DESIG_ENG) | grepl("ional Park", DESIG_ENG)) %>% 
  filter(!grepl("arine", DESIG_ENG)) %>% 
  filter(!grepl("Hunting Area", DESIG_ENG))


shp3.raw <- read_sf("data/spatialData/protectedAreas/rawPAs/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_2/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

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
mam.to.raw <- read_sf("data/spatialData/mammalRanges/iucn/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

mam.to <- mam.to.raw %>% 
  mutate(nativeness = "IUCN_Native") %>% 
  dplyr::select(sci_name, category, nativeness, geometry)
mapview(mam.to[mam.to$sci_name == "Lynx lynx",])

## erick's ranges ------------

mam.ejl <- read_sf("data/spatialData/mammalRanges/lundgren/Native_Intro_Mammals_Updated_Taxonomy.gpkg") %>%
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

st_write(pas.loop, "data/spatialData/protectedAreas/pasWithMegafaunaColumns.gpkg", append = FALSE)



#### get the species stuff sorted -------------------------


species.cols <- names(pas.loop %>% dplyr::select(-c(unique_id, NAME, WDPA_PID, DESIG_ENG, IUCN_CAT, STATUS_YR, GIS_AREA, geometry, sovereignt, continent)))
species.cols <- species.cols[!grepl("geometry", species.cols)]
dt.pas <- pas.loop %>% 
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

unique(dt.pas[grepl("Lille Vildmose", NAME), ]$species)

unique(dt.pas[grepl("Ilulissat Isfjord", NAME), ]$species)
unique(dt.pas[country == "Denmark", ]$NAME)
unique(dt.pas$country)

## Rangifer tarandus range is kinda bullshit

### calculate functional groups etc... 


# Load megafauna traits and groups----------

mam.traits <- fread("data/processedData/dataFragments/mammal_traits.csv")

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
pas.traits2[is.na(guild_only_herbivory)]$WDPA_PID

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
pas.herbi <- pas.traits2[herbi_or_carni == "Herbivore", ]
pas.herbi[, n_herbi_sp_pa := uniqueN(species, na.rm = TRUE), by = unique_id]
quantile(pas.herbi$n_herbi_sp_pa)



pas.herbi2 <- pas.herbi %>% mutate(
  regular_predator_present = ifelse(mass_kg < largest_accessible_prey_size_pa, 1, 0), 
  occasional_predator_present =  ifelse(mass_kg < maximum_prey_size_pa, 1, 0), 
) %>% filter(!n_herbi_sp_pa < 1) %>% unique()


pas.herbi2[, species_w_occasional_predator := sum(occasional_predator_present), by = unique_id]
pas.herbi2[, species_w_regular_predator := sum(regular_predator_present), by = unique_id]
quantile(pas.herbi2$species_w_regular_predator, na.rm = T )
quantile(pas.herbi2$species_w_occasional_predator, na.rm = T)

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



fwrite(x = dt.final, file = "data/processedData/cleanData/pasWithMegafauna.csv")

