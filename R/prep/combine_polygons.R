### combine different polygons to create a shapefile of layers of interest
library(sf)
library(tidyverse)


# protected areas --------------------------

# controls 

cont_dt <- fread("data/processedData/cleanData/controls_with_covs.csv")

cont_raw <- read_sf("data/spatialData/protectedAreas/controls_for_strict_pas.gpkg")

cont <- cont_raw %>%
  left_join(cont_dt) %>%
  mutate(iucn_cat = NA,
         og_layer = "controls", 
         NAME = NA, 
         STATUS_YR = NA, 
         DESIG_ENG = NA, 
         Country = NA, 
         PaAge = NA)

sf_use_s2(FALSE)
cont$area_km2 <- st_area(cont)/1000000

# pas
pas_raw <- read_sf("data/spatialData/protectedAreas/strict_pas.gpkg")

pas <- pas_raw %>%
  rename(iucn_cat = IUCN_CAT) %>% 
  mutate(og_layer = "protected_areas", 
         control_for = NA) %>% 
  filter(unique_id %in% c(unique(cont$control_for))) %>% 
  dplyr::select(-c(WDPAID, WDPA_PID, iucnCatOrd, GIS_AREA, seed))

setdiff(names(pas), names(cont))
setdiff(names(cont), names(pas))

pa_sf <- rbind(pas, cont)

write_sf(pa_sf, "data/spatialData/pas_and_controls.gpkg")


## grid ----------------

grid_raw <- read_sf("data/spatialData/protectedAreas/gridWithCovs.gpkg") %>% 
  filter(!is.na(FunctionalBiomeShort)) %>% 
  dplyr::select(-Country) 

coords <- grid_raw %>% st_centroid() %>% st_coordinates()

grid_raw$X <- coords[,1]
grid_raw$Y <- coords[,2]

set.seed(161)
grid_int <- grid_raw %>%
  rename(iucn_cat = iucnCat) %>% 
  group_by(gridID) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  #group_by(X, Y) %>% 
  #slice_sample(n = 1) %>% 
  #ungroup() %>% 
  mutate(productivity = case_when(
    grepl("L", FunctionalBiomeShort) ~ "low", 
    grepl("M", FunctionalBiomeShort) ~ "medium", 
    grepl("H", FunctionalBiomeShort) ~ "high"
    ), 
    ndvi_min = case_when(
    grepl("C", FunctionalBiomeShort) ~ "cold", 
    grepl("D", FunctionalBiomeShort) ~ "dry", 
    grepl("B", FunctionalBiomeShort) ~ "cold_and_dry", 
    grepl("N", FunctionalBiomeShort) ~ "non_seasonal"
    ), 
    iucn_cat = ifelse(is.na(iucn_cat), "NotProtected", iucn_cat), 
         protection_cat_broad = case_when(
           iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
           iucn_cat %in% c("III", "IV", "V", "VI", "UnknownOrNA") ~ "Mixed",
           iucn_cat == "NotProtected" ~ "Unprotected"
    ), 
    super_biome = case_when(
      grepl("C", FunctionalBiomeShort) & grepl("T", FunctionalBiomeShort) ~ "cold_tall", 
      grepl("C", FunctionalBiomeShort) & grepl("S", FunctionalBiomeShort) ~ "cold_short", 
      !grepl("C", FunctionalBiomeShort) & grepl("T", FunctionalBiomeShort) ~ "not_cold_tall", 
      !grepl("C", FunctionalBiomeShort) & grepl("S", FunctionalBiomeShort) ~ "not_cold_short"
    )
    ) %>% 
  dplyr::select(-unique_id) %>%
  rename(unique_id = gridID)

grid_sub_1 <- grid_int %>% 
  group_by(protection_cat_broad) %>% 
  filter(n() > 25000) %>% 
  sample_n(25000) %>% 
  ungroup()

grid_sub_2 <- grid_int %>% 
  group_by(super_biome) %>% 
  filter(n() > 25000) %>% 
  sample_n(25000) %>% 
  ungroup()



grid_sub <- rbind(grid_sub_1, grid_sub_2) %>%
  unique() %>% 
  distinct(X, Y, .keep_all = TRUE)


#coords <- grid_sub %>% st_centroid() %>% st_coordinates()

#grid_sub$X <- coords[,1]
#grid_sub$Y <- coords[,2]

write_sf(grid_sub, "data/spatialData/grid_sample.gpkg")

