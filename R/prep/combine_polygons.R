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

set.seed(161)
grid_int <- grid_raw %>%
  rename(iucn_cat = iucnCat) %>% 
  group_by(gridID) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  mutate(iucn_cat = ifelse(is.na(iucn_cat), "NotProtected", iucn_cat), 
         protection_cat_broad = case_when(
           iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
           iucn_cat %in% c("III", "IV", "V", "VI", "UnknownOrNA") ~ "Mixed",
           iucn_cat == "NotProtected" ~ "Unprotected"
         ), 
         prot_cat_biome = paste0(protection_cat_broad, "_", FunctionalBiomeShort)) %>% 
  dplyr::select(-unique_id) %>%
  group_by(prot_cat_biome) %>% 
  filter(n() >= 2500) %>% 
  sample_n(2500) %>% 
  ungroup() 

coords <- grid_int %>% st_centroid() %>% st_coordinates()

grid_int$X <- coords[,1]
grid_int$Y <- coords[,2]

write_sf(grid_int, "data/spatialData/grid_sample.gpkg")






grid <- grid_int %>% 
  rename(unique_id = gridID) %>% 
  dplyr::select(unique_id, iucn_cat, geom, prot_cat_biome) %>% 
  mutate(og_layer = "grid")



#combine
roi <- rbind(pas %>% dplyr::select(unique_id , iucn_cat, geom, og_layer) %>% mutate(prot_cat_biome = NA),
             cont %>% dplyr::select(unique_id , iucn_cat, geom, og_layer) %>% mutate(prot_cat_biome = NA),
             grid)

mapview::mapview(grid, zcol = "prot_cat_biome")

write_sf(roi, "data/spatialData/regions_of_interest.gpkg")
