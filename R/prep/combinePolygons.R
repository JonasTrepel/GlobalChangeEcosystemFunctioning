### combine different polygons to create a shapefile of layers of interest
library(sf)
library(tidyverse)
# protected areas 

pas_raw <- read_sf("data/spatialData/protectedAreas/paShapes.gpkg")

pas <- pas_raw %>%
  rename(iucnCat = IUCN_CAT) %>% 
  dplyr::select(unique_id , iucnCat, geom) %>% 
  mutate(og_layer = "protected_areas", prot_cat_biome = NA)

# controls 
    
cont_raw <- read_sf("data/spatialData/protectedAreas/controlsForPas.gpkg")

cont <- cont_raw %>%
  mutate(iucnCat = NA, 
         unique_id = paste0(controlFor, "_control")) %>% 
  dplyr::select(unique_id , iucnCat, geom)  %>% 
  mutate(og_layer = "controls", prot_cat_biome = NA)

## grids

grid_raw <- read_sf("data/spatialData/protectedAreas/gridWithCovs.gpkg") %>% 
  filter(!is.na(FunctionalBiomeShort)) %>% 
  dplyr::select(-Country) 

set.seed(161)
grid_int <- grid_raw %>%
  group_by(gridID) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  mutate(iucnCat = ifelse(is.na(iucnCat), "NotProtected", iucnCat), 
         protection_cat_broad = case_when(
           iucnCat %in% c("Ia", "Ib", "II") ~ "Strict", 
           iucnCat %in% c("III", "IV", "V", "VI", "UnknownOrNA") ~ "Mixed",
           iucnCat == "NotProtected" ~ "Unprotected"
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
  dplyr::select(unique_id, iucnCat, geom, prot_cat_biome) %>% 
  mutate(og_layer = "grid")

#combine
roi <- rbind(pas, cont, grid)

mapview::mapview(grid, zcol = "prot_cat_biome")

write_sf(roi, "data/spatialData/regions_of_interest.gpkg")
