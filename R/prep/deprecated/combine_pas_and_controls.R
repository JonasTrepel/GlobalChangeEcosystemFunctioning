### combine different pas and controls 
library(sf)
library(tidyverse)
library(tidylog)


# control areas --------------------------

cont_raw <- read_sf("data/spatialData/protectedAreas/controls_for_strict_pas.gpkg") %>%
  st_transform(4326)

cont <- cont_raw %>%
  mutate(iucn_cat = NA,
         protection_cat_broad = "control", 
         NAME = NA, 
         STATUS_YR = NA, 
         DESIG_ENG = NA, 
         Country = NA, 
         pa_age = NA)


dist <- cont %>%
  as.data.table() %>% 
  mutate(geom = NULL) %>%
  dplyr::select(control_for, control_within_dist) %>% 
  rename(unique_id = control_for)

sf_use_s2(FALSE)
cont$area_km2 <- st_area(cont)/1000000

# pas ----
pas_raw <- read_sf("data/spatialData/protectedAreas/strict_pas.gpkg")

pas <- pas_raw %>%
  rename(iucn_cat = IUCN_CAT, 
         pa_age = PaAge) %>% 
  mutate(protection_cat_broad = "strictly_protected", 
         control_for = NA) %>% 
  filter(unique_id %in% c(unique(cont$control_for))) %>% 
  dplyr::select(c(control_for, Continent, functional_biome, unique_id, geom, 
                  iucn_cat, protection_cat_broad, NAME, STATUS_YR, DESIG_ENG, Country, pa_age, area_km2)) %>% 
  left_join(dist)

setdiff(names(pas), names(cont))
setdiff(names(cont), names(pas))

pa_sf <- rbind(pas, cont)

sf_use_s2(FALSE)
coords <- st_coordinates(st_centroid(pa_sf))

pa_sf$lon_pa <- coords[,1]
pa_sf$lat_pa <- coords[,2]

write_sf(pa_sf, "data/spatialData/pa_and_control_shapes.gpkg")

