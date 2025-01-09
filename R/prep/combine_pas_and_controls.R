### combine different pas and controls 
library(sf)
library(tidyverse)
library(tidylog)


# protected areas --------------------------

# controls 

cont_dt <- fread("data/processedData/cleanData/controls_with_covs.csv")

dist <- cont_dt %>%
  dplyr::select(control_for, control_within_dist) %>% 
  rename(unique_id = control_for)

cont_raw <- read_sf("data/spatialData/protectedAreas/controls_for_strict_pas.gpkg") %>%
  st_transform(4326)

cont <- cont_raw %>%
  left_join(cont_dt) %>%
  mutate(iucn_cat = NA,
         protection_cat_broad = "controls", 
         NAME = NA, 
         STATUS_YR = NA, 
         DESIG_ENG = NA, 
         Country = NA, 
         pa_age = NA, control_within_dist = NA)

sf_use_s2(FALSE)
cont$area_km2 <- st_area(cont)/1000000

# pas
pas_raw <- read_sf("data/spatialData/protectedAreas/strict_pas.gpkg")

pas <- pas_raw %>%
  rename(iucn_cat = IUCN_CAT, 
         pa_age = PaAge) %>% 
  mutate(protection_cat_broad = "strictly_protected", 
         control_for = NA) %>% 
  filter(unique_id %in% c(unique(cont$control_for))) %>% 
  dplyr::select(-c(WDPAID, WDPA_PID, iucnCatOrd, GIS_AREA, seed, own_pa_id)) %>% 
  left_join(dist)

setdiff(names(pas), names(cont))
setdiff(names(cont), names(pas))

pa_sf <- rbind(pas, cont)

write_sf(pa_sf, "data/spatialData/pas_and_controls.gpkg")

st_crs(pas)
