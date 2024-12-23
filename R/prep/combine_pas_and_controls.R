### combine different pas and controls 
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

