## join WDPA 


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
  filter(MARINE == 0) 

nrow(pas.largest[pas.largest$STATUS_YR <= 2013])
nrow(pas.largest[pas.largest$STATUS_YR <= 2003])
