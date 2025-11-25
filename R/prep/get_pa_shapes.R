## join WDPA 


library(sf)
library(tidyverse)
library(data.table)
library(terra)

shp1.raw <- read_sf("data/spatial_data/protected_areas/raw_pas/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_0/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

table(shp1.raw$IUCN_CAT)
shp1 <- shp1.raw %>% 
 # filter(IUCN_CAT %in% c("Ia", "Ib", "II") | grepl("ional park", DESIG_ENG) | grepl("ional Park", DESIG_ENG)) %>% 
  filter(!grepl("arine", DESIG_ENG)) #%>% filter(!grepl("Hunting Area", DESIG_ENG))


shp2.raw <- read_sf("data/spatial_data/protected_areas/raw_pas/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_1/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

shp2 <- shp2.raw %>% 
 # filter(IUCN_CAT %in% c("Ia", "Ib", "II") | grepl("ional park", DESIG_ENG) | grepl("ional Park", DESIG_ENG)) %>% 
  filter(!grepl("arine", DESIG_ENG)) #%>% filter(!grepl("Hunting Area", DESIG_ENG))


shp3.raw <- read_sf("data/spatial_data/protected_areas/raw_pas/WDPA_WDOECM_May2024_Public_all_shp/WDPA_WDOECM_May2024_Public_all_shp_2/WDPA_WDOECM_May2024_Public_all_shp-polygons.shp")

shp3 <- shp3.raw %>% 
 # filter(IUCN_CAT %in% c("Ia", "Ib", "II") | grepl("ional park", DESIG_ENG) | grepl("ional Park", DESIG_ENG)) %>% 
  filter(!grepl("arine", DESIG_ENG)) #%>%  filter(!grepl("Hunting Area", DESIG_ENG))

pas <- rbind(shp1, shp2, shp3)
pas <- pas[!duplicated(data.frame(pas)),]

names(pas)


world <- rnaturalearth::ne_countries() %>% dplyr::select(continent, sovereignt, geometry)

sf_use_s2(FALSE)
pas_cont <- st_join(pas, world)
pas_cont$area_km2 <- st_area(pas_cont)/1000000

set.seed(161)
table(pas_cont$MARINE)
pas_largest <- pas_cont %>% 
  mutate(iucn_cat_ord = case_when(
    IUCN_CAT == "Ia" ~  1,  
    IUCN_CAT == "Ib" ~  2,  
    IUCN_CAT == "II" ~  3,  
    IUCN_CAT == "III" ~  4,  
    IUCN_CAT == "IV" ~  5,  
    IUCN_CAT == "V" ~  6,  
    IUCN_CAT == "VI" ~  7,  
    IUCN_CAT %in% c("Not Applicable", "Not Reported","Not Assigned") ~  8,  
  ), 
  unique_id = paste0(sovereignt, "_", WDPA_PID)) %>% 
  filter(MARINE == 0) %>% 
  group_by(WDPA_PID) %>%
  slice_max(area_km2) %>% 
  rename(country = sovereignt, 
         continent = continent) %>% 
  dplyr::select(WDPAID, WDPA_PID, NAME, unique_id, area_km2, STATUS_YR, iucn_cat_ord, IUCN_CAT, DESIG_ENG, GIS_AREA, continent, country) %>%
  filter(GIS_AREA > 1) %>%   
  group_by(WDPA_PID) %>%
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  mutate(own_pa_id = 1:nrow(.))
  


nrow(pas_largest[pas_largest$STATUS_YR <= 2013,])
nrow(pas_largest[pas_largest$STATUS_YR <= 2003,])


write_sf(pas_largest, "data/spatial_data/protected_areas/pa_shapes.gpkg", append = FALSE)
