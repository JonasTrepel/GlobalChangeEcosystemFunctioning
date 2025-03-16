### get grid subset 

library(terra)
library(sf)
library(data.table)
library(tidyverse)
library(tidylog)

grid_sf <- read_sf("data/spatialData/protectedAreas/gridWithCovs.gpkg")

dt_grid_raw <- fread("data/processedData/cleanData/gridWithCovs.csv")

dt_grid <- dt_grid_raw %>% 
  filter(!is.na(functional_biome)) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "Mixed",
      iucn_cat == "unprotected" ~ "Unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    )) %>%
  rename(unique_id = gridID)  %>%
  unique() %>% 
  distinct(X, Y, .keep_all = TRUE)


columns_to_check <- c("functional_biome", "X", "Y", "lon", "lat", "unique_id",
                      "nitrogen_depo", "human_modification", "super_biome", "protection_cat_broad")
dt_grid <- dt_grid[complete.cases(dt_grid[, ..columns_to_check])]

table(dt_grid$protection_cat_broad)
table(dt_grid$super_biome)

set.seed(161)
grid_sub_biome <- dt_grid %>% 
  filter(fract_prot == 1 | fract_prot == 0) %>% 
  group_by(super_biome) %>% 
  filter(n() > 100000) %>% 
  sample_n(100000) %>% 
  ungroup()

table(grid_sub_biome$super_biome)
table(grid_sub_biome$protection_cat_broad)

keep_biome <- grid_sub_biome %>% dplyr::select(unique_id) %>% pull()

grid_biome <- grid_sf %>% 
  rename(unique_id = gridID) %>% 
  filter(unique_id %in% keep_biome) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "Mixed",
      iucn_cat == "unprotected" ~ "Unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    ))

write_sf(grid_biome, "data/spatialData/grid_sample.gpkg", append = FALSE)

### USA subset -----------------------
library(mapview)

#outline from https://public.opendatasoft.com/explore/dataset/us-state-boundaries/export/
usa_all <- st_read("data/rawData/raw_time_series/usa_n_depo/usa_boundaries/us-state-boundaries.shp")
mapview(usa_all)
st_crs(usa_all)


usa <- usa_all %>% 
  filter(name %in% c(
    "Alabama", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", 
    "Delaware", "Florida", "Georgia", "Idaho", "Illinois", "Indiana", "Iowa", 
    "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", 
    "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", 
    "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York", 
    "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", 
    "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", 
    "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", 
    "West Virginia", "Wisconsin", "Wyoming"
  )) %>% summarize() %>% st_transform(crs = st_crs(grid_sf))

st_crs(grid_sf) == st_crs(usa)

grid_usa <- grid_sf %>% 
  left_join(dt_grid) %>% 
  filter(fract_prot == 1 | fract_prot == 0) %>% 
  filter(lengths(st_intersects(., usa)) > 0) %>% 
  sample_n(100000) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "Mixed",
      iucn_cat == "unprotected" ~ "Unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    ),
    unique_id = NULL) %>% 
  rename(unique_id = gridID)

#mapview::mapview(grid_usa)

write_sf(grid_usa, "data/spatialData/grid_usa.gpkg", append = FALSE)


### Europe subset ---------------------------

r <- rast("data/rawData/raw_time_series/europe_n_depo/europe_n_depo_2020.tif")

extent_shp <- as.polygons(ext(r), crs = crs(r)) %>% 
  st_as_sf() %>% 
  st_transform(crs = st_crs(grid_sf))

set.seed(161)
grid_europe <- grid_sf %>% 
  left_join(dt_grid) %>% 
  filter(fract_prot == 1 | fract_prot == 0) %>% 
  filter(lengths(st_intersects(., extent_shp)) > 0) %>% 
  sample_n(100000) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "Mixed",
      iucn_cat == "unprotected" ~ "Unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    ),
    unique_id = NULL) %>% 
  rename(unique_id = gridID)

#mapview::mapview(grid_europe)

write_sf(grid_europe, "data/spatialData/grid_europe.gpkg", append = FALSE)


### Protection cat subset (if wanted)
set.seed(161)
grid_sub_prot <- dt_grid %>% 
  filter(fract_prot == 1 | fract_prot == 0) %>% 
  group_by(protection_cat_broad) %>% 
  filter(n() > 50000) %>% 
  sample_n(50000) %>% 
  ungroup()

table(grid_sub_prot$protection_cat_broad)

keep_prot <- grid_sub_biome %>% dplyr::select(unique_id) %>% pull()

grid_sub_prot <- grid_sf %>% 
  rename(unique_id = gridID) %>% 
  filter(unique_id %in% keep_prot) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "Mixed",
      iucn_cat == "unprotected" ~ "Unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    ))

write_sf(grid_sub_prot, "data/spatialData/grid_sub_prot_sample.gpkg")


### Protection cat subset (if wanted)
set.seed(161)
grid_sub_prot <- dt_grid %>% 
  filter(fract_prot == 1 | fract_prot == 0) %>% 
  group_by(protection_cat_broad) %>% 
  filter(n() > 50000) %>% 
  sample_n(50000) %>% 
  ungroup()

table(grid_sub_prot$protection_cat_broad)

keep_prot <- grid_sub_biome %>% dplyr::select(unique_id) %>% pull()

grid_sub_prot <- grid_sf %>% 
  rename(unique_id = gridID) %>% 
  filter(unique_id %in% keep_prot) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "Mixed",
      iucn_cat == "unprotected" ~ "Unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    ))

write_sf(grid_sub_prot, "data/spatialData/grid_sub_prot_sample.gpkg")
