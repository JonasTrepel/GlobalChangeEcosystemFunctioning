### get grid subset 

library(terra)
library(sf)
library(data.table)
library(tidyverse)

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

write_sf(grid_biome, "data/spatialData/grid_sample.gpkg")

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
