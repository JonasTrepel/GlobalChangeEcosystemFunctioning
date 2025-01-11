### get grid subset 

library(terra)
library(sf)
library(data.table)
library(tidyverse)


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
grid_sub_1 <- dt_grid %>% 
  filter(fract_prot == 1 | fract_prot == 0) %>% 
  group_by(protection_cat_broad) %>% 
  filter(n() > 50000) %>% 
  sample_n(50000) %>% 
  ungroup()

table(grid_sub_1$protection_cat_broad)

set.seed(161)
grid_sub_2 <- dt_grid %>% 
  filter(fract_prot == 1 | fract_prot == 0) %>% 
  group_by(super_biome) %>% 
  filter(n() > 50000) %>% 
  sample_n(50000) %>% 
  ungroup()

table(grid_sub_2$super_biome)


grid_sub <- rbind(grid_sub_1, grid_sub_2) %>%
  unique() %>% 
  distinct(X, Y, .keep_all = TRUE)

table(grid_sub$protection_cat_broad)
table(grid_sub$super_biome)

keep <- grid_sub %>% dplyr::select(unique_id) %>% pull()


grid_sf <- read_sf("data/spatialData/protectedAreas/gridWithCovs.gpkg") %>% 
  rename(unique_id = gridID) %>% 
  filter(unique_id %in% keep) %>% 
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


write_sf(grid_sf, "data/spatialData/grid_sample.gpkg")

