### Combine analysis ready datasets


library(terra)
library(sf)
library(data.table)
library(tidyverse)
library(tidylog)

#global grid -----------------

dt_grid_covs <- fread("data/processed_data/grid_with_covariates.csv")
glimpse(dt_grid_covs)

dt_grid_trends <- fread("data/processed_data/grid_with_trends.csv") %>% 
  dplyr::select(-c(x_moll, y_moll, lon, lat))

dt_grid <- dt_grid_covs %>% 
  filter(!is.na(functional_biome)) %>% 
  distinct(lon, lat, .keep_all = TRUE) %>% 
  left_join(dt_grid_trends) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "mixed",
      iucn_cat == "unprotected" ~ "unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    )) %>%
  unique() %>% 
  dplyr::select(-c(iucn_cat_ord)) %>% 
  filter(complete.cases(.))

fwrite(dt_grid, "data/processed_data/analysis_ready_world_grid.csv")

#### Strict Protected area dataset -------------------------------------------------
dt_pas_trends <- fread("data/processed_data/pas_with_trends.csv")

# overlap removal message: removed 3387 polgyons (24.2% of the data) - 10583 are left
dt_pas_covs_raw <- fread("data/processed_data/pas_with_covariates.csv") %>% 
  filter(GIS_AREA > 1) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "mixed",
      iucn_cat == "unprotected" ~ "unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    )) %>% 
  filter(protection_cat_broad == "strict") %>% 
  mutate(control_for = "none", 
         pair_id = paste0("pair_", unique_id), 
         control_within_dist = "na") %>% 
  left_join(dt_pas_trends)


dt_pa_controls_treds <- fread("data/processed_data/pa_controls_with_trends.csv")

dt_pa_controls_covs <- fread("data/processed_data/pa_controls_with_covariates.csv") %>% 
  filter(control_for %in% unique(dt_pas_covs_raw$unique_id)) %>%
  mutate(
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    )) %>% 
  mutate(pair_id = paste0("pair_", control_for), 
         WDPA_PID = "na", 
         area_km2 = "na",
         WDPAID = "na", 
         IUCN_CAT = "na",
         NAME = "na", 
         STATUS_YR = "na",
         DESIG_ENG = "na",
         GIS_AREA = "na", 
         country = "na", 
         own_pa_id = "na", 
         iucn_cat = "na", 
         protection_cat_broad = "unprotected") %>% 
  left_join(dt_pa_controls_treds)

nrow(dt_pa_controls_covs)/10583

dt_pas_covs <- dt_pas_covs_raw %>% 
  filter(unique_id %in% unique(dt_pa_controls_covs$control_for))


setdiff(names(dt_pas_covs), names(dt_pa_controls_covs))
setdiff(names(dt_pa_controls_covs), names(dt_pas_covs))


dt_pa_all <- rbind(dt_pas_covs, dt_pa_controls_covs, fill = F) %>% 
  dplyr::select(-c(iucn_cat_ord)) %>% 
  filter(complete.cases(.)) %>% 
  group_by(pair_id) %>% 
  filter(n() == 2) %>% 
  ungroup()
  

fwrite(dt_pa_all, "data/processed_data/analysis_ready_pa_dataset.csv")


dt_pa_all %>% 
  pivot_longer(cols = c(fire_frequency, hmi_change, evi_coef, mat_coef, nitrogen_depo, 
                        max_temp_coef, prec_coef)) %>% 
  ggplot() +
  geom_density_ridges(aes(y = protection_cat_broad, x = value)) +
  facet_wrap(~name, scales = "free")
