
library(tidyverse)
library(data.table)

evi_trend <- fread("data/processedData/dataFragments/pa_evi_trends.csv")
burned_area_trend <- fread("data/processedData/dataFragments/pa_burned_area_trends.csv")
greenup_trend <- fread("data/processedData/dataFragments/pa_greenup_trends.csv")


dt <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  group_by(FunctionalBiome) %>% 
  mutate(n_per_functional_biome = n()) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  # filter(n_per_functional_biome > 1000) %>% 
  rename(functional_biome = FunctionalBiome) %>% 
  mutate(
    productivity = case_when(
      grepl("L", functional_biome) ~ "low", 
      grepl("M", functional_biome) ~ "medium", 
      grepl("H", functional_biome) ~ "high"
    ), 
    ndvi_min = case_when(
      grepl("C", functional_biome) ~ "cold", 
      grepl("D", functional_biome) ~ "dry", 
      grepl("B", functional_biome) ~ "cold_and_dry", 
      grepl("N", functional_biome) ~ "non_seasonal"
    ), 
    nitrogen_depo = scale(NitrogenDepo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(HumanModification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = scale(log(PaAge + 0.0001))) %>% 
  left_join(evi_trend) %>% 
  left_join(burned_area_trend) %>% 
  left_join(greenup_trend)


#select only pas

dt_pa <- dt %>% filter(og_layer == "protected_areas")

dt_diff <- data.table()

for(i in 1:nrow(dt_pa)){
  
  pa_id <- unique(dt_pa[i,]$unique_id)
  
  dt_pa_sub <- dt_pa[i,]
  
  dt_cont_sub <- dt %>% filter(control_for == pa_id)
  
  dt_tmp <- data.table(
    evi_coef_diff = dt_pa_sub$evi_coef - dt_cont_sub$evi_coef,
    burned_area_coef_diff = dt_pa_sub$burned_area_coef - dt_cont_sub$burned_area_coef,
    greenup_coef_diff = dt_pa_sub$greenup_coef - dt_cont_sub$greenup_coef, 
    functional_biome = dt_pa_sub$functional_biome, 
    pa_age = dt_pa_sub$PaAge, 
    pa_age_log = dt_pa_sub$pa_age_log, 
    area_km2 = dt_pa_sub$area_km2,
    area_km2_log = dt_pa_sub$area_km2_log,
    productivity = dt_pa_sub$productivity, 
    ndvi_min = dt_pa_sub$ndvi_min
    )
  
  dt_diff <- rbind(dt_diff, dt_tmp)
  
  print(paste0(i, "/", nrow(dt_pa)))
}

names(dt_diff) <- gsub(".V1", "", names(dt_diff))



dt_long <- dt_diff %>% 
  pivot_longer(cols = c(evi_coef_diff, burned_area_coef_diff, greenup_coef_diff), names_to = "trend_name", values_to = "trend_value")


dt_long %>%
  group_by(trend_name) %>%
  mutate(trend_value = ifelse(
    trend_value > quantile(trend_value, 0.975, na.rm = TRUE) | 
      trend_value < quantile(trend_value, 0.025, na.rm = TRUE), 
    NA, trend_value)) %>%
  ungroup() %>%
  ggplot() +
  geom_boxplot(aes(x = productivity, y = trend_value)) +
  facet_wrap(~trend_name, scales = "free")


dt_long %>%
  group_by(trend_name) %>%
  mutate(trend_value = ifelse(
    trend_value > quantile(trend_value, 0.975, na.rm = TRUE) | 
      trend_value < quantile(trend_value, 0.025, na.rm = TRUE), 
    NA, trend_value)) %>%
  ungroup() %>%
  ggplot() +
  geom_boxplot(aes(x = ndvi_min, y = trend_value)) +
  facet_wrap(~trend_name, scales = "free")


dt_long %>%
  group_by(trend_name) %>%
  mutate(trend_value = ifelse(
    trend_value > quantile(trend_value, 0.975, na.rm = TRUE) | 
      trend_value < quantile(trend_value, 0.025, na.rm = TRUE), 
    NA, trend_value)) %>%
  ungroup() %>%
  filter(pa_age_log > -4) %>% 
  ggplot() +
  geom_point(aes(x = pa_age_log, y = trend_value)) +
  geom_smooth(aes(x = pa_age_log, y = trend_value)) +
  facet_wrap(~trend_name, scales = "free")

dt_long %>%
  group_by(trend_name) %>%
  mutate(trend_value = ifelse(
    trend_value > quantile(trend_value, 0.975, na.rm = TRUE) | 
      trend_value < quantile(trend_value, 0.025, na.rm = TRUE), 
    NA, trend_value)) %>%
  ungroup() %>%
 # filter(pa_age_log > -4) %>% 
  ggplot() +
  geom_point(aes(x = area_km2_log, y = trend_value)) +
  geom_smooth(aes(x = area_km2_log, y = trend_value)) +
  facet_wrap(~trend_name, scales = "free")
