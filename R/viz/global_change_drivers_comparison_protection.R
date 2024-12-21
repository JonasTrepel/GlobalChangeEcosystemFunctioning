library(tidyverse)
library(data.table)
library(remotePARTS)


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
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    ),
    nitrogen_depo = NitrogenDepo,
    mat_coef = mat_coef,
    max_temp_coef = max_temp_coef,
    map_coef = map_coef,
    human_modification = HumanModification, 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = scale(log(PaAge + 0.0001))) 

p_comp <- dt %>% 
  pivot_longer(cols = c(human_modification, nitrogen_depo, mat_coef, max_temp_coef, map_coef), 
               values_to = "global_change_value", names_to = "global_change_driver") %>% 
  mutate(
    clean_driver = case_when(
      global_change_driver == "human_modification" ~ "Human Modification", 
      global_change_driver == "nitrogen_depo" ~ "Nitrogen Deposition", 
      global_change_driver == "mat_coef" ~ "MAT Trend", 
      global_change_driver == "max_temp_coef" ~ "Max Temp Trend", 
      global_change_driver == "map_coef" ~ "MAP Trend"), 
    clean_prot = case_when(
      og_layer == "controls" ~ "Control", 
      og_layer == "protected_areas" ~ "Protected")) %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_jitter(aes(x = clean_prot, y = global_change_value), alpha = 0.25, size = 0.5) +
  geom_boxplot(aes(x = clean_prot, y = global_change_value), outlier.size = 0.25, outlier.shape = NA) +
  facet_wrap(~clean_driver, scales = "free") +
  labs(x = "", y = "Driver Value") +
  theme_bw()
p_comp

ggsave(plot = p_comp, "builds/plots/global_change_drivers_comparison_protection.png", dpi = 600)
