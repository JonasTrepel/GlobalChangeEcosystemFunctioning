library(tidyverse)
library(data.table)
library(remotePARTS)
library(tidylog)
library(scico)
#load trends 

evi_trends <- fread("data/processedData/dataFragments/pa_evi_trends.csv")
burned_area_trends <- fread("data/processedData/dataFragments/pa_burned_area_trends.csv")
greenup_trends <- fread("data/processedData/dataFragments/pa_greenup_trends.csv")


dt <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  group_by(functional_biome) %>% 
  mutate(n_per_functional_biome = n()) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
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
    nitrogen_depo = nitrogen_depo,
    mat_coef = mat_coef,
    max_temp_coef = max_temp_coef,
    map_coef = map_coef,
    human_modification = human_modification, 
    area_km2_log = log(area_km2 + 0.0001),
    pa_age_log = log(pa_age + 0.0001)) %>%
  left_join(evi_trends) %>% 
  left_join(burned_area_trends) %>% 
  left_join(greenup_trends)

dt$protection_cat_broad

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
      protection_cat_broad == "control" ~ "Control", 
      protection_cat_broad == "strictly_protected" ~ "Protected")) %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_jitter(aes(x = clean_prot, y = global_change_value), alpha = 0.25, size = 0.5, color = "grey50") +
  geom_boxplot(aes(x = clean_prot, y = global_change_value), outlier.size = 0.25, outlier.shape = NA, alpha = 0.75) +
  facet_wrap(~clean_driver, scales = "free") +
  labs(x = "", y = "Driver Value") +
  theme_bw()
p_comp

ggsave(plot = p_comp, "builds/plots/global_change_drivers_comparison_protection.png", dpi = 600, height = 4, width = 7)

#### calculate differences 

dt_pa <- dt %>% 
  filter(protection_cat_broad == "strictly_protected")

dt_diff <- data.table()

for(pa_id in unique(dt_pa$unique_id)){
  
  tmp_diff <- data.table()
  
  tmp_diff$evi_diff <- dt[dt$unique_id %in% pa_id, ]$evi_coef - dt[dt$control_for %in% pa_id, ]$evi_coef
  tmp_diff$burned_area_diff <- dt[dt$unique_id %in% pa_id, ]$burned_area_coef - dt[dt$control_for %in% pa_id, ]$burned_area_coef
  tmp_diff$greenup_diff <- dt[dt$unique_id %in% pa_id, ]$greenup_coef - dt[dt$control_for %in% pa_id, ]$greenup_coef
  tmp_diff$human_modification_diff <- dt[dt$unique_id %in% pa_id, ]$human_modification - dt[dt$control_for %in% pa_id, ]$human_modification
  tmp_diff$nitrogen_depo_diff <- dt[dt$unique_id %in% pa_id, ]$nitrogen_depo - dt[dt$control_for %in% pa_id, ]$nitrogen_depo
  tmp_diff$mat_diff <- dt[dt$unique_id %in% pa_id, ]$mat_coef - dt[dt$control_for %in% pa_id, ]$mat_coef
  tmp_diff$max_temp_diff <- dt[dt$unique_id %in% pa_id, ]$max_temp_coef - dt[dt$control_for %in% pa_id, ]$max_temp_coef
  tmp_diff$map_diff <- dt[dt$unique_id %in% pa_id, ]$map_coef - dt[dt$control_for %in% pa_id, ]$map_coef
  
  tmp_diff$unique_id <- pa_id
  tmp_diff$functional_biome <- dt[dt$unique_id %in% pa_id, ]$functional_biome
  tmp_diff$super_biome <- dt[dt$unique_id %in% pa_id, ]$super_biome
  tmp_diff$pa_age <- dt[dt$unique_id %in% pa_id, ]$pa_age
  tmp_diff$pa_age_log <- dt[dt$unique_id %in% pa_id, ]$pa_age_log
  tmp_diff$area_km2_log <- dt[dt$unique_id %in% pa_id, ]$area_km2_log
  tmp_diff$area_km2 <- dt[dt$unique_id %in% pa_id, ]$area_km2
  

  dt_diff <- rbind(tmp_diff, dt_diff) 
  
  print(paste0(pa_id, " done"))
}

p_diff <- dt_diff %>%
  mutate(clean_biome = case_when(
    super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
    super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
    super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
    super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation"
  )) %>% 
  pivot_longer(cols = c("evi_diff", "greenup_diff", "burned_area_diff",
                        "human_modification_diff", "nitrogen_depo_diff", 
                        "mat_diff", "map_diff", "max_temp_diff"), 
               names_to = "diff_name", values_to = "diff_value") %>% 
  mutate(clean_diff = case_when(
    diff_name == "evi_diff" ~ "EVI\nTrend Diff.", 
    diff_name == "greenup_diff" ~ "Greenup\nTrend Diff.", 
    diff_name == "burned_area_diff" ~ "Burned Area\nTrend Diff.",
    diff_name == "human_modification_diff" ~ "Human Modification\nDiff.", 
    diff_name == "nitrogen_depo_diff" ~ "Nitrogen Depo\n Diff.", 
    diff_name == "mat_diff" ~ "MAT\nTrend Diff.",
    diff_name == "max_temp_diff" ~ "Max Temp\nTrend Diff.", 
    diff_name == "map_diff" ~ "MAP\nTrend Diff."
  )) %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_jitter(aes(x = clean_biome, y = diff_value, color = clean_biome), alpha = 0.25, size = 0.5, width = 0.2) +
  geom_boxplot(aes(x = clean_biome, y = diff_value), outlier.size = 0.25, outlier.shape = NA, alpha = 0.75) + 
  scale_color_scico_d("batlowK") +
  facet_wrap(~clean_diff, scales = "free", ncol = 4) +
  theme_bw() +
  labs(x = "", y = "", subtitle = "positive values indicate variable in PA > variable in control") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        legend.position = "none")
p_diff

ggsave(plot = p_diff, "builds/plots/driver_diff_pas.png", dpi = 600)


summary(lm(evi_diff ~ 1, data = dt_diff)); AIC(lm(evi_diff ~ 1, data = dt_diff))#66521.04
summary(lm(evi_diff ~ pa_age, data = dt_diff)); AIC(lm(evi_diff ~ pa_age, data = dt_diff)) #63816.09
summary(lm(evi_diff ~ pa_age_log, data = dt_diff)); AIC(lm(evi_diff ~ pa_age_log, data = dt_diff)) #63820.1
summary(lm(evi_diff ~ area_km2, data = dt_diff)); AIC(lm(evi_diff ~ area_km2, data = dt_diff)) #66522.99
summary(lm(evi_diff ~ area_km2_log, data = dt_diff)); AIC(lm(evi_diff ~ area_km2_log, data = dt_diff)) #66521.53
summary(lm(evi_diff ~ scale(area_km2_log) + scale(pa_age_log), data = dt_diff)); AIC(lm(evi_diff ~ scale(area_km2_log) + scale(pa_age_log), data = dt_diff)) #63820.64
summary(lm(evi_diff ~ scale(area_km2) + scale(pa_age), data = dt_diff)); AIC(lm(evi_diff ~ scale(area_km2) + scale(pa_age), data = dt_diff)) #63818.04
#best is PA age unlogged, second best is logged 

summary(lm(burned_area_diff ~ 1, data = dt_diff)); AIC(lm(burned_area_diff ~ 1, data = dt_diff))#-65565.77
summary(lm(burned_area_diff ~ pa_age, data = dt_diff)); AIC(lm(burned_area_diff ~ pa_age, data = dt_diff)) #-62799.76
summary(lm(burned_area_diff ~ pa_age_log, data = dt_diff)); AIC(lm(burned_area_diff ~ pa_age_log, data = dt_diff)) #-62797.45
summary(lm(burned_area_diff ~ area_km2, data = dt_diff)); AIC(lm(burned_area_diff ~ area_km2, data = dt_diff)) #-65563.79
summary(lm(burned_area_diff ~ area_km2_log, data = dt_diff)); AIC(lm(burned_area_diff ~ area_km2_log, data = dt_diff)) #-65565.38
summary(lm(burned_area_diff ~ scale(area_km2_log) + scale(pa_age_log), data = dt_diff)); AIC(lm(burned_area_diff ~ scale(area_km2_log) + scale(pa_age_log), data = dt_diff)) #-62796.11
summary(lm(burned_area_diff ~ scale(area_km2) + scale(pa_age), data = dt_diff)); AIC(lm(burned_area_diff ~ scale(area_km2) + scale(pa_age), data = dt_diff)) #-62797.9
#best seems to be intercept only 

summary(lm(greenup_diff ~ 1, data = dt_diff)); AIC(lm(greenup_diff ~ 1, data = dt_diff))#18449.42
summary(lm(greenup_diff ~ pa_age, data = dt_diff)); AIC(lm(greenup_diff ~ pa_age, data = dt_diff)) #17799.91
summary(lm(greenup_diff ~ pa_age_log, data = dt_diff)); AIC(lm(greenup_diff ~ pa_age_log, data = dt_diff)) #17798.39
summary(lm(greenup_diff ~ area_km2, data = dt_diff)); AIC(lm(greenup_diff ~ area_km2, data = dt_diff)) #18451.34
summary(lm(greenup_diff ~ area_km2_log, data = dt_diff)); AIC(lm(greenup_diff ~ area_km2_log, data = dt_diff)) #18451.23
summary(lm(greenup_diff ~ scale(area_km2_log) + scale(pa_age_log), data = dt_diff)); AIC(lm(greenup_diff ~ scale(area_km2_log) + scale(pa_age_log), data = dt_diff)) #17800.39
summary(lm(greenup_diff ~ scale(area_km2) + scale(pa_age), data = dt_diff)); AIC(lm(greenup_diff ~ scale(area_km2) + scale(pa_age), data = dt_diff)) #17801.85
#best is PA logged


dt_diff %>%
  mutate(clean_biome = case_when(
    super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
    super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
    super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
    super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation"
  )) %>% 
  pivot_longer(cols = c("evi_diff", "greenup_diff", "burned_area_diff"), 
               names_to = "diff_name", values_to = "diff_value") %>% 
  mutate(clean_diff = case_when(
    diff_name == "evi_diff" ~ "EVI\nTrend Diff.", 
    diff_name == "greenup_diff" ~ "Greenup\nTrend Diff.", 
    diff_name == "burned_area_diff" ~ "Burned Area\nTrend Diff."
  )) %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_point(aes(x = pa_age, y = diff_value), alpha = 0.25, size = 0.5, color = "grey50", width = 0.2) +
  geom_smooth(aes(x = pa_age, y = diff_value), method = "lm") +
  facet_wrap(~clean_diff, scales = "free") +
  theme_bw() +
  labs(x = "PA age (years)", y = "") +
  theme(axis.text.x = element_text(angle = 0))  

dt_diff %>%
  mutate(clean_biome = case_when(
    super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
    super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
    super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
    super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation"
  )) %>% 
  pivot_longer(cols = c("evi_diff", "greenup_diff", "burned_area_diff"), 
               names_to = "diff_name", values_to = "diff_value") %>% 
  mutate(clean_diff = case_when(
    diff_name == "evi_diff" ~ "EVI\nTrend Diff.", 
    diff_name == "greenup_diff" ~ "Greenup\nTrend Diff.", 
    diff_name == "burned_area_diff" ~ "Burned Area\nTrend Diff."
  )) %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_point(aes(x = area_km2, y = diff_value), alpha = 0.25, size = 0.5, color = "grey50", width = 0.2) +
  geom_smooth(aes(x = area_km2, y = diff_value), method = "lm") +
  facet_wrap(~clean_diff, scales = "free") +
  theme_bw() +
  labs(x = "PA area (km^2)", y = "") +
  #scale_x_log10() +
  theme(axis.text.x = element_text(angle = 0))  

