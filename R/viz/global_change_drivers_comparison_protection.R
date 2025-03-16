library(tidyverse)
library(data.table)
library(remotePARTS)
library(tidylog)
library(scico)
library(ggridges)
#load trends 


evi_trend <- fread("data/processedData/dataFragments/pas_mean_evi_trends.csv") %>% 
  mutate(unique_pa_id = gsub("_[^_]+$", "", unique_id)) %>% 
  group_by(unique_pa_id) %>% 
  summarize(mean_evi_coef = median(mean_evi_coef), 
            abs_mean_evi_coef = median(abs_mean_evi_coef))
burned_area_trend <- fread("data/processedData/dataFragments/pas_burned_area_trends.csv") %>% 
  mutate(unique_pa_id = gsub("_[^_]+$", "", unique_id)) %>% 
  group_by(unique_pa_id) %>% 
  summarize(burned_area_coef = median(burned_area_coef), 
            abs_burned_area_coef = median(abs_burned_area_coef))

greenup_trend <- fread("data/processedData/dataFragments/pas_greenup_trends.csv") %>% 
  mutate(unique_pa_id = gsub("_[^_]+$", "", unique_id)) %>% 
  group_by(unique_pa_id) %>% 
  summarize(greenup_coef = median(greenup_coef), 
            abs_greenup_coef = median(abs_greenup_coef))

get_mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

f_biome <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  dplyr::select(unique_pa_id, functional_biome) %>%
  group_by(unique_pa_id) %>% 
  summarize(functional_biome = get_mode(functional_biome)) %>% 
  rename(unique_id = unique_pa_id) %>% 
  mutate(
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    )
  )

#### include super biome!
dt_sum <- climate_trends %>%
  group_by(unique_pa_id) %>% 
  summarize(mat_coef = mean(mat_coef, na.rm = T),
            map_coef = mean(map_coef, na.rm = T),
            max_temp_coef = mean(max_temp_coef, na.rm = T),
            human_modification = mean(human_modification, na.rm = T),
            nitrogen_depo = mean(nitrogen_depo, na.rm = T),
            protection_cat_broad = unique(protection_cat_broad), 
            lon_pa = unique(lon_pa), 
            lat_pa = unique(lat_pa)) %>% 
  left_join(evi_trend) %>% 
  left_join(burned_area_trend) %>% 
  left_join(greenup_trend) %>% 
  mutate(mean_evi_coef = mean_evi_coef /100, 
         burned_area_coef = burned_area_coef*100) %>% 
  rename(unique_id = unique_pa_id) 


dt <- dt_sum %>% mutate(
  control_for = gsub("_control", "", unique_id)
) %>% left_join(f_biome)

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
  geom_density_ridges(aes(y = clean_prot, x = global_change_value, color = clean_prot, fill = clean_prot), alpha = 0.75, size = 0.5) +
 # geom_boxplot(aes(x = clean_prot, y = global_change_value), outlier.size = 0.25, outlier.shape = NA, alpha = 0.75) +
  facet_wrap(~clean_driver, scales = "free", ncol = 5) +
  scale_color_scico_d(palette = "bamako", begin = 0.25, end = 0.75) +
  scale_fill_scico_d(palette = "bamako", begin = 0.25, end = 0.75) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "none")
p_comp

ggsave(plot = p_comp, "builds/plots/global_change_drivers_comparison_protection.png", dpi = 600, height = 3, width = 10)

#### calculate differences 

dt_pa <- dt %>% 
  filter(protection_cat_broad == "strictly_protected")

dt_diff <- data.table()

for(pa_id in unique(dt_pa$unique_id)){
  
  tmp_diff <- data.table()
  
  tmp_diff$evi_diff <- dt[dt$unique_id %in% pa_id, ]$mean_evi_coef - dt[dt$control_for %in% pa_id, ]$mean_evi_coef
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

## plot. 

q_025_evi <- as.numeric(quantile(dt_diff$evi_diff, c(.025), na.rm = T)) 
q_975_evi <- as.numeric(quantile(dt_diff$evi_diff, c(.975), na.rm = T))
q_025_greenup <- as.numeric(quantile(dt_diff$greenup_diff, c(.025), na.rm = T)) 
q_975_greenup <- as.numeric(quantile(dt_diff$greenup_diff, c(.975), na.rm = T))
q_025_burned_area <- as.numeric(quantile(dt_diff$burned_area_diff, c(.025), na.rm = T)) 
q_975_burned_area <- as.numeric(quantile(dt_diff$burned_area_diff, c(.975), na.rm = T))
q_025_human_modification <- as.numeric(quantile(dt_diff$human_modification_diff, c(.025), na.rm = T)) 
q_975_human_modification <- as.numeric(quantile(dt_diff$human_modification_diff, c(.975), na.rm = T))
q_025_nitrogen_depo <- as.numeric(quantile(dt_diff$nitrogen_depo_diff, c(.025), na.rm = T)) 
q_975_nitrogen_depo <- as.numeric(quantile(dt_diff$nitrogen_depo_diff, c(.975), na.rm = T))
q_025_mat <- as.numeric(quantile(dt_diff$mat_diff, c(.025), na.rm = T)) 
q_975_mat <- as.numeric(quantile(dt_diff$mat_diff, c(.975), na.rm = T))
q_025_map <- as.numeric(quantile(dt_diff$map_diff, c(.025), na.rm = T)) 
q_975_map <- as.numeric(quantile(dt_diff$map_diff, c(.975), na.rm = T))
q_025_max_temp <- as.numeric(quantile(dt_diff$max_temp_diff, c(.025), na.rm = T)) 
q_975_max_temp <- as.numeric(quantile(dt_diff$max_temp_diff, c(.975), na.rm = T))



p_diff <- dt_diff %>%
  mutate(clean_biome = case_when(
    super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
    super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
    super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
    super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation"),
    evi_diff = ifelse(evi_diff > q_975_evi, NA, evi_diff),
    evi_diff = ifelse(evi_diff < q_025_evi, NA, evi_diff),
    burned_area_diff = ifelse(burned_area_diff > q_975_burned_area, NA, burned_area_diff),
    burned_area_diff = ifelse(burned_area_diff < q_025_burned_area, NA, burned_area_diff),
    greenup_diff = ifelse(greenup_diff > q_975_greenup, NA, greenup_diff),
    greenup_diff = ifelse(greenup_diff < q_025_greenup, NA, greenup_diff),
    human_modification_diff = ifelse(human_modification_diff > q_975_human_modification, NA, human_modification_diff),
    human_modification_diff = ifelse(human_modification_diff < q_025_human_modification, NA, human_modification_diff),
    nitrogen_depo_diff = ifelse(nitrogen_depo_diff > q_975_nitrogen_depo, NA, nitrogen_depo_diff),
    nitrogen_depo_diff = ifelse(nitrogen_depo_diff < q_025_nitrogen_depo, NA, nitrogen_depo_diff),
    mat_diff = ifelse(mat_diff > q_975_mat, NA, mat_diff),
    mat_diff = ifelse(mat_diff < q_025_mat, NA, mat_diff),
    map_diff = ifelse(map_diff > q_975_map, NA, map_diff),
    map_diff = ifelse(map_diff < q_025_map, NA, map_diff),
    max_temp_diff = ifelse(max_temp_diff > q_975_max_temp, NA, max_temp_diff),
    max_temp_diff = ifelse(max_temp_diff < q_025_max_temp, NA, max_temp_diff)) %>% 
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
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey25") +
  geom_density_ridges(aes(y = clean_biome, x = diff_value, color = clean_biome, fill = clean_biome), alpha = 0.5, size = 0.5, width = 0.2) +
  #geom_boxplot(aes(x = clean_biome, y = diff_value), outlier.size = 0.25, outlier.shape = NA, alpha = 0.75) + 
  scale_color_scico_d("batlowK") +
  scale_fill_scico_d("batlowK") +
  facet_wrap(~clean_diff, scales = "free_x", ncol = 4) +
  theme_minimal() +
  labs(x = "", y = "", subtitle = "positive values indicate variable in PA > variable in control") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "none")
p_diff

ggsave(plot = p_diff, "builds/plots/driver_diff_pas.png", dpi = 600, height = 6, width = 9)


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


p_diff_pa_age <- dt_diff %>%
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
  geom_point(aes(x = pa_age, y = diff_value), alpha = 0.25, size = 0.5, color = "grey50", width = 0.2) +
  geom_smooth(aes(x = pa_age, y = diff_value), method = "lm") +
  facet_wrap(~clean_diff, scales = "free_y", ncol = 4) +
  theme_minimal() +
  scale_x_log10() +
  labs(x = "PA age (years)\n(log-scale)", y = "", title = "a)") +
  theme(axis.text.x = element_text(angle = 0))  
p_diff_pa_age

p_diff_pa_area <- dt_diff %>%
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
  geom_point(aes(x = area_km2, y = diff_value), alpha = 0.25, size = 0.5, color = "grey50", width = 0.2) +
  geom_smooth(aes(x = area_km2, y = diff_value), method = "lm") +
  facet_wrap(~clean_diff, scales = "free_y", ncol = 4) +
  scale_x_log10() +
  theme_minimal() +
  labs(x = "PA area (km^2)\n(log-scale)", y = "", title = "b)") +
  #scale_x_log10() +
  theme(axis.text.x = element_text(angle = 0))  
p_diff_pa_area


p <- gridExtra::grid.arrange(p_diff_pa_age, p_diff_pa_area)
ggsave(plot = p, "builds/plots/pas_diff_vs_area_and_age.png", dpi = 600, height = 10, width = 10)
