### make nice maps ####

library(tidyverse)
library(data.table)
library(RColorBrewer)
library(scico)
library(rnaturalearth)
library(sf)
library(gridExtra)
library(tidylog)
library(RColorBrewer)

evi_trend <- fread("data/processedData/dataFragments/grid_mean_evi_trends.csv") 
burned_area_trend <- fread("data/processedData/dataFragments/grid_burned_area_trends.csv")
greenup_trend <- fread("data/processedData/dataFragments/grid_greenup_trends.csv")

climate_trends <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") %>% 
  dplyr::select(unique_id, mat_coef, map_coef, max_temp_coef, mat_p_value, map_p_value, max_temp_p_value)

sf_use_s2(FALSE)
world <- rnaturalearth::ne_countries() %>% filter(!name_en == "Antarctica") %>%
  st_transform(crs = 'ESRI:54030') %>% 
  mutate(geometry = st_make_valid(geometry), 
         world = "world") %>% 
  group_by(world) %>% 
  summarize()


raw_shapes <- read_sf("data/spatialData/grid_sample.gpkg")

shapes <- raw_shapes %>%
  st_transform(crs = 'ESRI:54030') %>% 
  left_join(evi_trend) %>% 
  left_join(burned_area_trend) %>% 
  left_join(greenup_trend) %>% 
  left_join(climate_trends) %>% 
  mutate(mean_evi_coef = mean_evi_coef /100, 
         burned_area_coef = burned_area_coef*100)



sum(is.na(shapes$mean_evi_coef))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
##################################     MAPS    #######################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

## evi maps -------------
quantile(shapes$mean_evi_coef, c(.025, .975), na.rm = T)

q_025_evi <- as.numeric(quantile(shapes$mean_evi_coef, c(.025), na.rm = T))
q_975_evi <- as.numeric(quantile(shapes$mean_evi_coef, c(.975), na.rm = T))

p_evi_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = shapes %>%
            filter(!is.na(mean_evi_coef)) %>% 
            mutate(mean_evi_coef = ifelse(mean_evi_coef > q_975_evi, q_975_evi, mean_evi_coef),
                   mean_evi_coef = ifelse(mean_evi_coef < q_025_evi, q_025_evi, mean_evi_coef)),
          aes(color = mean_evi_coef, fill = mean_evi_coef)) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  labs(color = "EVI\nTrend", fill = "EVI\nTrend") +
  theme_void() +
  theme(legend.position = "right", 
        # legend.key.width = unit(1, "cm"),
        # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 0))
p_evi_shapes

ggsave(plot = p_evi_shapes, "builds/plots/evi_grid_shapes_map.png", dpi = 600)

## burned area maps-------
quantile(shapes$burned_area_coef, c(.025, .975), na.rm = T)

q_025_burned_area <- as.numeric(quantile(shapes$burned_area_coef, c(.025), na.rm = T))
q_975_burned_area <- as.numeric(quantile(shapes$burned_area_coef, c(.975), na.rm = T))

p_burned_area_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = shapes %>%
            filter(!is.na(burned_area_coef)) %>% 
            mutate(burned_area_coef = ifelse(burned_area_coef > q_975_burned_area, q_975_burned_area, burned_area_coef),
                   burned_area_coef = ifelse(burned_area_coef < q_025_burned_area, q_025_burned_area, burned_area_coef)),
          aes(color = burned_area_coef, fill = burned_area_coef)) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  labs(color = "Burned\nArea\nTrend", fill = "Burned\nArea\nTrend") +
  theme_void() +
  theme(legend.position = "right", 
        # legend.key.width = unit(1, "cm"),
        # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 0))
p_burned_area_shapes

ggsave(plot = p_burned_area_shapes, "builds/plots/burned_area_grid_shapes_map.png", dpi = 600)

## greenup maps ---------------------
quantile(shapes$greenup_coef, c(.025, .975), na.rm = T)

q_025_greenup <- as.numeric(quantile(shapes$greenup_coef, c(.025), na.rm = T))
q_975_greenup <- as.numeric(quantile(shapes$greenup_coef, c(.975), na.rm = T))

p_greenup_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = shapes %>%
            filter(!is.na(greenup_coef)) %>% 
            mutate(greenup_coef = ifelse(greenup_coef > q_975_greenup, q_975_greenup, greenup_coef),
                   greenup_coef = ifelse(greenup_coef < q_025_greenup, q_025_greenup, greenup_coef)),
          aes(color = greenup_coef, fill = greenup_coef)) +
  scale_color_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  labs(color = "Vegetation\nGreen-Up\nTrend", fill = "Vegetation\nGreen-Up\nTrend") +
  theme_void() +
  theme(legend.position = "right", 
        # legend.key.width = unit(1, "cm"),
        # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 0))
p_greenup_shapes

ggsave(plot = p_greenup_shapes, "builds/plots/greenup_grid_shapes_map.png", dpi = 600)

p_all_trends_map_shapes <- gridExtra::grid.arrange(p_evi_shapes, p_burned_area_shapes, p_greenup_shapes, ncol = 1)
ggsave(plot = p_all_trends_map_shapes, "builds/plots/grid_all_trends_map_shapes.png", dpi = 600, height = 10, width = 8)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################     ESTIMATES     ###################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

evi_est <- fread("builds/model_estimates/envi_grid_estimates.csv") %>% mutate(response = "evi", 
                                                                             estimate = estimate/100, 
                                                                             std_error = std_error/100, 
                                                                             ci_lb = ci_lb/100, 
                                                                             ci_ub = ci_ub/100)
burned_area_est <- fread("builds/model_estimates/burned_area_grid_estimates.csv") %>% mutate(response = "burned_area", 
                                                                                             estimate = estimate*100, 
                                                                                             std_error = std_error*100, 
                                                                                             ci_lb = ci_lb*100, 
                                                                                             ci_ub = ci_ub*100)
greenup_est <- fread("builds/model_estimates/greenup_grid_estimates.csv") %>% mutate(response = "greenup")

dt_est <- rbind(evi_est, burned_area_est, greenup_est) %>% 
  as.data.table() %>%
  mutate(clean_term = case_when(
    grepl("Intercept", term) ~ "Intercept", 
    grepl("nitrogen_depo", term) ~ "Nitrogen Deposition",
    grepl("mat_coef", term) ~ "MAT Trend", 
    grepl("map_coef", term) ~ "MAP Trend",
    grepl("max_temp_coef", term) ~ "Max Temp Trend",
    grepl("human_modification", term) ~ "Human Modification",
    term == "super_biomecold_short" ~ "Cold Limited\nShort Vegetation",
    term == "super_biomecold_tall" ~ "Cold Limited\nTall Vegetation",
    term == "super_biomenot_cold_short" ~ "Not Cold Limited\nShort Vegetation",
    term == "super_biomenot_cold_tall" ~ "Not Cold Limited\nTall Vegetation",
    term == "protection_cat_broadStrict" ~ "Strict", 
    term == "protection_cat_broadMixed" ~ "Mixed", 
    term == "protection_cat_broadUnprotected" ~ "Unprotected", 
    term == "area_km2_log" ~ "PA Area",
    term == "pa_age_log" ~ "PA Age"),
    sig_pn = case_when(
      estimate < 0 & p_value < 0.05 ~ "Sig. Negative", 
      estimate > 0 & p_value < 0.05 ~ "Sig. Positive", 
      p_value >= 0.05 ~ "Non-Significant"
    ),
    facet_label = case_when(
      model == "H1" ~ "H1: Trend ~\nIntercept", 
      model == "H2" ~ "Global dataset", 
      model == "H3" ~ "H3: Trend ~\nBiome", 
      model == "H4" ~ "H4: Trend ~\nProtection", 
      model == "H4.1" ~ "H4.1: Change ~\nProtection", 
      model == "H5" ~ "H5: Trend ~\nPA Characteristics", 
      model == "H5.1" ~ "H5.1: Change ~\nPA Characteristics", 
      model == "cold_short" ~ "Cold Limited\nShort Vegetation",
      model == "cold_tall" ~ "Cold Limited\nTall Vegetation",
      model == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
      model == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation",
    )
  ) 

fwrite(dt_est %>% 
         dplyr::select(-facet_label, -sig, -sig_pn, -term, -t.stat, -std_error) %>% 
         filter(!model == "H4.1" & !model == "H5.1") %>% 
         mutate(p_value = round(p_value, 4), 
                estimate = round(estimate, 3), 
                ci_lb = round(ci_lb, 3),
                ci_ub = round(ci_ub, 3), 
                dataset = "Grid") %>% 
         dplyr::select(Response = response, Model = model,
                       Variable = clean_term, Estimate = estimate,
                       `lower CI` = ci_lb, `upper CI` = ci_ub, p = p_value, dataset), 
       "builds/model_estimates/combined_clean_estimates_grid.csv")

# Evi estimates -----
scico(palette = "bam", n = 10)
#"#65014B" "#9E3C85" "#C86FB1" "#E4ADD6" "#F4E3EF" "#EFF3E5" "#C0D9A1" "#7BA755" "#457B2A" "#0C4C00"

p_est_evi <- dt_est %>%
  mutate(clean_term = factor(clean_term, levels = c(
    "Intercept", 
    "MAP Trend", 
    "Max Temp Trend", 
    "MAT Trend", 
    "Human Modification",
    "Nitrogen Deposition"))) %>%
  mutate(facet_label = factor(facet_label, levels = c(
    "Global dataset",
    "Cold Limited\nShort Vegetation",
    "Cold Limited\nTall Vegetation",
    "Not Cold Limited\nShort Vegetation",
    "Not Cold Limited\nTall Vegetation"
  ))) %>% 
  filter(response == "evi" & model %in% c("H2", "cold_tall", "not_cold_tall", "cold_short", "not_cold_short")) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free_x", ncol = 5) +
  scale_color_manual(values = c("Sig. Negative" = "#9E3C85",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "b)", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_est_evi

fig_evi <- grid.arrange(p_evi_shapes + labs(title = "a)"), 
                        p_est_evi,
                        heights = c(1.3, 1))
ggsave(plot = fig_evi, "builds/plots/grid_evi_figure.png", height = 8, width = 9.5, dpi = 600)

# burned area estimates -----
scico(palette = "vik", n = 10)
#"#001260" "#023E7D" "#1D6E9C" "#71A7C4" "#C9DDE7" "#EACEBE" "#D29773" "#BD6432" "#8B2706" "#590007"

p_est_burned_area <- dt_est %>%
  mutate(clean_term = factor(clean_term, levels = c(
    "Intercept", 
    "MAP Trend", 
    "Max Temp Trend", 
    "MAT Trend", 
    "Human Modification",
    "Nitrogen Deposition"))) %>%
  mutate(facet_label = factor(facet_label, levels = c(
    "Global dataset",
    "Cold Limited\nShort Vegetation",
    "Cold Limited\nTall Vegetation",
    "Not Cold Limited\nShort Vegetation",
    "Not Cold Limited\nTall Vegetation"
  ))) %>% 
  filter(response == "burned_area" & model %in% c("H2", "cold_tall", "not_cold_tall", "cold_short", "not_cold_short")) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free_x", ncol = 5) +
  scale_color_manual(values = c("Sig. Negative" = "#023E7D",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#8B2706"
  )) +
  labs(title = "b)", y = NULL, x = "Burned Area Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 22.5, hjust = 1),
        plot.title = element_text(size = 12))

p_est_burned_area

fig_burned_area <- grid.arrange(p_burned_area_shapes + labs(title = "a)"), 
                        p_est_burned_area,
                        heights = c(1.3, 1))
ggsave(plot = fig_burned_area, "builds/plots/grid_burned_area_figure.png", height = 8, width = 9.5, dpi = 600)


# Greenup estimates -----
scico(palette = "cork", n = 10, direction = -1, begin = 0.1, end = 0.9)
#"#195715" "#3D7D3C" "#6C9C6B" "#A2C0A1" "#D8E5D9" "#D2DDE7" "#97AFC8" "#6388AD" "#386695" "#284074"

p_est_greenup <- dt_est %>%
  mutate(clean_term = factor(clean_term, levels = c(
    "Intercept", 
    "MAP Trend", 
    "Max Temp Trend", 
    "MAT Trend", 
    "Human Modification",
    "Nitrogen Deposition"))) %>%
  mutate(facet_label = factor(facet_label, levels = c(
    "Global dataset",
    "Cold Limited\nShort Vegetation",
    "Cold Limited\nTall Vegetation",
    "Not Cold Limited\nShort Vegetation",
    "Not Cold Limited\nTall Vegetation"
  ))) %>% 
  filter(response == "greenup" & model %in% c("H2", "cold_tall", "not_cold_tall", "cold_short", "not_cold_short")) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free_x", ncol = 5) +
  scale_color_manual(values = c("Sig. Negative" = "#3D7D3C",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#386695"
  )) +
  labs(title = "b)", y = NULL, x = "Vegetation Green-Up Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
       # axis.text.x = element_text(angle = 22.5, hjust = 1),
        plot.title = element_text(size = 12))

p_est_greenup


fig_greenup <- grid.arrange(p_greenup_shapes + labs(title = "a)"), 
                            p_est_greenup, 
                            heights = c(1.3, 1))
ggsave(plot = fig_greenup, "builds/plots/grid_greenup_figure.png", height = 8, width = 9.5, dpi = 600)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################     CORRELATION     ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

dt_corr <- shapes %>% 
  left_join(climate_trends) %>% 
  mutate(pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA)) %>%
  rename(pa_area = area_km2) %>% 
  as.data.table() %>% 
  dplyr::select(
    mean_evi_coef, burned_area_coef, greenup_coef,
    mat_coef, map_coef, max_temp_coef, 
    nitrogen_depo, human_modification, super_biome) %>% 
  filter(complete.cases(.)) %>% 
  rename(`EVI Trend` = mean_evi_coef, 
         `Burned Area Trend` = burned_area_coef, 
         `Greenup Trend` = greenup_coef, 
         `MAT Trend` = mat_coef,
         `MAP Trend` = map_coef, 
         `Max Temp Trend` = max_temp_coef, 
         `Nitrogen Deposition` = nitrogen_depo, 
         `Human Modification` = human_modification)

dt_corr_full <- dt_corr %>% dplyr::select(-super_biome)

library(ggcorrplot)
corr <- round(cor(dt_corr_full), 1)
p_corr_full <- ggcorrplot(corr, hc.order = TRUE, type = "lower",
                     lab = TRUE)
p_corr_full
#ggsave(plot = p_corr_full, "builds/plots/grid_variable_correlations.png", dpi = 600, height = 8, width = 8)


#### Biome specific correlations ---
## cold_short
dt_corr_cold_short <- dt_corr %>%
  filter(super_biome == "cold_short") %>%
  dplyr::select(-super_biome)

library(ggcorrplot)
corr_cold_short <- round(cor(dt_corr_cold_short), 1)
p_corr_cold_short <- ggcorrplot(corr_cold_short, hc.order = TRUE, type = "lower",
                     lab = TRUE)
p_corr_cold_short

## cold_tall
dt_corr_cold_tall <- dt_corr %>%
  filter(super_biome == "cold_tall") %>%
  dplyr::select(-super_biome)

library(ggcorrplot)
corr_cold_tall <- round(cor(dt_corr_cold_tall), 1)
p_corr_cold_tall <- ggcorrplot(corr_cold_tall, hc.order = TRUE, type = "lower",
                                lab = TRUE)
p_corr_cold_tall

## not_cold_short
dt_corr_not_cold_short <- dt_corr %>%
  filter(super_biome == "not_cold_short") %>%
  dplyr::select(-super_biome)

library(ggcorrplot)
corr_not_cold_short <- round(cor(dt_corr_not_cold_short), 1)
p_corr_not_cold_short <- ggcorrplot(corr_not_cold_short, hc.order = TRUE, type = "lower",
                                lab = TRUE)
p_corr_not_cold_short

## not_cold_tall
dt_corr_not_cold_tall <- dt_corr %>%
  filter(super_biome == "not_cold_tall") %>%
  dplyr::select(-super_biome)

library(ggcorrplot)
corr_not_cold_tall <- round(cor(dt_corr_not_cold_tall), 1)
p_corr_not_cold_tall <- ggcorrplot(corr_not_cold_tall, hc.order = TRUE, type = "lower",
                                lab = TRUE)
p_corr_not_cold_tall

p_biome_corr <- gridExtra::grid.arrange(
  p_corr_cold_short + labs(title = "Cold limited, short vegetation"),
  p_corr_cold_tall + labs(title = "Cold limited, tall vegetation"),
  p_corr_not_cold_short + labs(title = "Not cold limited, short vegetation"),
  p_corr_not_cold_tall + labs(title = "Not cold limited, tall vegetation")
)

#ggsave(plot = p_biome_corr, "builds/plots/grid_variable_correlations_biome_spec.png", dpi = 600, height = 12, width = 12)

p_corr <- gridExtra::grid.arrange(
  p_corr_full + labs(title = "Full Dataset"),
  p_corr_cold_short + labs(title = "Cold limited, short vegetation"),
  p_corr_cold_tall + labs(title = "Cold limited, tall vegetation"),
  p_corr_not_cold_short + labs(title = "Not cold limited, short vegetation"),
  p_corr_not_cold_tall + labs(title = "Not cold limited, tall vegetation")
)

ggsave(plot = p_corr, "builds/plots/grid_variable_correlations_all.png", dpi = 600, height = 13, width = 11)


### correlations in strictly protected areas 
dt_corr_pa <- shapes %>% 
  left_join(climate_trends) %>% 
  mutate(pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA)) %>%
  rename(pa_area = area_km2) %>% 
  as.data.table() %>% 
  dplyr::select(
    mean_evi_coef, burned_area_coef, greenup_coef,
    mat_coef, map_coef, max_temp_coef, 
    nitrogen_depo, human_modification, super_biome, pa_age, pa_area) %>% 
  filter(complete.cases(.)) %>% 
  rename(`EVI Trend` = mean_evi_coef, 
         `Burned Area Trend` = burned_area_coef, 
         `Greenup Trend` = greenup_coef, 
         `MAT Trend` = mat_coef,
         `MAP Trend` = map_coef, 
         `Max Temp Trend` = max_temp_coef, 
         `Nitrogen Deposition` = nitrogen_depo, 
         `Human Modification` = human_modification, 
         `PA Age` = pa_age, 
         `PA Area` = pa_area)

dt_corr_pa_full <- dt_corr_pa %>% dplyr::select(-super_biome)

library(ggcorrplot)
corr_pa <- round(cor(dt_corr_pa_full), 1)
p_corr_pa_full <- ggcorrplot(corr_pa, hc.order = TRUE, type = "lower",
                     lab = TRUE)
p_corr_pa_full
#ggsave(plot = p_corr_pa_full, "builds/plots/grid_variable_correlations_strict_pas.png", dpi = 600, height = 8, width = 8)


#### Biome specific correlations ---
## cold_short
dt_corr_pa_cold_short <- dt_corr_pa %>%
  filter(super_biome == "cold_short") %>%
  dplyr::select(-super_biome)

library(ggcorrplot)
corr_pa_cold_short <- round(cor(dt_corr_pa_cold_short), 1)
p_corr_pa_cold_short <- ggcorrplot(corr_pa_cold_short, hc.order = TRUE, type = "lower",
                                lab = TRUE)
p_corr_pa_cold_short

## cold_tall
dt_corr_pa_cold_tall <- dt_corr_pa %>%
  filter(super_biome == "cold_tall") %>%
  dplyr::select(-super_biome)

library(ggcorrplot)
corr_pa_cold_tall <- round(cor(dt_corr_pa_cold_tall), 1)
p_corr_pa_cold_tall <- ggcorrplot(corr_pa_cold_tall, hc.order = TRUE, type = "lower",
                               lab = TRUE)
p_corr_pa_cold_tall

## not_cold_short
dt_corr_pa_not_cold_short <- dt_corr_pa %>%
  filter(super_biome == "not_cold_short") %>%
  dplyr::select(-super_biome)

library(ggcorrplot)
corr_pa_not_cold_short <- round(cor(dt_corr_pa_not_cold_short), 1)
p_corr_pa_not_cold_short <- ggcorrplot(corr_pa_not_cold_short, hc.order = TRUE, type = "lower",
                                    lab = TRUE)
p_corr_pa_not_cold_short

## not_cold_tall
dt_corr_pa_not_cold_tall <- dt_corr_pa %>%
  filter(super_biome == "not_cold_tall") %>%
  dplyr::select(-super_biome)

library(ggcorrplot)
corr_pa_not_cold_tall <- round(cor(dt_corr_pa_not_cold_tall), 1)
p_corr_pa_not_cold_tall <- ggcorrplot(corr_pa_not_cold_tall, hc.order = TRUE, type = "lower",
                                   lab = TRUE)
p_corr_pa_not_cold_tall

dt_biome_corr_pa <- gridExtra::grid.arrange(
  p_corr_pa_cold_short + labs(title = "Cold limited, short vegetation"),
  p_corr_pa_cold_tall + labs(title = "Cold limited, tall vegetation"),
  p_corr_pa_not_cold_short + labs(title = "Not cold limited, short vegetation"),
  p_corr_pa_not_cold_tall + labs(title = "Not cold limited, tall vegetation")
)


#ggsave(plot = dt_biome_corr_pa, "builds/plots/grid_variable_correlations_biome_spec_strict_pas.png", dpi = 600, height = 12, width = 12)

p_corr_pa <- gridExtra::grid.arrange(
  p_corr_pa_full + labs(title = "Full Dataset"),
  p_corr_pa_cold_short + labs(title = "Cold limited, short vegetation"),
  p_corr_pa_cold_tall + labs(title = "Cold limited, tall vegetation"),
  p_corr_pa_not_cold_short + labs(title = "Not cold limited, short vegetation"),
  p_corr_pa_not_cold_tall + labs(title = "Not cold limited, tall vegetation")
)

ggsave(plot = p_corr_pa, "builds/plots/grid_variable_correlations_strict_pas_all.png", dpi = 600, height = 13, width = 11)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
############################     TREND DISTRIBUTION     ##############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(ggridges)

dt_ridges <- shapes %>% 
  as.data.table() %>%
  mutate(geom = NULL, 
         pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
         area_km2_log = scale(log(area_km2 + 0.0001)),
         pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA),
         mean_evi_coef = ifelse(mean_evi_coef > q_975_evi, NA, mean_evi_coef),
         mean_evi_coef = ifelse(mean_evi_coef < q_025_evi, NA, mean_evi_coef), 
         burned_area_coef = ifelse(burned_area_coef > q_975_burned_area, NA, burned_area_coef),
         burned_area_coef = ifelse(burned_area_coef < q_025_burned_area, NA, burned_area_coef), 
         greenup_coef = ifelse(greenup_coef > q_975_greenup, NA, greenup_coef),
         greenup_coef = ifelse(greenup_coef < q_025_greenup, NA, greenup_coef), 
         protection_status = case_when(
           protection_cat_broad == "Strict" ~ "Strict",
           protection_cat_broad == "Unprotected" ~ "Unprotected",
           protection_cat_broad == "Mixed" ~ "Mixed"), 
         biome_clean = case_when(
           super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
           super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
           super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
           super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation"))

#evi ridges
p_evi_prot <- ggplot() +
  geom_density_ridges_gradient(data = dt_ridges %>%
                                 filter(!is.na(mean_evi_coef)),
                               aes(x = mean_evi_coef, y = protection_status, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(y = "", x = "EVI Trend") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_evi_prot

p_evi_biome <- ggplot() +
  geom_density_ridges_gradient(data = dt_ridges %>%
                                 filter(!is.na(mean_evi_coef)),
                               aes(x = mean_evi_coef, y = biome_clean, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(y = "", x = "EVI Trend") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_evi_biome

p_evi_dens <- gridExtra::grid.arrange(p_evi_prot, p_evi_biome, ncol = 1)


# burned area ridges
p_burned_area_prot <- ggplot() +
  geom_density_ridges_gradient(data = dt_ridges %>%
                                 filter(!is.na(burned_area_coef)),
                               aes(x = burned_area_coef, y = protection_status, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(y = "", x = "Burned Area Trend") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_burned_area_prot

p_burned_area_biome <- ggplot() +
  geom_density_ridges_gradient(data = dt_ridges %>%
                                 filter(!is.na(burned_area_coef)),
                               aes(x = burned_area_coef, y = biome_clean, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(y = "", x = "Burned Area Trend") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_burned_area_biome

p_burned_area_dens <- gridExtra::grid.arrange(p_burned_area_prot, p_burned_area_biome, ncol = 1)

## greenup ridges 
p_greenup_prot <- ggplot() +
  geom_density_ridges_gradient(data = dt_ridges %>%
                                 filter(!is.na(greenup_coef)),
                               aes(x = greenup_coef, y = protection_status, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(y = "", x = "Vegetation Green-Up Trend") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_greenup_prot

p_greenup_biome <- ggplot() +
  geom_density_ridges_gradient(data = dt_ridges %>%
                                 filter(!is.na(greenup_coef)),
                               aes(x = greenup_coef, y = biome_clean, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(y = "", x = "Vegetation Green-Up Trend") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_greenup_biome

p_greenup_dens <- gridExtra::grid.arrange(p_greenup_prot, p_greenup_biome, ncol = 1)

p_ridges <- gridExtra::grid.arrange(p_evi_dens, p_burned_area_dens, p_greenup_dens, ncol = 3)
ggsave(plot = p_ridges, "builds/plots/grid_ridges.png", dpi = 600, height = 6, width = 12)


p_pa_age <- ggplot() +
  geom_density_ridges(data = dt_ridges,
                               aes(x = pa_age_log, y = biome_clean, fill = biome_clean, color = biome_clean), alpha = 0.7) +
  scale_color_scico_d("batlowK") +
  scale_fill_scico_d("batlowK") +
  theme_bw() +
  labs(y = "", x = "PA Age (log)", title = "a)") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_pa_age

p_pa_area <- ggplot() +
  geom_density_ridges(data = dt_ridges,
                               aes(x = area_km2_log, y = biome_clean, fill = biome_clean, color = biome_clean), alpha = 0.7) +
  scale_color_scico_d("batlowK") +
  scale_fill_scico_d("batlowK") +
  theme_bw() +
  labs(y = "", x = "PA Area (log)", title = "b)") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_pa_area

p_pa_area_trend <- grid.arrange(p_pa_age, p_pa_area, ncol = 2)
ggsave(plot = p_pa_area_trend, "builds/plots/grid_pa_area_and_age_dist.png", dpi = 600, height = 4, width = 8)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################     BIOME MAP     ####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

p_biome_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = shapes %>% 
            filter(!is.na(mean_evi_coef)) %>%
            mutate(biome_clean = case_when(
              super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
              super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
              super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
              super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation")),
          aes(color = biome_clean, fill = biome_clean), alpha = 1) +
  scale_color_scico_d(palette = "batlowK") +
  scale_fill_scico_d(palette = "batlowK") +
  theme_void() +
  labs(color = "Biome", fill = "Biome") +
  theme(axis.title = element_blank())
p_biome_shapes

ggsave(plot = p_biome_shapes, "builds/plots/grid_biome_maps_shapes.png", dpi = 600)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   PROTECTION MAP   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

p_pa_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = shapes %>% 
            filter(!is.na(mean_evi_coef)),
          aes(color = protection_cat_broad, fill = protection_cat_broad), alpha = 1) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  theme_void() +
  labs(color = "Protection\nStatus", fill = "Protection\nStatus") +
  theme(axis.title = element_blank())
p_pa_shapes

ggsave(plot = p_pa_shapes, "builds/plots/grid_protection_map_shapes.png", dpi = 600)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
############################     GLOBAL CHANGE MAPS     ##############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
dt_env_raw <- shapes %>% 
  left_join(climate_trends) 

q_025_mat <- as.numeric(quantile(dt_env_raw$mat_coef, c(.025), na.rm = T))
q_975_mat <- as.numeric(quantile(dt_env_raw$mat_coef, c(.975), na.rm = T))
q_025_map <- as.numeric(quantile(dt_env_raw$map_coef, c(.025), na.rm = T))
q_975_map <- as.numeric(quantile(dt_env_raw$map_coef, c(.975), na.rm = T))
q_025_max_temp <- as.numeric(quantile(dt_env_raw$max_temp_coef, c(.025), na.rm = T))
q_975_max_temp <- as.numeric(quantile(dt_env_raw$max_temp_coef, c(.975), na.rm = T))
q_025_nitrogen_depo <- as.numeric(quantile(dt_env_raw$nitrogen_depo, c(.025), na.rm = T))
q_975_nitrogen_depo <- as.numeric(quantile(dt_env_raw$nitrogen_depo, c(.975), na.rm = T))
q_025_human_modification <- as.numeric(quantile(dt_env_raw$human_modification, c(.025), na.rm = T))
q_975_human_modification <- as.numeric(quantile(dt_env_raw$human_modification, c(.975), na.rm = T))

dt_env <-  dt_env_raw %>%
  mutate(mat_coef = ifelse(mat_coef > q_975_mat, q_975_mat, mat_coef),
         mat_coef = ifelse(mat_coef < q_025_mat, q_025_mat, mat_coef),
         map_coef = ifelse(map_coef > q_975_map, q_975_map, map_coef),
         map_coef = ifelse(map_coef < q_025_map, q_025_map, map_coef),
         max_temp_coef = ifelse(max_temp_coef > q_975_max_temp, q_975_max_temp, max_temp_coef),
         max_temp_coef = ifelse(max_temp_coef < q_025_max_temp, q_025_max_temp, max_temp_coef),
         nitrogen_depo = ifelse(nitrogen_depo > q_975_nitrogen_depo, q_975_nitrogen_depo, nitrogen_depo),
         nitrogen_depo = ifelse(nitrogen_depo < q_025_nitrogen_depo, q_025_nitrogen_depo, nitrogen_depo),
         human_modification = ifelse(human_modification > q_975_human_modification, q_975_human_modification, human_modification),
         human_modification = ifelse(human_modification < q_025_human_modification, q_025_human_modification, human_modification)
  )

p_mat_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)),
          aes(color = mat_coef, fill = mat_coef), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "MAT\nTrend", fill = "MAT\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_mat_shapes

ggsave(plot = p_mat_shapes, "builds/plots/grid_mat_map_shapes.png", dpi = 600)

p_map_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)),
          aes(color = map_coef, fill = map_coef), alpha = 1) +
  scale_color_scico(palette = "broc", midpoint = 0, direction = -1, end = 0.9, begin = 0.1) +
  scale_fill_scico(palette = "broc", midpoint = 0, direction = -1, end = 0.9, begin = 0.1) +
  theme_void() +
  labs(color = "MAP\nTrend", fill = "MAP\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_map_shapes

ggsave(plot = p_map_shapes, "builds/plots/grid_map_map_shapes.png", dpi = 600)


p_max_temp_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)),
          aes(color = max_temp_coef, fill = max_temp_coef), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "Max\nTemp\nTrend", fill = "Max\nTemp\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_max_temp_shapes

ggsave(plot = p_max_temp_shapes, "builds/plots/grid_max_temp_map_shapes.png", dpi = 600)


p_n_depo_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)),
          aes(color = nitrogen_depo, fill = nitrogen_depo), alpha = 1) +
  scale_color_scico(palette = "batlow") +
  scale_fill_scico(palette = "batlow") +
  theme_void() +
  labs(color = "Nitrogen\nDeposition", fill = "Nitrogen\nDeposition") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_n_depo_shapes

ggsave(plot = p_n_depo_shapes, "builds/plots/grid_n_depo_map_shapes.png", dpi = 600)

p_human_modification_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)),
          aes(color = human_modification, fill = human_modification), alpha = 1) +
  scale_color_scico(palette = "bamako") +
  scale_fill_scico(palette = "bamako") +
  theme_void() +
  labs(color = "Human\nModification", fill = "Human\nModification") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_human_modification_shapes

ggsave(plot = p_human_modification_shapes, "builds/plots/grid_human_modification_map_shapes.png", dpi = 600)

gc_maps <- grid.arrange(p_n_depo_shapes, p_human_modification_shapes, 
                         p_mat_shapes, p_max_temp_shapes, 
                         p_map_shapes, ncol = 2)

ggsave(plot = gc_maps, "builds/plots/grid_global_change_map_shapes.png", dpi = 600)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#################################     QUANTILES     ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

dt_env_raw <- shapes %>% 
  left_join(climate_trends) 



quantile(dt_env_raw[!is.na(dt_env_raw$mean_evi_coef),]$mat_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#-0.045576069  0.004695984  0.012743080  0.022046219  0.029451839  0.040541147  0.079413342
quantile(dt_env_raw[!is.na(dt_env_raw$mean_evi_coef),]$max_temp_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#-0.081249289 -0.004047223  0.008281520  0.017152691  0.025692758  0.041528039  0.219463404 
quantile(dt_env_raw[!is.na(dt_env_raw$mean_evi_coef),]$map_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#-50.77784118  -4.59155461  -0.07938528   0.47506487   1.21686432   4.28346692  43.26811849 
quantile(dt_env_raw[!is.na(dt_env_raw$mean_evi_coef),]$nitrogen_depo, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#        0%         5%        25%        50%        75%        95%       100% 
#   8.93624   38.26200  108.45508  215.97137  446.29657 1136.52124 6848.16553 
quantile(dt_env_raw[!is.na(dt_env_raw$mean_evi_coef),]$human_modification, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#0.000000e+00 1.911390e-06 5.226702e-03 5.079506e-02 1.635810e-01 4.231209e-01 9.732789e-01 
quantile(dt_env_raw[!is.na(dt_env_raw$mean_evi_coef),]$mean_evi_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#-163.8214168  -14.7408541   -0.7717181    4.5138879   13.1196674   29.2534735  180.7633833 
quantile(dt_env_raw[!is.na(dt_env_raw$mean_evi_coef),]$burned_area_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#-0.049467534 -0.004512144  0.000000000  0.000000000  0.000000000  0.002043215  0.048151842 
quantile(dt_env_raw[!is.na(dt_env_raw$mean_evi_coef),]$greenup_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#           0%            5%           25%           50%           75%           95%          100%
#-11.590909090  -0.922264706  -0.246465032  -0.009843419   0.238165137   1.726039631  11.247583560 


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#########################     SIG. GLOBAL CHANGE MAPS     ############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
dt_env_raw <- shapes %>% 
  left_join(climate_trends)  

q_025_mat <- as.numeric(quantile(dt_env_raw$mat_coef, c(.025), na.rm = T))
q_975_mat <- as.numeric(quantile(dt_env_raw$mat_coef, c(.975), na.rm = T))
q_025_map <- as.numeric(quantile(dt_env_raw$map_coef, c(.025), na.rm = T))
q_975_map <- as.numeric(quantile(dt_env_raw$map_coef, c(.975), na.rm = T))
q_025_max_temp <- as.numeric(quantile(dt_env_raw$max_temp_coef, c(.025), na.rm = T))
q_975_max_temp <- as.numeric(quantile(dt_env_raw$max_temp_coef, c(.975), na.rm = T))
q_025_nitrogen_depo <- as.numeric(quantile(dt_env_raw$nitrogen_depo, c(.025), na.rm = T))
q_975_nitrogen_depo <- as.numeric(quantile(dt_env_raw$nitrogen_depo, c(.975), na.rm = T))
q_025_human_modification <- as.numeric(quantile(dt_env_raw$human_modification, c(.025), na.rm = T))
q_975_human_modification <- as.numeric(quantile(dt_env_raw$human_modification, c(.975), na.rm = T))
q_025_evi <- as.numeric(quantile(dt_env_raw$mean_evi_coef, c(.025), na.rm = T))
q_975_evi <- as.numeric(quantile(dt_env_raw$mean_evi_coef, c(.975), na.rm = T))
q_025_burned_area <- as.numeric(quantile(dt_env_raw$burned_area_coef, c(.025), na.rm = T))
q_975_burned_area <- as.numeric(quantile(dt_env_raw$burned_area_coef, c(.975), na.rm = T))
q_025_greenup <- as.numeric(quantile(dt_env_raw$greenup_coef, c(.025), na.rm = T))
q_975_greenup <- as.numeric(quantile(dt_env_raw$greenup_coef, c(.975), na.rm = T))

dt_env <-  dt_env_raw %>%
  mutate(mat_coef = ifelse(mat_coef > q_975_mat, q_975_mat, mat_coef),
         mat_coef = ifelse(mat_coef < q_025_mat, q_025_mat, mat_coef),
         map_coef = ifelse(map_coef > q_975_map, q_975_map, map_coef),
         map_coef = ifelse(map_coef < q_025_map, q_025_map, map_coef),
         max_temp_coef = ifelse(max_temp_coef > q_975_max_temp, q_975_max_temp, max_temp_coef),
         max_temp_coef = ifelse(max_temp_coef < q_025_max_temp, q_025_max_temp, max_temp_coef),
         nitrogen_depo = ifelse(nitrogen_depo > q_975_nitrogen_depo, q_975_nitrogen_depo, nitrogen_depo),
         nitrogen_depo = ifelse(nitrogen_depo < q_025_nitrogen_depo, q_025_nitrogen_depo, nitrogen_depo),
         human_modification = ifelse(human_modification > q_975_human_modification, q_975_human_modification, human_modification),
         human_modification = ifelse(human_modification < q_025_human_modification, q_025_human_modification, human_modification),
         greenup_coef = ifelse(greenup_coef > q_975_greenup, q_975_greenup, greenup_coef),
         greenup_coef = ifelse(greenup_coef < q_025_greenup, q_025_greenup, greenup_coef),
         burned_area_coef = ifelse(burned_area_coef > q_975_burned_area, q_975_burned_area, burned_area_coef),
         burned_area_coef = ifelse(burned_area_coef < q_025_burned_area, q_025_burned_area, burned_area_coef),
         mean_evi_coef = ifelse(mean_evi_coef > q_975_evi, q_975_evi, mean_evi_coef),
         mean_evi_coef = ifelse(mean_evi_coef < q_025_evi, q_025_evi, mean_evi_coef)
  )

p_mat_shapes_sig <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey85") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)) %>% 
            filter(mat_p_value >= 0.05), color = "grey70", alpha = 1) +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)) %>% 
            filter(mat_p_value < 0.05),
          aes(color = mat_coef, fill = mat_coef), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "MAT\nTrend", fill = "MAT\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_mat_shapes_sig

ggsave(plot = p_mat_shapes_sig, "builds/plots/grid_mat_map_shapes_sig.png", dpi = 600)

p_map_shapes_sig <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey85") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)) %>% 
            filter(map_p_value >= 0.05), color = "grey70", alpha = 1) +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)) %>% 
            filter(map_p_value < 0.05),
          aes(color = map_coef, fill = map_coef), alpha = 1) +
  scale_color_scico(palette = "broc", midpoint = 0, direction = -1, end = 0.9, begin = 0.1) +
  scale_fill_scico(palette = "broc", midpoint = 0, direction = -1, end = 0.9, begin = 0.1) +
  theme_void() +
  labs(color = "MAP\nTrend", fill = "MAP\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_map_shapes_sig

ggsave(plot = p_map_shapes_sig, "builds/plots/grid_map_map_shapes_sig.png", dpi = 600)


p_max_temp_shapes_sig <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey85") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)) %>% 
            filter(max_temp_p_value >= 0.05), color = "grey70", alpha = 1) +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)) %>% 
            filter(max_temp_p_value < 0.05),
          aes(color = max_temp_coef, fill = max_temp_coef), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "Max\nTemp\nTrend", fill = "Max\nTemp\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_max_temp_shapes_sig

ggsave(plot = p_max_temp_shapes_sig, "builds/plots/grid_max_temp_map_shapes_sig.png", dpi = 600)

p_evi_shapes_sig <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey85") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)) %>% 
            filter(mean_evi_p_value >= 0.05), color = "grey70", alpha = 1) +
  geom_sf(data = dt_env %>% 
            filter(!is.na(mean_evi_coef)) %>% 
            filter(mean_evi_p_value < 0.05),
          aes(color = mean_evi_coef, fill = mean_evi_coef), alpha = 1) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  theme_void() +
  labs(color = "EVI\nTrend", fill = "EVI\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_evi_shapes_sig

ggsave(plot = p_evi_shapes_sig, "builds/plots/grid_evi_map_shapes_sig.png", dpi = 600)

p_burned_area_shapes_sig <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey85") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(burned_area_coef)) %>% 
            filter(burned_area_p_value >= 0.05), color = "grey70", alpha = 1) +
  geom_sf(data = dt_env %>% 
            filter(!is.na(burned_area_coef)) %>% 
            filter(burned_area_p_value < 0.05),
          aes(color = burned_area_coef, fill = burned_area_coef), alpha = 1) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  theme_void() +
  labs(color = "Burned\nArea\nTrend", fill = "Burned\nArea\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_burned_area_shapes_sig

ggsave(plot = p_burned_area_shapes_sig, "builds/plots/grid_burned_area_map_shapes_sig.png", dpi = 600)

p_greenup_shapes_sig <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey85") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(greenup_coef)) %>% 
            filter(greenup_p_value >= 0.05), color = "grey70", alpha = 1) +
  geom_sf(data = dt_env %>% 
            filter(!is.na(greenup_coef)) %>% 
            filter(greenup_p_value < 0.05),
          aes(color = greenup_coef, fill = greenup_coef), alpha = 1) +
  scale_color_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  theme_void() +
  labs(color = "Vegetation\nGreen-Up\nTrend", fill = "Vegetation\nGreen-Up\nTrend") +
  theme(axis.title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines"))
p_greenup_shapes_sig

ggsave(plot = p_greenup_shapes_sig, "builds/plots/grid_greenup_map_shapes_sig.png", dpi = 600)

sig_maps <- grid.arrange(p_mat_shapes_sig, p_max_temp_shapes_sig, 
                         p_map_shapes_sig, p_evi_shapes_sig, 
                         p_burned_area_shapes_sig, p_greenup_shapes_sig, ncol = 2)

ggsave(plot = sig_maps, "builds/plots/grid_sig_map_shapes.png", dpi = 600, height = 7, width = 9)

sig_maps_climate <- grid.arrange(p_mat_shapes_sig, p_max_temp_shapes_sig, 
                         p_map_shapes_sig, ncol = 2)

ggsave(plot = sig_maps_climate, "builds/plots/grid_sig_climate_map_shapes.png", dpi = 600, width = 10, height = 4.5) 

### P VALUE vs MAGNITUDE --------------

dt_p <- shapes %>% 
  left_join(climate_trends) %>% 
  mutate(abs_mean_evi_coef = abs_mean_evi_coef/100, 
         abs_burned_area_coef = abs_burned_area_coef*100, 
         abs_mat_coef = abs(mat_coef), 
         abs_map_coef = abs(map_coef), 
         abs_max_temp_coef = abs(max_temp_coef))


# Evi Trend 
evi_cor <- cor.test(dt_p$abs_mean_evi_coef, dt_p$mean_evi_p_value, method = "s")
evi_e <- round(as.numeric(evi_cor$estimate), 2)
evi_p <- round(as.numeric(evi_cor$p.value), 4)
evi_label <- paste0("Rho = ", evi_e, "; p = ", evi_p)

evi_cor_p <- dt_p %>% 
  ggplot() +
  geom_point(aes(y = mean_evi_p_value, x = abs_mean_evi_coef), alpha = 0.5) +
  annotate("text", x = 1, y = 0.75, label = evi_label) +
  labs(x = "Abs. EVI Estimate", y = "P-Value") +
  theme_classic()
evi_cor_p

# Burned Area Trend 
burned_area_cor <- cor.test(dt_p$abs_burned_area_coef, dt_p$burned_area_p_value, method = "s")
burned_area_e <- round(as.numeric(burned_area_cor$estimate), 2)
burned_area_p <- round(as.numeric(burned_area_cor$p.value), 4)
burned_area_label <- paste0("Rho = ", burned_area_e, "; p = ", burned_area_p)

burned_area_cor_p <- dt_p %>% 
  ggplot() +
  geom_point(aes(y = burned_area_p_value, x = abs_burned_area_coef), alpha = 0.5) +
  annotate("text", x = 3, y = 0.75, label = burned_area_label) +
  labs(x = "Abs. Burned Area Estimate", y = "P-Value") +
  theme_classic()
burned_area_cor_p


# Greenup Trend 
greenup_cor <- cor.test(dt_p$abs_greenup_coef, dt_p$greenup_p_value, method = "s")
greenup_e <- round(as.numeric(greenup_cor$estimate), 2)
greenup_p <- round(as.numeric(greenup_cor$p.value), 4)
greenup_label <- paste0("Rho = ", greenup_e, "; p = ", greenup_p)

greenup_cor_p <- dt_p %>% 
  ggplot() +
  geom_point(aes(y = greenup_p_value, x = abs_greenup_coef), alpha = 0.5) +
  annotate("text", x = 6, y = 0.75, label = greenup_label) +
  labs(x = "Abs. Veg-Greenup Estimate", y = "P-Value") +
  theme_classic()
greenup_cor_p


# MAT Trend 
mat_cor <- cor.test(dt_p$abs_mat_coef, dt_p$mat_p_value, method = "s")
mat_e <- round(as.numeric(mat_cor$estimate), 2)
mat_p <- round(as.numeric(mat_cor$p.value), 4)
mat_label <- paste0("Rho = ", mat_e, "; p = ", mat_p)

mat_cor_p <- dt_p %>% 
  ggplot() +
  geom_point(aes(y = mat_p_value, x = abs_mat_coef), alpha = 0.5) +
  annotate("text", x = 0.06, y = 0.75, label = mat_label) +
  labs(x = "Abs. MAT Estimate", y = "P-Value") +
  theme_classic()
mat_cor_p

#Max Temp 
max_temp_cor <- cor.test(dt_p$abs_max_temp_coef, dt_p$max_temp_p_value, method = "s")
max_temp_e <- round(as.numeric(max_temp_cor$estimate), 2)
max_temp_p <- round(as.numeric(max_temp_cor$p.value), 4)
max_temp_label <- paste0("Rho = ", max_temp_e, "; p = ", max_temp_p)

max_temp_cor_p <- dt_p %>% 
  ggplot() +
  geom_point(aes(y = max_temp_p_value, x = abs_max_temp_coef), alpha = 0.5) +
  annotate("text", x = 0.1, y = 0.75, label = max_temp_label) +
  labs(x = "Abs. Max Temp. Estimate", y = "P-Value") +
  theme_classic()
max_temp_cor_p

#MAP
map_cor <- cor.test(dt_p$abs_map_coef, dt_p$map_p_value, method = "s")
map_e <- round(as.numeric(map_cor$estimate), 2)
map_p <- round(as.numeric(map_cor$p.value), 4)
map_label <- paste0("Rho = ", map_e, "; p = ", map_p)

map_cor_p <- dt_p %>% 
  ggplot() +
  geom_point(aes(y = map_p_value, x = abs_map_coef), alpha = 0.5) +
  annotate("text", x = 40, y = 0.75, label = map_label) +
  labs(x = "Abs. MAP Estimate", y = "P-Value") +
  theme_classic()
map_cor_p

p_cor_p <- grid.arrange(evi_cor_p, burned_area_cor_p, greenup_cor_p, 
                        mat_cor_p, max_temp_cor_p, map_cor_p, ncol = 3)
ggsave(plot = p_cor_p, "builds/plots/grid_trends_vs_p.png", dpi = 600, height = 8, width = 12)
