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

evi_trend <- fread("data/processedData/dataFragments/pas_mean_evi_trends.csv") 
burned_area_trend <- fread("data/processedData/dataFragments/pas_burned_area_trends.csv")
greenup_trend <- fread("data/processedData/dataFragments/pas_greenup_trends.csv")

climate_trends <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  dplyr::select(unique_id, mat_coef, map_coef, max_temp_coef, mat_p_value, map_p_value, max_temp_p_value)

sf_use_s2(FALSE)
world <- rnaturalearth::ne_countries() %>% filter(!name_en == "Antarctica") %>%
  st_transform(crs = 'ESRI:54030') %>% 
  mutate(geometry = st_make_valid(geometry), 
         world = "world") %>% 
  group_by(world) %>% 
  summarize()


raw_shapes <- read_sf("data/spatialData/pas_and_controls.gpkg")

shapes <- raw_shapes %>%
  st_transform(crs = 'ESRI:54030') %>% 
  left_join(evi_trend) %>% 
  left_join(burned_area_trend) %>% 
  left_join(greenup_trend) %>% 
  mutate(mean_evi_coef = mean_evi_coef /100, 
         burned_area_coef = burned_area_coef*100) %>%  # to convert to %/year
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
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "strictly_protected", 
      is.na(iucn_cat) ~ "control"))


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

ggsave(plot = p_evi_shapes, "builds/plots/evi_pas_shapes_map.png", dpi = 600)

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

ggsave(plot = p_burned_area_shapes, "builds/plots/burned_area_pas_shapes_map.png", dpi = 600)

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

ggsave(plot = p_greenup_shapes, "builds/plots/greenup_pas_shapes_map.png", dpi = 600)

p_all_trends_map_shapes <- gridExtra::grid.arrange(p_evi_shapes, p_burned_area_shapes, p_greenup_shapes, ncol = 1)
ggsave(plot = p_all_trends_map_shapes, "builds/plots/pas_all_trends_map_shapes.png", dpi = 600, height = 10, width = 8)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################     ESTIMATES     ###################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

evi_est <- fread("builds/model_estimates/mean_evi_pas_estimates.csv") %>% mutate(response = "evi", 
                                                                              estimate = estimate/100, 
                                                                              std_error = std_error/100, 
                                                                              ci_lb = ci_lb/100, 
                                                                              ci_ub = ci_ub/100)
burned_area_est <- fread("builds/model_estimates/burned_area_pas_estimates.csv") %>% mutate(response = "burned_area", 
                                                                                             estimate = estimate*100, 
                                                                                             std_error = std_error*100, 
                                                                                             ci_lb = ci_lb*100, 
                                                                                             ci_ub = ci_ub*100)
greenup_est <- fread("builds/model_estimates/greenup_pas_estimates.csv") %>% mutate(response = "greenup")

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
    term == "protection_cat_broadstrictly_protected" ~ "Strictly\nProtected", 
    term == "protection_cat_broadcontrol" ~ "Control", 
    term == "area_km2_log" ~ "PA Area",
    term == "pa_age_log" ~ "PA Age"),
    sig_pn = case_when(
      estimate < 0 & p_value < 0.05 ~ "Sig. Negative", 
      estimate > 0 & p_value < 0.05 ~ "Sig. Positive", 
      p_value >= 0.05 ~ "Non-Significant"
    ),
    facet_label = case_when(
      model == "H1" ~ "H1: Trend ~\nIntercept", 
      model == "H2" ~ "H2: Trend ~\nGlobal Change", 
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

# Evi estimates -----
scico(palette = "bam", n = 10)
#"#65014B" "#9E3C85" "#C86FB1" "#E4ADD6" "#F4E3EF" "#EFF3E5" "#C0D9A1" "#7BA755" "#457B2A" "#0C4C00"

p_est_evi_b <- dt_est %>% 
  filter(response == "evi" & grepl("H", model)) %>% 
  filter(!model == "H4.1" & !model == "H5.1" & !model == "H1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 4) +
  scale_color_manual(values = c("Sig. Negative" = "#9E3C85",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "b)", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        #panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_est_evi_b

p_est_evi_c <- dt_est %>% 
  filter(response == "evi" & !grepl("H", model)) %>% 
  filter(!model == "H4.1" & !model == "H5.1" & !model == "H1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free_x", ncol = 4) +
  scale_color_manual(values = c("Sig. Negative" = "#9E3C85",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "c)", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_est_evi_c

fig_evi <- grid.arrange(p_evi_shapes + labs(title = "a)"), 
                        p_est_evi_b, p_est_evi_c, 
                        heights = c(1.3, 0.8, 1))
ggsave(plot = fig_evi, "builds/plots/pas_evi_figure.png", height = 10.5, width = 9.5, dpi = 600)

# burned area estimates -----
scico(palette = "vik", n = 10)
#"#001260" "#023E7D" "#1D6E9C" "#71A7C4" "#C9DDE7" "#EACEBE" "#D29773" "#BD6432" "#8B2706" "#590007"

p_est_burned_area_b <- dt_est %>% 
  filter(response == "burned_area" & grepl("H", model)) %>% 
  filter(!model == "H4.1" & !model == "H5.1" & !model == "H1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 4) +
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

p_est_burned_area_b

p_est_burned_area_c <- dt_est %>% 
  filter(response == "burned_area" & !grepl("H", model)) %>% 
  filter(!model == "H4.1" & !model == "H5.1" & !model == "H1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free_x", ncol = 4) +
  scale_color_manual(values = c("Sig. Negative" = "#023E7D",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#8B2706"
  )) +
  labs(title = "c)", y = NULL, x = "Burned Area Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 22.5, hjust = 1),
        plot.title = element_text(size = 12))

p_est_burned_area_c

fig_burned_area <- grid.arrange(p_burned_area_shapes + labs(title = "a)"), 
                                p_est_burned_area_b, p_est_burned_area_c, 
                                heights = c(1.3, 0.8, 1))
ggsave(plot = fig_burned_area, "builds/plots/pas_burned_area_figure.png", height = 10.5, width = 9.5, dpi = 600)

# Greenup estimates -----
scico(palette = "cork", n = 10, direction = -1, begin = 0.1, end = 0.9)
#"#195715" "#3D7D3C" "#6C9C6B" "#A2C0A1" "#D8E5D9" "#D2DDE7" "#97AFC8" "#6388AD" "#386695" "#284074"

p_est_greenup_b <- dt_est %>% 
  filter(response == "greenup" & grepl("H", model)) %>% 
  filter(!model == "H4.1" & !model == "H5.1" & !model == "H1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 4) +
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
        plot.title = element_text(size = 12))

p_est_greenup_b

p_est_greenup_c <- dt_est %>% 
  filter(response == "greenup" & !grepl("H", model)) %>% 
  filter(!model == "H4.1" & !model == "H5.1" & !model == "H1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free_x", ncol = 4) +
  scale_color_manual(values = c("Sig. Negative" = "#3D7D3C",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#386695"
  )) +
  labs(title = "c)", y = NULL, x = "Vegetation Green-Up Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_est_greenup_c

fig_greenup <- grid.arrange(p_greenup_shapes + labs(title = "a)"), 
                            p_est_greenup_b, p_est_greenup_c, 
                            heights = c(1.3, 0.8, 1))
ggsave(plot = fig_greenup, "builds/plots/pas_greenup_figure.png", height = 10.5, width = 9.5, dpi = 600)


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
    nitrogen_depo, human_modification) %>% 
  filter(complete.cases(.)) %>% 
  rename(`EVI Trend` = mean_evi_coef, 
         `Burned Area Trend` = burned_area_coef, 
         `Greenup Trend` = greenup_coef, 
         `MAT Trend` = mat_coef,
         `MAP Trend` = map_coef, 
         `Max Temp Trend` = max_temp_coef, 
         `Nitrogen Deposition` = nitrogen_depo, 
         `Human Modification` = human_modification)

library(ggcorrplot)
corr <- round(cor(dt_corr), 1)
p_corr <- ggcorrplot(corr, hc.order = TRUE, type = "lower",
                     lab = TRUE)
p_corr
ggsave(plot = p_corr, "builds/plots/pas_variable_correlations.png", dpi = 600, height = 8, width = 8)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
############################     TREND DISTRIBUTION     ##############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(ggridges)

dt_ridges <- shapes %>% 
  as.data.table() %>%
  mutate(geom = NULL, 
         mean_evi_coef = ifelse(mean_evi_coef > q_975_evi, NA, mean_evi_coef),
         mean_evi_coef = ifelse(mean_evi_coef < q_025_evi, NA, mean_evi_coef), 
         burned_area_coef = ifelse(burned_area_coef > q_975_burned_area, NA, burned_area_coef),
         burned_area_coef = ifelse(burned_area_coef < q_025_burned_area, NA, burned_area_coef), 
         greenup_coef = ifelse(greenup_coef > q_975_greenup, NA, greenup_coef),
         greenup_coef = ifelse(greenup_coef < q_025_greenup, NA, greenup_coef), 
         protection_status = case_when(
           protection_cat_broad == "strictly_protected" ~ "Strictly\nProtected",
           protection_cat_broad == "control" ~ "Control"),
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
        axis.text.y = element_text(size = 12))

p_greenup_biome

p_greenup_dens <- gridExtra::grid.arrange(p_greenup_prot, p_greenup_biome, ncol = 1)

p_ridges <- gridExtra::grid.arrange(p_evi_dens, p_burned_area_dens, p_greenup_dens, ncol = 3)
ggsave(plot = p_ridges, "builds/plots/pas_ridges.png", dpi = 600, height = 6, width = 12)


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

ggsave(plot = p_biome_shapes, "builds/plots/pas_biome_maps_shapes.png", dpi = 600)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   PROTECTION MAP   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
scico(palette = "bamako", n = 10)
#[1] "#003A46" "#0E433F" "#1F4E34" "#355E26" "#527014" "#728202" "#988C02" "#BEA82E" "#E1C76D" "#FFE5AC"

p_pa_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = shapes %>% 
            mutate(protection_cat_broad = ifelse(protection_cat_broad == "control", "Control", "Protected")) %>% 
            filter(!is.na(mean_evi_coef)),
          aes(color = protection_cat_broad, fill = protection_cat_broad), alpha = 1) +
  scale_color_scico_d(palette = "bamako", begin = 0.25, end = 0.75) +
  scale_fill_scico_d(palette = "bamako", begin = 0.25, end = 0.75) +
  theme_void() +
  labs(color = "Protection\nStatus", fill = "Protection\nStatus") +
  theme(axis.title = element_blank())
p_pa_shapes

ggsave(plot = p_pa_shapes, "builds/plots/pas_protection_map_shapes.png", dpi = 600)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################ STOP HERE ##########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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

ggsave(plot = p_mat_shapes, "builds/plots/pas_mat_map_shapes.png", dpi = 600)

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

ggsave(plot = p_map_shapes, "builds/plots/pas_map_map_shapes.png", dpi = 600)


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

ggsave(plot = p_max_temp_shapes, "builds/plots/pas_max_temp_map_shapes.png", dpi = 600)


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

ggsave(plot = p_n_depo_shapes, "builds/plots/pas_n_depo_map_shapes.png", dpi = 600)

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

ggsave(plot = p_human_modification_shapes, "builds/plots/pas_human_modification_map_shapes.png", dpi = 600)

gc_maps <- grid.arrange(p_n_depo_shapes, p_human_modification_shapes, 
                        p_mat_shapes, p_max_temp_shapes, 
                        p_map_shapes, ncol = 2)

ggsave(plot = gc_maps, "builds/plots/pas_global_change_map_shapes.png", dpi = 600)

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

ggsave(plot = p_mat_shapes_sig, "builds/plots/pas_mat_map_shapes_sig.png", dpi = 600)

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

ggsave(plot = p_map_shapes_sig, "builds/plots/pas_map_map_shapes_sig.png", dpi = 600)


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

ggsave(plot = p_max_temp_shapes_sig, "builds/plots/pas_max_temp_map_shapes_sig.png", dpi = 600)

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

ggsave(plot = p_evi_shapes_sig, "builds/plots/pas_evi_map_shapes_sig.png", dpi = 600)

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

ggsave(plot = p_burned_area_shapes_sig, "builds/plots/pas_burned_area_map_shapes_sig.png", dpi = 600)

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

ggsave(plot = p_greenup_shapes_sig, "builds/plots/pas_greenup_map_shapes_sig.png", dpi = 600)

sig_maps <- grid.arrange(p_mat_shapes_sig, p_max_temp_shapes_sig, 
                         p_map_shapes_sig, p_evi_shapes_sig, 
                         p_burned_area_shapes_sig, p_greenup_shapes_sig, ncol = 2)

ggsave(plot = sig_maps, "builds/plots/pas_sig_map_shapes.png", dpi = 600, height = 7, width = 9)

sig_maps_climate <- grid.arrange(p_mat_shapes_sig, p_max_temp_shapes_sig, 
                                 p_map_shapes_sig, ncol = 2)

ggsave(plot = sig_maps_climate, "builds/plots/pas_sig_climate_map_shapes.png", dpi = 600, width = 10, height = 4.5) 

