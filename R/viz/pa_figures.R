### make nice maps ####

library(tidyverse)
library(data.table)
library(RColorBrewer)
library(scico)
library(rnaturalearth)
library(sf)
library(gridExtra)

evi_trend <- fread("data/processedData/dataFragments/pa_evi_trends.csv")
burned_area_trend <- fread("data/processedData/dataFragments/pa_burned_area_trends.csv")
greenup_trend <- fread("data/processedData/dataFragments/pa_greenup_trends.csv")


climate_trends <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  dplyr::select(unique_id, mat_coef, map_coef, max_temp_coef)

sf_use_s2(FALSE)
world <- rnaturalearth::ne_countries() %>% filter(!name_en == "Antarctica") %>%
  st_transform(crs = 'ESRI:54030') %>% 
  mutate(geometry = st_make_valid(geometry)) %>% 
  group_by(continent) %>% 
  summarize()


raw_shapes <- read_sf("data/spatialData/pas_and_controls.gpkg")

shapes <- raw_shapes %>%
  st_transform(crs = 'ESRI:54030') %>% 
  left_join(evi_trend) %>% 
  left_join(burned_area_trend) %>% 
  left_join(greenup_trend) %>% 
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
    ))

coords <- st_coordinates(st_centroid(shapes))
shapes$X <- coords[,1]
shapes$Y <- coords[,2]

sum(is.na(shapes$evi_coef))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
##################################     MAPS    #######################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

## evi maps -------------
quantile(shapes$evi_coef, c(.025, .975), na.rm = T)

q_025_evi <- as.numeric(quantile(shapes$evi_coef, c(.025), na.rm = T))
q_975_evi <- as.numeric(quantile(shapes$evi_coef, c(.975), na.rm = T))

p_evi_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = shapes %>%
            filter(!is.na(evi_coef)) %>% 
            mutate(evi_coef = ifelse(evi_coef > q_975_evi, q_975_evi, evi_coef),
                   evi_coef = ifelse(evi_coef < q_025_evi, q_025_evi, evi_coef)),
          aes(color = evi_coef, fill = evi_coef)) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  labs(color = "EVI\nTrend", fill = "EVI\nTrend") +
  theme_void() +
  theme(legend.position = "right", 
       # legend.key.width = unit(1, "cm"),
       # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 25))
p_evi_shapes

ggsave(plot = p_evi_shapes, "builds/plots/evi_pa_shapes_map.png", dpi = 600)

p_evi_points <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_point(data = shapes %>%
            filter(!is.na(evi_coef)) %>% 
            mutate(evi_coef = ifelse(evi_coef > q_975_evi, q_975_evi, evi_coef),
                   evi_coef = ifelse(evi_coef < q_025_evi, q_025_evi, evi_coef)),
          aes(x = X, y = Y, color = evi_coef, fill = evi_coef)) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  labs(color = "EVI\nTrend", fill = "EVI\nTrend") +
  theme_void() +
  theme(legend.position = "right", 
        # legend.key.width = unit(1, "cm"),
        # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 25))
p_evi_points

ggsave(plot = p_evi_points, "builds/plots/evi_pa_points_map.png", dpi = 600)

## burned area maps-------
quantile(shapes$burned_area_coef, c(.025, .975), na.rm = T)

q_025_burned_area <- as.numeric(quantile(shapes$burned_area_coef, c(.025), na.rm = T))
q_975_burned_area <- as.numeric(quantile(shapes$burned_area_coef, c(.975), na.rm = T))

p_burned_area_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
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
        legend.text = element_text(angle = 25))
p_burned_area_shapes

ggsave(plot = p_burned_area_shapes, "builds/plots/burned_area_pa_shapes_map.png", dpi = 600)

p_burned_area_points <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_point(data = shapes %>%
               filter(!is.na(burned_area_coef)) %>% 
               mutate(burned_area_coef = ifelse(burned_area_coef > q_975_burned_area, q_975_burned_area, burned_area_coef),
                      burned_area_coef = ifelse(burned_area_coef < q_025_burned_area, q_025_burned_area, burned_area_coef)),
             aes(x = X, y = Y, color = burned_area_coef, fill = burned_area_coef)) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  labs(color = "Burned\nArea\nTrend", fill = "Burned Area\nTrend") +
  theme_void() +
  theme(legend.position = "right", 
        # legend.key.width = unit(1, "cm"),
        # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 25))
p_burned_area_points

ggsave(plot = p_burned_area_points, "builds/plots/burned_area_pa_points_map.png", dpi = 600)

## greenup maps ---------------------
quantile(shapes$greenup_coef, c(.025, .975), na.rm = T)

q_025_greenup <- as.numeric(quantile(shapes$greenup_coef, c(.025), na.rm = T))
q_975_greenup <- as.numeric(quantile(shapes$greenup_coef, c(.975), na.rm = T))

p_greenup_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = shapes %>%
            filter(!is.na(greenup_coef)) %>% 
            mutate(greenup_coef = ifelse(greenup_coef > q_975_greenup, q_975_greenup, greenup_coef),
                   greenup_coef = ifelse(greenup_coef < q_025_greenup, q_025_greenup, greenup_coef)),
          aes(color = greenup_coef, fill = greenup_coef)) +
  scale_color_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  labs(color = "Greenup\nTrend", fill = "Greenup\nTrend") +
  theme_void() +
  theme(legend.position = "right", 
        # legend.key.width = unit(1, "cm"),
        # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 25))
p_greenup_shapes

ggsave(plot = p_greenup_shapes, "builds/plots/greenup_pa_shapes_map.png", dpi = 600)

p_greenup_points <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_point(data = shapes %>%
               filter(!is.na(greenup_coef)) %>% 
               mutate(greenup_coef = ifelse(greenup_coef > q_975_greenup, q_975_greenup, greenup_coef),
                      greenup_coef = ifelse(greenup_coef < q_025_greenup, q_025_greenup, greenup_coef)),
             aes(x = X, y = Y, color = greenup_coef, fill = greenup_coef)) +
  scale_color_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_scico(palette = "cork", midpoint = 0, direction = -1, begin = 0.1, end = 0.9) +
  labs(color = "Greenup\nTrend", fill = "Greenup\nTrend") +
  theme_void() +
  theme(legend.position = "right", 
        # legend.key.width = unit(1, "cm"),
        # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 25))
p_greenup_points

ggsave(plot = p_greenup_points, "builds/plots/greenup_pa_points_map.png", dpi = 600)

p_all_trends_map_points <- gridExtra::grid.arrange(p_evi_points, p_burned_area_points, p_greenup_points, ncol = 1)
ggsave(plot = p_all_trends_map_points, "builds/plots/pa_all_trends_map_points.png", dpi = 600, height = 10, width = 8)

p_all_trends_map_shapes <- gridExtra::grid.arrange(p_evi_shapes, p_burned_area_shapes, p_greenup_shapes, ncol = 1)
ggsave(plot = p_all_trends_map_shapes, "builds/plots/pa_all_trends_map_shapes.png", dpi = 600, height = 10, width = 8)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################     ESTIMATES     ###################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

evi_est <- fread("builds/model_estimates/evi_pa_estimates.csv") %>% mutate(response = "evi")
burned_area_est <- fread("builds/model_estimates/burned_area_pa_estimates.csv") %>% mutate(response = "burned_area")
greenup_est <- fread("builds/model_estimates/greenup_pa_estimates.csv") %>% mutate(response = "greenup")

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
    term == "og_layercontrols" ~ "Control", 
    term == "og_layerprotected_areas" ~ "Protected", 
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
                                "Non-Significant" = "grey", 
                                "Sig. Positive" = "#457B2A"
                                )) +
  labs(title = "b)", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
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
                                "Non-Significant" = "grey", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "c)", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        plot.title = element_text(size = 12))

p_est_evi_c

fig_evi <- grid.arrange(p_evi_shapes + labs(title = "a)"), 
                        p_est_evi_b, p_est_evi_c, 
                        heights = c(1.3, 0.8, 1))
ggsave(plot = fig_evi, "builds/plots/pa_evi_figure.png", height = 10.5, width = 9, dpi = 600)

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
                                "Non-Significant" = "grey", 
                                "Sig. Positive" = "#8B2706"
  )) +
  labs(title = "b)", y = NULL, x = "Burned Area Trend Estimate", color = "Significance") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
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
                                "Non-Significant" = "grey", 
                                "Sig. Positive" = "#8B2706"
  )) +
  labs(title = "c)", y = NULL, x = "Burned Area Trend Estimate", color = "Significance") +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        plot.title = element_text(size = 12))

p_est_burned_area_c

fig_burned_area <- grid.arrange(p_burned_area_shapes + labs(title = "a)"), 
                        p_est_burned_area_b, p_est_burned_area_c, 
                        heights = c(1.3, 0.8, 1))
ggsave(plot = fig_burned_area, "builds/plots/pa_burned_area_figure.png", height = 10.5, width = 9, dpi = 600)

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
                                "Non-Significant" = "grey", 
                                "Sig. Positive" = "#386695"
  )) +
  labs(title = "b)", y = NULL, x = "Burned Area Trend Estimate", color = "Significance") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
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
                                "Non-Significant" = "grey", 
                                "Sig. Positive" = "#386695"
  )) +
  labs(title = "c)", y = NULL, x = "Burned Area Trend Estimate", color = "Significance") +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.grid = element_blank(),
        plot.title = element_text(size = 12))

p_est_greenup_c

fig_greenup <- grid.arrange(p_greenup_shapes + labs(title = "a)"), 
                                p_est_greenup_b, p_est_greenup_c, 
                                heights = c(1.3, 0.8, 1))
ggsave(plot = fig_greenup, "builds/plots/pa_greenup_figure.png", height = 10.5, width = 9, dpi = 600)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################     CORRELATION     ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

dt_corr <- shapes %>% 
  left_join(climate_trends) %>% 
  rename(pa_age = PaAge, 
         pa_area = area_km2, 
         nitrogen_depo = NitrogenDepo, 
         human_modification = HumanModification
         ) %>% 
  as.data.table() %>% 
  dplyr::select(
                evi_coef, burned_area_coef, greenup_coef,
                mat_coef, map_coef, max_temp_coef, 
                nitrogen_depo, human_modification) %>% 
  filter(complete.cases(.)) %>% 
  rename(`EVI Trend` = evi_coef, 
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
ggsave(plot = p_corr, "builds/plots/pa_variable_correlations.png", dpi = 600, height = 8, width = 8)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
############################     TREND DISTRIBUTION     ##############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
 
library(ggridges)

dt_ridges <- shapes %>% 
  as.data.table() %>%
  mutate(geom = NULL, 
         evi_coef = ifelse(evi_coef > q_975_evi, NA, evi_coef),
         evi_coef = ifelse(evi_coef < q_025_evi, NA, evi_coef), 
         burned_area_coef = ifelse(burned_area_coef > q_975_burned_area, NA, burned_area_coef),
         burned_area_coef = ifelse(burned_area_coef < q_025_burned_area, NA, burned_area_coef), 
         greenup_coef = ifelse(greenup_coef > q_975_greenup, NA, greenup_coef),
         greenup_coef = ifelse(greenup_coef < q_025_greenup, NA, greenup_coef), 
         protection_status = case_when(
           og_layer == "protected_areas" ~ "Protected",
           og_layer == "controls" ~ "Control"), 
         biome_clean = case_when(
           super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
           super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
           super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
           super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation"))

#evi ridges
p_evi_prot <- ggplot() +
  geom_density_ridges_gradient(data = dt_ridges %>%
                                 filter(!is.na(evi_coef)),
                               aes(x = evi_coef, y = protection_status, fill = ..x..), alpha = 0.7) +
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
                                 filter(!is.na(evi_coef)),
                               aes(x = evi_coef, y = biome_clean, fill = ..x..), alpha = 0.7) +
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
  labs(y = "", x = "Greenup Trend") + 
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
  labs(y = "", x = "Greenup Trend") + 
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 12))

p_greenup_biome

p_greenup_dens <- gridExtra::grid.arrange(p_greenup_prot, p_greenup_biome, ncol = 1)

p_ridges <- gridExtra::grid.arrange(p_evi_dens, p_burned_area_dens, p_greenup_dens, ncol = 3)
ggsave(plot = p_ridges, "builds/plots/pa_ridges.png", dpi = 600, height = 6, width = 12)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################     BIOME MAP     ####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

p_biome_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = shapes %>% 
            filter(!is.na(evi_coef)) %>%
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

ggsave(plot = p_biome_shapes, "builds/plots/pa_biome_maps_shapes.png", dpi = 600)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   PROTECTION MAP   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

p_prot_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = shapes %>% 
               filter(!is.na(evi_coef)) %>%
               mutate(prot_clean = case_when(
                 og_layer == "protected_areas" ~ "Protected",
                 og_layer == "controls" ~ "Control")),
             aes(color = prot_clean, fill = prot_clean), alpha = 1) +
  scale_color_scico_d(palette = "batlow") +
  scale_fill_scico_d(palette = "batlow") +
  theme_void() +
  labs(color = "Protection\nStatus", fill = "Protection\nStatus") +
  theme(axis.title = element_blank())
p_prot_shapes


ggsave(plot = p_prot_shapes, "builds/plots/pa_protection_map_points.png", dpi = 600)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
############################     GLOBAL CHANGE MAPS     ##############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
dt_env_raw <- shapes %>% 
  left_join(climate_trends)  %>% 
  rename(nitrogen_deposition = NitrogenDepo, 
         human_modification = HumanModification)

q_025_mat <- as.numeric(quantile(dt_env_raw$mat_coef, c(.025), na.rm = T))
q_975_mat <- as.numeric(quantile(dt_env_raw$mat_coef, c(.975), na.rm = T))
q_025_map <- as.numeric(quantile(dt_env_raw$map_coef, c(.025), na.rm = T))
q_975_map <- as.numeric(quantile(dt_env_raw$map_coef, c(.975), na.rm = T))
q_025_max_temp <- as.numeric(quantile(dt_env_raw$max_temp_coef, c(.025), na.rm = T))
q_975_max_temp <- as.numeric(quantile(dt_env_raw$max_temp_coef, c(.975), na.rm = T))
q_025_nitrogen_deposition <- as.numeric(quantile(dt_env_raw$nitrogen_deposition, c(.025), na.rm = T))
q_975_nitrogen_deposition <- as.numeric(quantile(dt_env_raw$nitrogen_deposition, c(.975), na.rm = T))
q_025_human_modification <- as.numeric(quantile(dt_env_raw$human_modification, c(.025), na.rm = T))
q_975_human_modification <- as.numeric(quantile(dt_env_raw$human_modification, c(.975), na.rm = T))

dt_env <-  dt_env_raw %>%
  mutate(mat_coef = ifelse(mat_coef > q_975_mat, q_975_mat, mat_coef),
         mat_coef = ifelse(mat_coef < q_025_mat, q_025_mat, mat_coef),
         map_coef = ifelse(map_coef > q_975_map, q_975_map, map_coef),
         map_coef = ifelse(map_coef < q_025_map, q_025_map, map_coef),
         max_temp_coef = ifelse(max_temp_coef > q_975_max_temp, q_975_max_temp, max_temp_coef),
         max_temp_coef = ifelse(max_temp_coef < q_025_max_temp, q_025_max_temp, max_temp_coef),
         nitrogen_deposition = ifelse(nitrogen_deposition > q_975_nitrogen_deposition, q_975_nitrogen_deposition, nitrogen_deposition),
         nitrogen_deposition = ifelse(nitrogen_deposition < q_025_nitrogen_deposition, q_025_nitrogen_deposition, nitrogen_deposition),
         human_modification = ifelse(human_modification > q_975_human_modification, q_975_human_modification, human_modification),
         human_modification = ifelse(human_modification < q_025_human_modification, q_025_human_modification, human_modification)
  )

p_mat_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(evi_coef)),
          aes(color = mat_coef, fill = mat_coef), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "MAT\nTrend", fill = "MAT\nTrend") +
  theme(axis.title = element_blank())
p_mat_shapes

ggsave(plot = p_mat_shapes, "builds/plots/pa_mat_map_shapes.png", dpi = 600)

p_map_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(evi_coef)),
          aes(color = map_coef, fill = map_coef), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "MAP\nTrend", fill = "MAP\nTrend") +
  theme(axis.title = element_blank())
p_map_shapes

ggsave(plot = p_map_shapes, "builds/plots/pa_map_map_shapes.png", dpi = 600)


p_max_temp_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(evi_coef)),
          aes(color = max_temp_coef, fill = max_temp_coef), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "Max\nTemp\nTrend", fill = "Max\nTemp\nTrend") +
  theme(axis.title = element_blank())
p_max_temp_shapes

ggsave(plot = p_max_temp_shapes, "builds/plots/pa_max_temp_map_shapes.png", dpi = 600)

p_n_depo_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(evi_coef)),
          aes(color = nitrogen_deposition, fill = nitrogen_deposition), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "Nitrogen\nDeposition", fill = "Nitrogen\nDeposition") +
  theme(axis.title = element_blank())
p_n_depo_shapes

ggsave(plot = p_n_depo_shapes, "builds/plots/pa_n_depo_map_shapes.png", dpi = 600)

p_human_modification_shapes <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey75") +
  geom_sf(data = dt_env %>% 
            filter(!is.na(evi_coef)),
          aes(color = human_modification, fill = human_modification), alpha = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  labs(color = "Human\nModification", fill = "Human\nModification") +
  theme(axis.title = element_blank())
p_human_modification_shapes

ggsave(plot = p_human_modification_shapes, "builds/plots/pa_human_modification_map_shapes.png", dpi = 600)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#################################     QUANTILES     ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

dt_env_raw <- shapes %>% 
  left_join(climate_trends)  %>% 
  rename(nitrogen_deposition = NitrogenDepo, 
         human_modification = HumanModification)



quantile(dt_env_raw[!is.na(dt_env_raw$evi_coef),]$mat_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#-0.040871137  0.002113063  0.011041658  0.017943548  0.025513459  0.030837586  0.061575663 
quantile(dt_env_raw[!is.na(dt_env_raw$evi_coef),]$max_temp_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#-0.066161953 -0.005806273  0.007165408  0.015833825  0.024810994  0.036342111  0.139455726 
quantile(dt_env_raw[!is.na(dt_env_raw$evi_coef),]$map_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100% 
#-29.67555174  -2.07683237   0.09181447   0.79981597   1.60869519   4.49253768  15.60217714 
quantile(dt_env_raw[!is.na(dt_env_raw$evi_coef),]$nitrogen_deposition, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#        0%         5%        25%        50%        75%        95%       100% 
#  19.41463  100.25536  215.19000  332.70844  579.09265 1179.13915 4738.10938 
quantile(dt_env_raw[!is.na(dt_env_raw$evi_coef),]$human_modification, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#         0%          5%         25%         50%         75%         95%        100%
#0.000000000 0.002115423 0.039965609 0.113280658 0.230818849 0.449265328 0.885243416 
quantile(dt_env_raw[!is.na(dt_env_raw$evi_coef),]$evi_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#         0%          5%         25%         50%         75%         95%        100%
#-98.4236617 -13.7772218   0.1529625   7.2733991  15.5606324  31.1578728 113.6631148 
quantile(dt_env_raw[!is.na(dt_env_raw$evi_coef),]$burned_area_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#          0%           5%          25%          50%          75%          95%         100%  
#-0.038989831 -0.001550148  0.000000000  0.000000000  0.000000000  0.001209549  0.044689271 
quantile(dt_env_raw[!is.na(dt_env_raw$evi_coef),]$greenup_coef, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)
#         0%          5%         25%         50%         75%         95%        100%
#-7.92594632 -1.11387837 -0.22116058 -0.03874859  0.13539879  1.30038490  6.70514126 