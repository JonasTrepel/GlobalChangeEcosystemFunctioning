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

climate_trends <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") 


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
  filter(unique_id %in% c(unique(dt_sum$unique_id))) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "strictly_protected", 
      is.na(iucn_cat) ~ "control")) %>% 
  left_join(dt_sum)


sum(is.na(shapes$mean_evi_coef))

sum(shapes[shapes$protection_cat_broad == "strictly_protected", ]$area_km2, na.rm = T)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   PROTECTION MAP   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
scico(palette = "bamako", n = 10)
#[1] "#003A46" "#0E433F" "#1F4E34" "#355E26" "#527014" "#728202" "#988C02" "#BEA82E" "#E1C76D" "#FFE5AC"

p_pa_shapes <- ggplot() +
  geom_sf(data = world, fill = "white", color = "grey75") +
  geom_sf(data = shapes %>% 
            mutate(protection_cat_broad = ifelse(protection_cat_broad == "control", "Control", "Protected")),
          aes(color = protection_cat_broad, fill = protection_cat_broad), alpha = 1) +
  scale_color_scico_d(palette = "bamako", begin = 0.25, end = 0.75) +
  scale_fill_scico_d(palette = "bamako", begin = 0.25, end = 0.75) +
  theme_void() +
  labs(color = "Protection\nStatus", fill = "Protection\nStatus") +
  theme(axis.title = element_blank(), 
        legend.position = "bottom")
p_pa_shapes

ggsave(plot = p_pa_shapes, "builds/plots/pas_protection_map_shapes.png", dpi = 600)

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
      model == "H4" ~ "Trend ~\nProtection",  #PA mean
      model == "H4.1" ~ "Abs. Trend ~\nProtection", #PA mean 
      model == "H5" ~ "Trend ~\nPA Characteristics",  #PA mean
      model == "H5.1" ~ "Abs. Trend ~\nPA Characteristics",  #PA mean
      model == "H6" ~ "Trend ~\nPA Characteristics", #PA median
      model == "H6.1" ~ "Abs. Trend ~\nPA Characteristics", #PA median 
      model == "H7" ~ "Trend ~\nProtection", #PA median 
      model == "H7.1" ~ "Abs. Trend ~\nProtection" #PA median
    )
  )

#H4 - mean protection 
#H5- mean prot characteristics  
#H6: median protection 
#H7: median prot characteristics 

### Pa protection estimtes --------
# evi 
p_pa_est_evi <- dt_est %>% 
  filter(response == "evi" & grepl("H", model)) %>% 
  filter(model == "H4" | model == "H6") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 12) +
  scale_color_manual(values = c("Sig. Negative" = "#9E3C85",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "a)", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow", color = "snow"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))
p_pa_est_evi

p_pa_est_abs_evi <- dt_est %>% 
  filter(response == "evi" & grepl("H", model)) %>% 
  filter(model == "H4.1" | model == "H6.1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 12) +
  scale_color_manual(values = c("Sig. Negative" = "#9E3C85",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "b)", y = NULL, x = "EVI Trend Estimate (absolute)", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow", color = "snow"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_pa_est_abs_evi

# burned area
p_pa_est_burned_area <- dt_est %>% 
  filter(response == "burned_area" & grepl("H", model)) %>% 
  filter(model == "H4" | model == "H6") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 12) +
  scale_color_manual(values = c("Sig. Negative" = "#023E7D",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#8B2706"
  )) +
  labs(title = "c)", y = NULL, x = "Burned Area Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow", color = "snow"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12, color = "black"))
p_pa_est_burned_area

p_pa_est_abs_burned_area <- dt_est %>% 
  filter(response == "burned_area" & grepl("H", model)) %>% 
  filter(model == "H4.1" | model == "H6.1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 12) +
  scale_color_manual(values = c("Sig. Negative" = "#023E7D",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#8B2706"
  )) +
  labs(title = "d)", y = NULL, x = "Burned Area Trend Estimate (absolute)", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow", color = "snow"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12, color = "black"))

p_pa_est_abs_burned_area


# greenup
p_pa_est_greenup <- dt_est %>% 
  filter(response == "greenup" & grepl("H", model)) %>% 
  filter(model == "H4" | model == "H6") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 12) +
  scale_color_manual(values = c("Sig. Negative" = "#3D7D3C",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#386695"
  )) +
  labs(title = "e)", y = NULL, x = "Vegetation Green-Up Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow", color = "snow"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12, color = "black"))
p_pa_est_greenup

p_pa_est_abs_greenup <- dt_est %>% 
  filter(response == "greenup" & grepl("H", model)) %>% 
  filter(model == "H4.1" | model == "H5.1") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free", ncol = 12) +
  scale_color_manual(values = c("Sig. Negative" = "#3D7D3C",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#386695"
  )) +
  labs(title = "f)", y = NULL, x = "Vegetation Green-Up Trend Estimate (absolute)", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow", color = "snow"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12, color = "black"))

p_pa_est_abs_greenup

empty_plot <- ggplot() +
  theme_void()

p_pa_a <- gridExtra::grid.arrange(p_pa_est_evi, p_pa_est_burned_area, p_pa_est_greenup)
p_pa_b <- gridExtra::grid.arrange(p_pa_est_abs_evi, p_pa_est_abs_burned_area, p_pa_est_abs_greenup)

p_pa_est <- gridExtra::grid.arrange(p_pa_a, empty_plot, p_pa_b, ncol = 3, widths = c(1, 0.1, 1)) 
ggsave(plot = p_pa_est, "builds/plots/pa_protection_estimates.png", dpi = 600, height = 6, width = 9)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
############################     TREND DISTRIBUTION     ##############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(ggridges)

dt_ridges <- dt_sum %>% 
  as.data.table() %>%
  mutate(
         mean_evi_coef = ifelse(mean_evi_coef > q_975_evi, NA, mean_evi_coef),
         mean_evi_coef = ifelse(mean_evi_coef < q_025_evi, NA, mean_evi_coef), 
         burned_area_coef = ifelse(burned_area_coef > q_975_burned_area, NA, burned_area_coef),
         burned_area_coef = ifelse(burned_area_coef < q_025_burned_area, NA, burned_area_coef), 
         greenup_coef = ifelse(greenup_coef > q_975_greenup, NA, greenup_coef),
         greenup_coef = ifelse(greenup_coef < q_025_greenup, NA, greenup_coef))


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

p_evi_dens <- gridExtra::grid.arrange(p_evi_prot, p_evi_biome, ncol = 12)


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

p_burned_area_dens <- gridExtra::grid.arrange(p_burned_area_prot, p_burned_area_biome, ncol = 12)

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

p_greenup_dens <- gridExtra::grid.arrange(p_greenup_prot, p_greenup_biome, ncol = 12)

p_ridges <- gridExtra::grid.arrange(p_evi_dens, p_burned_area_dens, p_greenup_dens, ncol = 3)
ggsave(plot = p_ridges, "builds/plots/pas_ridges.png", dpi = 600, height = 6, width = 12)


#### area and size distribution 
p_pa_age <- ggplot() +
  geom_density_ridges(data = dt_ridges %>%
                        filter(!pa_age == 0), 
                      aes(x = pa_age, y = biome_clean, fill = biome_clean, color = biome_clean), alpha = 0.7) +
  scale_color_scico_d("batlowK") +
  scale_fill_scico_d("batlowK") +
  scale_x_log10() +
  theme_bw() +
  labs(y = "", x = "PA Age (Years)", title = "a)") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_pa_age

p_pa_area <- ggplot() +
  geom_density_ridges(data = dt_ridges,
                      aes(x = area_km2, y = biome_clean, fill = biome_clean, color = biome_clean), alpha = 0.7) +
  scale_color_scico_d("batlowK") +
  scale_fill_scico_d("batlowK") +
  scale_x_log10(labels = scales::label_comma()) +
  theme_bw() +
  labs(y = "", x = "PA Area (km^2)", title = "b)") + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12))

p_pa_area

p_pa_area_trend <- grid.arrange(p_pa_age, p_pa_area, ncol = 2)
ggsave(plot = p_pa_area_trend, "builds/plots/pas_pa_area_and_age_dist.png", dpi = 600, height = 4, width = 8)


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

ggsave(plot = p_corr_pa, "builds/plots/pas_variable_correlations_strict_pas_all.png", dpi = 600, height = 13, width = 11)





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

