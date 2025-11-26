library(tidyverse)
library(data.table)
library(sf)
library(ggridges)
library(rnaturalearth)
library(scico)
library(effectsize)
library(patchwork)

#### Load all data ------------

dt <- fread("data/processed_data/analysis_ready_pa_dataset.csv")
dt_pred <-  fread("builds/model_outputs/pa_dataset_global_change_predictions.csv")
dt_bm = fread("builds/model_outputs/sdmtmb_results_pa_dataset.csv")

sf_strict = st_read("data/spatial_data/protected_areas/pa_shapes_no_overlap.gpkg") %>% 
  mutate(protection_cat_broad = "protected") %>% 
  dplyr::select(unique_id, protection_cat_broad)
sf_control = st_read("data/spatial_data/protected_areas/controls_for_pas.gpkg") %>% 
  mutate(protection_cat_broad = "control") %>% 
  dplyr::select(unique_id, protection_cat_broad)

sf_pa = rbind(sf_strict, sf_control) %>% 
  filter(unique_id %in% unique(dt$unique_id))

#### Map 

world <- rnaturalearth::ne_countries() %>% filter(!name_en == "Antarctica") %>%
  st_transform(crs = 'ESRI:54009') %>% 
  mutate(geometry = st_make_valid(geometry), 
         world = "world") %>% 
  group_by(world) %>% 
  summarize()



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   PROTECTION MAP   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
scico(palette = "managua", n = 10)
c("#FFCE66","#DC9954","#B96C46","#92463A", "#662A3C","#4E305D","#4D5492","#5B80BC","#6CB0DD","#80E6FF")


world_crop <- st_crop(world, st_bbox(sf_pa))

plot(world_crop)

p_pa_shapes <- ggplot() +
  geom_sf(data = world_crop, fill = "white", color = "grey75") +
  geom_sf(data = sf_pa,
                    aes(color = protection_cat_broad, fill = protection_cat_broad), alpha = 1) +
  scale_fill_manual(values = c("#DC9954", "#5B80BC")) +
  scale_color_manual(values = c("#DC9954", "#5B80BC")) +
  theme_void() +
  labs(color = "Protection\nStatus", fill = "Protection\nStatus") +
  theme(axis.title = element_blank(), 
        legend.position = "bottom")
p_pa_shapes

ggsave(plot = p_pa_shapes, "builds/plots/pas_protection_map_shapes.png", dpi = 900)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   Trend comparison   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

p_pred <- dt_pred%>% 
  filter(x_unscaled > q05_unscaled,  x_unscaled < q95_unscaled) %>% 
  ggplot() +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_ribbon(aes(x = x_unscaled, ymin = conf.low, ymax = conf.high, 
                  fill = group_clean), alpha = 0.5) +
  geom_line(aes(x = x_unscaled, y = predicted, color = group_clean), linewidth = 1) +
  facet_wrap(~var_clean, scales = "free_x", ncol = 6) +
  labs(y = "EVI Trend", x = "Predictor Value", color = "", fill = "") +
  scale_fill_manual(values = c("#DC9954", "#5B80BC")) +
  scale_color_manual(values = c("#DC9954", "#5B80BC")) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 22.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"))

p_pred

p_pa <- (p_pa_shapes / p_pred) + plot_layout(heights = c(3, 1)) +
  plot_annotation(tag_levels = 'A')
ggsave(plot = p_pa, "builds/plots/pas_map_n_preds.png", dpi = 900, height = 8, width = 10)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   Trend comparison   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#### Ridges -----
dt_long <- dt %>% pivot_longer(
  cols = c("nitrogen_depo", "mat_coef",
           "prec_coef", "hmi_change", 
           "max_temp_coef", "fire_frequency"), 
  names_to = "var_name", 
  values_to = "var_value") %>% 
  mutate(var_clean = case_when(
    grepl("nitrogen_depo", var_name) ~ "N Deposition",
    grepl("mat_coef", var_name) ~ "MAT Trend",
    grepl("prec_coef", var_name) ~ "Prec. Trend",
    grepl("max_temp_coef", var_name) ~ "Max. Temp. Trend",
    grepl("hmi_change", var_name) ~ "HMI Change",
    grepl("fire_frequency", var_name) ~ "Fire frequency"), 
    protection_cat_broad = ifelse(protection_cat_broad == "strict", "Protected", "Control")
  )


p_ridge = dt_long %>% 
  ggplot() +
  geom_density_ridges_gradient(aes(x = var_value, y = protection_cat_broad,
                                   fill = protection_cat_broad, color = protection_cat_broad), alpha = 0.7) +
  scale_fill_manual(values = c("#DC9954", "#5B80BC")) +
  scale_color_manual(values = c("#DC9954", "#5B80BC")) +
  facet_wrap(~var_clean, scale = "free_x", ncol = 2) +
  theme_minimal() +
  labs(y = "", x = "Variable Value") + 
  theme(legend.position = "nonw", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"))

p_ridge



#### cohens D -----------
setDT(dt)
dt_prot <- dt %>% 
  filter(protection_cat_broad == "strict") %>% 
  dplyr::select(nitrogen_depo_protected = nitrogen_depo, mat_coef_protected = mat_coef,
                prec_coef_protected =  prec_coef,hmi_change_protected = hmi_change, 
                max_temp_coef_protected = max_temp_coef, fire_frequency_protected  = fire_frequency,
                pair_id)

dt_control <- dt %>% 
  filter(!protection_cat_broad == "strict") %>% 
  dplyr::select(nitrogen_depo_control = nitrogen_depo, mat_coef_control = mat_coef,
                prec_coef_control =  prec_coef,hmi_change_control = hmi_change, 
                max_temp_coef_control = max_temp_coef, fire_frequency_control  = fire_frequency,
                pair_id)

dt_comp = dt_prot %>% 
  left_join(dt_control)

dt_comp = dt_prot %>% 
  left_join(dt_control)

glimpse(dt_comp)

protected_vars <- names(dt_comp)[grepl("_protected$", names(dt_comp))]
control_vars   <- sub("_protected$", "_control", protected_vars)

dt_res <- data.frame()

for(var in protected_vars){
  
  ctrl <- sub("_protected", "_control", var)
  general_var <- sub("_protected", "", var)
  
  
  # t-test
  tt <- t.test(dt_comp[[var]],
               dt_comp[[ctrl]],
               paired = TRUE)
  
  # Cohen's d
  d_out <- cohens_d(var, ctrl,
                    data = dt_comp,
                    paired = TRUE,
                    ci = 0.95)
  
  tmp <- data.frame(
    variable  = general_var,
    mean_diff = tt$estimate,
    ci_low = tt$conf.int[1],
    ci_high = tt$conf.int[2],
    p_value = round(tt$p.value,4),
    cohens_d = d_out$Cohens_d,
    d_ci_lb = d_out$CI_low,
    d_ci_ub = d_out$CI_high
  )
  
  dt_res <- rbind(dt_res, tmp)
}

interpret_cohens_d(dt_res$cohens_d, rules = "cohen1988")
interpret_cohens_d(dt_res$cohens_d, rules = "sawilowsky2009")
interpret_cohens_d(dt_res$cohens_d, rules = "gignac2016")
interpret_cohens_d(dt_res$cohens_d, rules = "lovakov2021")

#

p_d = dt_res %>% 
  mutate(var_clean = case_when(
    grepl("nitrogen_depo", variable) ~ "N Deposition",
    grepl("mat_coef", variable) ~ "MAT Trend",
    grepl("prec_coef", variable) ~ "Prec. Trend",
    grepl("max_temp_coef", variable) ~ "Max. Temp. Trend",
    grepl("hmi_change", variable) ~ "HMI Change",
    grepl("fire_frequency", variable) ~ "Fire Frequency"),
    sig = d_ci_lb > 0 | d_ci_ub < 0   # CI excludes 0
  ) %>%
  ggplot(aes(x = cohens_d, y = var_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = d_ci_lb,
                      xmax = d_ci_ub,
                      colour = sig),
                  shape = 18,
                  size = 1.4,
                  linewidth = 1, 
                  alpha = 0.75) +
  scale_fill_manual(values = c("grey50", "black"), guide = "none") +
  scale_colour_manual(values = c("grey50", "black"), guide = "none") +
  annotate("text", x = min(dt_res$cohens_d, na.rm=TRUE) - 0.05,
           y = 6.1,
           label = "Higher in\ncontrols",
           hjust = 0, size = 2, color = "grey50", fontface = "italic") +
  annotate("text", x = max(dt_res$cohens_d, na.rm=TRUE) + 0.05,
           y = 6.1,
           label = "Higher in\nprotected areas",
           hjust = 1, size = 2, color = "grey50", fontface = "italic") +
  
  labs(x = "Cohen's d", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
  )

p_d

# estimates ------

p_me <- dt_bm %>% 
  filter(tier == "prot_importance") %>% 
  ggplot() +
  geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high,
                      y = clean_term, color = clean_term),
                  shape = 18,
                  size = 1.4,
                  linewidth = 1, 
                  alpha = 0.9) +
  scale_fill_manual(values = c("#DC9954", "#5B80BC")) +
  scale_color_manual(values = c("#DC9954", "#5B80BC")) +
  labs(x = "Model estimate (± 95 % CI)", y = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank())
p_me

(p_est = (p_d / p_me) + plot_layout(heights = c(2, 1)))

p_difference <- (p_ridge | p_est) + 
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = 'A')
ggsave(plot = p_difference, "builds/plots/global_change_difference_pas.png",
       dpi = 900, height = 6, width = 9)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   CONTROL DISTANCE   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

count_plot <- dt %>%
  filter(control_within_dist != "na") %>% 
  mutate(control_within_dist = factor(control_within_dist, 
                                      levels = c("100k", "500k", "1000k", "over_1000k"),
                                      labels = c("≤ 100 km", "≤ 500 km", "≤ 1000 km", "> 1000 km"))) %>%
  ggplot(aes(x = control_within_dist)) +
  geom_bar() +
  labs(x = "Distance", y = "Count", title = "Distance between control and PA") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"))
count_plot

ggsave(plot = count_plot, "builds/plots/pa_control_within_dist_ditribution.png", dpi = 600, height = 4.5, width = 4.5)

