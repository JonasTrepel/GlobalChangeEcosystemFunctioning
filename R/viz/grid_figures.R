library(tidyverse)
library(data.table)
library(sf)
library(ggridges)
library(rnaturalearth)
library(scico)
library(effectsize)
library(patchwork)

#### Load all data ------------

dt <- fread("data/processed_data/analysis_ready_world_grid.csv")
dt_pred <-  fread("builds/model_outputs/world_grid_predictions.csv")
dt_bm = fread("builds/model_outputs/sdmtmb_results_world_grid.csv") %>% 
  mutate(subset_group = case_when(
    grepl("climatic_region", tier) ~ "climatic_region", 
    grepl("full_dataset", tier) ~ "full_dataset",
    grepl("functional_biome", tier) ~ "functional_biome",
    grepl("olson_biome", tier) ~ "olson_biome",
    grepl("super_biome", tier) ~ "super_biome"),
    subset = case_when(
      grepl("climatic_region", tier) ~ gsub("climatic_region_", "", tier), 
      grepl("full_dataset", tier) ~ gsub("climatic_region_", "", tier),
      grepl("functional_biome", tier) ~ gsub("functional_biome_", "", tier),
      grepl("olson_biome", tier) ~ gsub("olson_biome_", "", tier),
      grepl("super_biome", tier) ~ gsub("super_biome_", "", tier)), 
    subset_clean = case_when(
      .default = subset, 
      subset == "cold_short" ~ "Cold Limited\nShort Vegetation",
      subset == "cold_tall" ~ "Cold Limited\nTall Vegetation",
      subset == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
      subset == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation"), 
    sig = ifelse(conf.low > 0 | conf.high < 0, "significant", "not significiant")
    ) %>% 
  as.data.table()

ct_tiers = unique(dt_bm[subset_group == "climatic_region",]$tier)

climatic_grid_ids = data.frame()
for(t in ct_tiers){

  m = readRDS(unique(dt_bm[tier == t,]$model_path))
  dat <- m$data

  ids = unique(dat$unique_id)

  tmp = data.frame(unique_id = ids,
                   tier = t)


  climatic_grid_ids = rbind(climatic_grid_ids, tmp)

}

sb_tiers = unique(dt_bm[subset_group == "super_biome",]$tier)

sb_grid_ids = data.frame()
for(t in sb_tiers){
  
  m = readRDS(unique(dt_bm[tier == t,]$model_path))
  dat <- m$data
  
  ids = unique(dat$unique_id)
  
  tmp = data.frame(unique_id = ids,
                   tier = t)
  
  
  sb_grid_ids = rbind(sb_grid_ids, tmp)
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   Full Dataset   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
m <- readRDS(unique(dt_bm[dt_bm$tier == "full_dataset_yes", ]$model_path))
dt_fm <- m$data

quantile(shapes$mean_evi_coef, c(.025, .975), na.rm = T)

q_025_evi <- as.numeric(quantile(dt_fm$evi_coef, c(.025), na.rm = T))
q_975_evi <- as.numeric(quantile(dt_fm$evi_coef, c(.975), na.rm = T))

p_evi_map <- dt_fm %>% 
  mutate(evi_coef = case_when(
    .default = evi_coef, 
    evi_coef < quantile(dt_fm$evi_coef, c(.025), na.rm = T) ~ quantile(dt_fm$evi_coef, c(.025), na.rm = T), 
    evi_coef > quantile(dt_fm$evi_coef, c(.975), na.rm = T) ~ quantile(dt_fm$evi_coef, c(.975), na.rm = T), 
  )) %>% 
  ggplot() +
  geom_tile(aes(x = x_moll, y = y_moll, fill= evi_coef, color = evi_coef)) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  labs(color = "EVI\nTrend", fill = "EVI\nTrend") +
  theme_void() +
  theme(legend.position = "bottom",
        #legend.position.inside = c(0.1, 0.1),
        # legend.key.width = unit(1, "cm"),
        # legend.key.height = unit(0.4, "cm"), 
        legend.text = element_text(angle = 0)) +
  coord_fixed()
p_evi_map

p_pred <- dt_pred %>% 
  filter(tier == "full_dataset_yes") %>% 
  filter(x_unscaled > q05_unscaled,  x_unscaled < q95_unscaled) %>% 
  left_join(dt_bm %>% 
              select(tier, term, sig) %>% 
              unique()) %>% 
  ggplot() +
 # geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_ribbon(aes(x = x_unscaled, ymin = conf.low, ymax = conf.high),
              fill = "#662A3C", alpha = 0.25) +
  geom_line(aes(x = x_unscaled, y = predicted, 
                linetype = sig), linewidth = 1, color = "#662A3C") +
  facet_wrap(~var_clean, scales = "free_x", ncol = 6) +
  labs(y = "Evi Trend", x = "Predictor Value", color = "", fill = "") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"))

p_pred

p_main = (p_evi_map / p_pred) + plot_layout(heights = c(3, 1))
ggsave(plot = p_main, "builds/plots/evi_map_and_preds.png", 
       dpi = 900, height = 8, width = 10)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   Main Subsets   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
scico(palette = "managua", n = 5)
c("#FFCE66","#B16243","#572948","#5775B3","#80E6FF")

p_ct_map = dt %>%
  filter(unique_id %in% unique(climatic_grid_ids$unique_id)) %>% 
 # slice_sample(prop = .25) %>% 
  ggplot() +
  geom_tile(aes(x = x_moll, y = y_moll, fill= climatic_region, color = climatic_region)) +
  theme_void() +
  scale_color_manual(values = c("Arid" = "#FFCE66",
                                "Temperate" = "#B16243",
                                "Tropical" = "#572948",
                                "Cold" = "#5775B3",
                                "Polar" = "#80E6FF"
                                )) +
  scale_fill_manual(values = c("Arid" = "#FFCE66",
                                "Temperate" = "#B16243",
                                "Tropical" = "#572948",
                                "Cold" = "#5775B3",
                                "Polar" = "#80E6FF"
  )) +
  theme(legend.position = "bottom") +
  labs(color = "", fill = "") +
  coord_fixed()
p_ct_map

p_ct_est <- dt_bm %>% 
  filter(subset_group == "super_biome" & !clean_term == "(Intercept)") %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = subset, group = subset), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85,
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Arid" = "#FFCE66",
                                "Temperate" = "#B16243",
                                "Tropical" = "#572948",
                                "Cold" = "#5775B3",
                                "Polar" = "#80E6FF"
  )) + 
  theme_minimal() +
  labs(y = "") +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
p_ct_est



scico(palette = "batlow", n = 10)
c("#001959", "#0F3F60", "#1B5961", "#3B6C55", "#687A3D", "#9C892A", "#D29243",
  "#F8A17B","#FDB6BB", "#F9CCF9")

p_sb_map = dt %>%
  filter(unique_id %in% unique(sb_grid_ids$unique_id)) %>% 
  mutate(    super_biome = case_when(
    super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
    super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
    super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
    super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation")) %>% 
  # slice_sample(prop = .25) %>% 
  ggplot() +
  geom_tile(aes(x = x_moll, y = y_moll, fill= super_biome, color = super_biome)) +
  theme_void() +
  scale_color_manual(values = c("Cold Limited\nShort Vegetation" = "#9C892A",
                                "Cold Limited\nTall Vegetation" = "#0F3F60",
                                "Not Cold Limited\nShort Vegetation" = "#F8A17B",
                                "Not Cold Limited\nTall Vegetation" = "#3B6C55"
  )) +
  scale_fill_manual(values = c("Cold Limited\nShort Vegetation" = "#9C892A",
                               "Cold Limited\nTall Vegetation" = "#0F3F60",
                               "Not Cold Limited\nShort Vegetation" = "#F8A17B",
                               "Not Cold Limited\nTall Vegetation" = "#3B6C55"
  )) +
  theme(legend.position = "bottom") +
  labs(color = "", fill = "") +
  coord_fixed()
p_sb_map

p_sb_est <- dt_bm %>% 
  filter(subset_group == "super_biome" & !clean_term == "(Intercept)") %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = subset_clean, group = subset_clean), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85,
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Cold Limited\nShort Vegetation" = "#9C892A",
                                "Cold Limited\nTall Vegetation" = "#0F3F60",
                                "Not Cold Limited\nShort Vegetation" = "#F8A17B",
                                "Not Cold Limited\nTall Vegetation" = "#3B6C55"
  )) + 
  theme_minimal() +
  labs(y = "") +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
p_sb_est

p_ct <- (p_ct_est | p_ct_map) + plot_layout(widths = c(1, 2))
p_ct

p_sb <- (p_sb_est | p_sb_map) + plot_layout(widths = c(1, 2))
p_sb

p_sub <- p_ct/p_sb
p_sub
ggsave(plot = p_sub, "builds/plots/main_subsets.png", dpi = 900, 
       height = 9, width = 10)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################   Smaller Subsets   ##################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
scico(palette = "bam", n = 9)
c("#65014B","#A4428B", "#D07EBB", "#EBC5E1" ,"#F5F0F0", "#D7E7C0", "#8CB464", "#4B802E" ,"#0C4C00")

p_fb_est <- dt_bm %>% 
  filter(subset_group == "functional_biome" & !clean_term == "(Intercept)") %>%
  mutate(sig = (conf.low > 0 | conf.high < 0)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = sig), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85) +
  scale_color_manual(values = c("#D7E7C0", "#0C4C00")) + 
  facet_wrap(~subset, scales = "free_x", ncol = 5) +
  theme_minimal() +
  labs(y = "", color = "significance") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 22.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"), 
        strip.text = element_text(face = "bold"))
p_fb_est
ggsave(plot = p_fb_est, "builds/plots/supplement/functional_biomes_estimates.png", 
       dpi = 900, height = 7.5, width = 10)

p_ob_est <- dt_bm %>% 
  filter(subset_group == "olson_biome" & !clean_term == "(Intercept)") %>%
  mutate(sig = (conf.low > 0 | conf.high < 0)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = sig), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85) +
  scale_color_manual(values = c("#D7E7C0", "#0C4C00")) + 
  facet_wrap(~subset, scales = "free_x", ncol = 3) +
  theme_minimal() +
  labs(y = "", color = "significance") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 22.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"), 
        strip.text = element_text(face = "bold"))
p_ob_est

ggsave(plot = p_ob_est, "builds/plots/supplement/olson_biome_estimates.png", 
       dpi = 900, height = 10, width = 10)
