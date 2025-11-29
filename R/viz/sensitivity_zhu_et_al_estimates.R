#### Visualize Zhu N depo 

library(data.table)
library(tidyverse)


dt_bm <- fread("builds/model_outputs/sdmtmb_results_world_grid_zhu_sensitivity.csv")  %>% 
  mutate(subset_group = case_when(
    grepl("climatic_region", tier) ~ "climatic_region", 
    grepl("full_dataset", tier) ~ "full_dataset",
    grepl("functional_biome", tier) ~ "functional_biome",
    grepl("olson_biome", tier) ~ "olson_biome",
    grepl("super_biome", tier) ~ "super_biome"),
    subset = case_when(
      grepl("climatic_region", tier) ~ gsub("climatic_region_", "", tier), 
      grepl("full_dataset", tier) ~ "Full Dataset",
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

scico::scico(palette = "lajolla", n = 10)


p_full_est <- dt_bm %>% 
  filter(subset_group == "full_dataset" & !clean_term == "(Intercept)") %>%
  mutate(sig = (conf.low > 0 | conf.high < 0)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = sig), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85) +
  scale_color_manual(values = c("#8E3F3D")) + 
  facet_wrap(~subset, scales = "free_x", ncol = 5) +
  theme_minimal() +
  labs(y = "", color = "significance", title = "Alternative N Deposition Product", 
       subtitle = "Zhu et al. (2025), Nature Communications") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 22.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"), 
        strip.text = element_text(face = "bold"))
p_full_est


p_cr_est <- dt_bm %>% 
  filter(subset_group == "climatic_region" & !clean_term == "(Intercept)") %>%
  mutate(sig = (conf.low > 0 | conf.high < 0)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = sig), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85) +
  scale_color_manual(values = c("#EEB554", "#8E3F3D")) + 
  facet_wrap(~subset, scales = "free_x", ncol = 5) +
  theme_minimal() +
  labs(y = "", color = "significance") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 22.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"), 
        strip.text = element_text(face = "bold"))
p_cr_est

p_sb_est <- dt_bm %>% 
  filter(subset_group == "super_biome" & !clean_term == "(Intercept)") %>%
  mutate(sig = (conf.low > 0 | conf.high < 0)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = sig), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85) +
  scale_color_manual(values = c("#EEB554", "#8E3F3D")) + 
  facet_wrap(~subset_clean, scales = "free_x", ncol = 5) +
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
p_sb_est

library(patchwork)
p_zhu_comb = (p_full_est / p_cr_est / p_sb_est) + plot_annotation(tag_levels = "A") 
p_zhu_comb
ggsave(plot = p_zhu_comb, "builds/plots/supplement/zhu_n_depo_main_tiers_estimates.png", 
       dpi = 900, height = 8, width = 10)


p_fb_est <- dt_bm %>% 
  filter(subset_group == "functional_biome" & !clean_term == "(Intercept)") %>%
  mutate(sig = (conf.low > 0 | conf.high < 0)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = sig), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85) +
  scale_color_manual(values = c("#EEB554", "#8E3F3D")) + 
  facet_wrap(~subset, scales = "free_x", ncol = 5) +
  theme_minimal() +
  labs(y = "", color = "significance", title = "Alternative N Deposition Product", 
       subtitle = "Zhu et al. (2025), Nature Communications") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 22.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"), 
        strip.text = element_text(face = "bold"))
p_fb_est

ggsave(plot = p_fb_est, "builds/plots/supplement/zhu_n_depo_functional_biomes_estimates.png", 
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
  scale_color_manual(values = c("#EEB554", "#8E3F3D")) + 
  facet_wrap(~subset, scales = "free_x", ncol = 3) +
  theme_minimal() +
  labs(y = "", color = "significance", title = "Alternative N Deposition Product", 
       subtitle = "Zhu et al. (2025), Nature Communications") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 22.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"), 
        strip.text = element_text(face = "bold"))
p_ob_est

ggsave(plot = p_ob_est, "builds/plots/supplement/zhu_n_depo_olson_biomes_estimates.png", 
       dpi = 900, height = 10, width = 10)
