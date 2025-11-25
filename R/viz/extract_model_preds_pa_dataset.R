library(data.table)
library(patchwork)
library(tidyverse)
library(sdmTMB)
library(patchwork)
library(groupdata2)
library(future)
library(furrr)
library(ggeffects)

#### extract predictions for protected area dataset 

#### subsets ----------

dt_bm_subset <- fread("builds/model_outputs/sdmtmb_results_pa_dataset.csv") %>% 
  dplyr::select(tier, model_id, model_path, response, dev_explained_full, dev_explained_fixed) %>% 
  unique()



#### Extract variable specific predictions -----


tiers <- "a"

vars = c("nitrogen_depo_scaled", 
         "prec_coef_scaled", 
         "hmi_change_scaled",
         "max_temp_coef_scaled", 
         "fire_frequency_scaled", 
         "mat_coef_scaled")

extr_guide_prot <- CJ(tier = tiers, 
                      vars = vars) %>% 
  left_join(dt_bm_subset)

plan(multisession, workers = 15)

for_results_pred <- future_map(
  1:nrow(extr_guide),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE),
  function(i) {
    
    response <- extr_guide[i, ]$response
    tier <- unique(extr_guide[i, ]$tier)
    var <- unique(extr_guide[i, ]$var)
    term <- unique(extr_guide[i, ]$term)
    
    
    
    m <- readRDS(extr_guide[i, ]$model_path)
    dat <- m$data
    
    
    var_us <- gsub("_scaled", "", var)
    mean_x <- mean(dat[[var_us]], na.rm = TRUE)
    sd_x   <- sd(dat[[var_us]], na.rm = TRUE)
    
    var_call <- paste0(var, " [all]")
    
    m_plot <- ggeffects::ggpredict(m, terms = c(var_call, "protection_cat_broad"))
    
    plot_data <- as.data.table(m_plot) %>%
      mutate(
        x_unscaled = round(x * sd_x + mean_x, 3),
        var_name = var,
        response_name = response,
        tier = tier,
        dev_explained_full = extr_guide[i,]$dev_explained_full,
        dev_explained_var = extr_guide[i,]$dev_explained_var,
        n = nrow(dat),
        q95_unscaled = as.numeric(quantile(dat[[var_us]], .95, na.rm = T)), 
        q05_unscaled = as.numeric(quantile(dat[[var_us]], .05, na.rm = T)), 
        q95 = as.numeric(quantile(dat[[var]], .95, na.rm = T)), 
        q05 = as.numeric(quantile(dat[[var]], .05, na.rm = T))
      )
    
    # Ensure confidence interval columns exist
    if (!any(grepl("conf", names(plot_data)))) {
      plot_data <- plot_data %>%
        mutate(conf.low = NA,
               conf.high = NA,
               std.error = NA)
    }
    
    
    #print(paste0(response, " done at: ", Sys.time()))
    
    # Clean up that mess
    rm(m)
    gc()
    
    return(plot_data)
  }
)

plan(sequential)
print(paste0("subsetp done ", Sys.time()))

dt_pred_comp <- rbindlist(for_results_pred) %>% 
  mutate(scale = "km2", 
         response_clean = case_when(
           response_name == "canopy_height_900m_coef" ~ "Canopy Height Trend",
           response_name == "tree_cover_1000m_coef" ~ "Tree Cover Trend"
         ),
         var_clean = case_when(
           var_name == "local_density_km2_scaled" ~ "Local Elephant Density",
           var_name == "months_extreme_drought_scaled" ~ "N Drought Months",
           var_name == "fire_frequency_scaled" ~ "Fire Frequency",
           var_name == "mat_coef_scaled" ~ "Temperature Change",
           var_name == "n_deposition_scaled" ~ "N Deposition",
         ),
         tier_clean = case_when(
           grepl("_lq", tier) ~  "Lower Third", 
           grepl("_mq", tier) ~  "Middle Third", 
           grepl("_uq", tier) ~  "Upper Third", 
           grepl("_unfenced", tier) ~  "Unfenced", 
           grepl("_fenced", tier) ~  "Fenced"
         ), 
         response_tier = paste0(response_clean, "\n", tier_clean))
unique(dt_pred_comp$response_name)
unique(dt_pred_comp$var_name)

fwrite(dt_pred_comp, "builds/model_outputs/subset_predictions_1000m.csv")

###### Plot - 1 facte_grid for each response, comparing different tiers? 


p_q <- dt_pred_comp %>% 
  filter(tier_clean %in% c("Lower Third", "Middle Third", "Upper Third")) %>% 
  ggplot() +
  # geom_point(data = dt_long, aes(x = var_value, y = response_value), alpha = 0.1, size = 0.1, color = "grey25") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_ribbon(aes(x = x_unscaled, ymin = conf.low, ymax = conf.high, 
                  fill = response_clean), alpha = 0.25) +
  geom_line(aes(x = x_unscaled, y = predicted, color = response_clean), linewidth = 1) +
  facet_grid(rows = vars(response_tier), cols = vars(var_clean), scales = "free") +
  labs(y = "Response Value", x = "Predictor Value") +
  scico::scale_color_scico_d(begin = .2, end = .8, palette = "batlow") +
  scico::scale_fill_scico_d(begin = .2, end = .8, palette = "batlow") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "snow"), 
        strip.background = element_rect(fill = "linen", color = "linen"))

p_q

ggsave(plot = p_q, "builds/plots/supplement/1000m_predictions_different_starting_conditions.png",
       dpi = 600, height = 10, width = 10)


p_f <- dt_pred_comp %>% 
  filter(!tier_clean %in% c("Lower Third", "Middle Third", "Upper Third")) %>% 
  ggplot() +
  # geom_point(data = dt_long, aes(x = var_value, y = response_value), alpha = 0.1, size = 0.1, color = "grey25") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_ribbon(aes(x = x_unscaled, ymin = conf.low, ymax = conf.high, 
                  fill = response_clean), alpha = 0.25) +
  geom_line(aes(x = x_unscaled, y = predicted, color = response_clean), linewidth = 1) +
  facet_grid(rows = vars(response_tier), cols = vars(var_clean), scales = "free") +
  labs(y = "Response Value", x = "Predictor Value") +
  scico::scale_color_scico_d(begin = .2, end = .8, palette = "batlow") +
  scico::scale_fill_scico_d(begin = .2, end = .8, palette = "batlow") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "snow"), 
        strip.background = element_rect(fill = "linen", color = "linen"))

p_f

ggsave(plot = p_f, "builds/plots/supplement/1000m_predictions_fences.png",
       dpi = 600, height = 7, width = 10)