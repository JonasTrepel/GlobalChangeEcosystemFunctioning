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


tiers <- "prot_global_change"

vars = c("nitrogen_depo_scaled", 
         "prec_coef_scaled", 
         "hmi_change_scaled",
         "max_temp_coef_scaled", 
         "fire_frequency_scaled", 
         "mat_coef_scaled")

extr_guide <- CJ(tier = tiers, 
                      vars = vars) %>% 
  left_join(dt_bm_subset)

plan(multisession, workers = 6)

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
    
    m_plot <- ggeffects::ggpredict(m, terms = c(var_call, "protection_cat_broad"), 
                                   type = "fixed")

    
    plot_data <- as.data.table(m_plot) %>%
      mutate(
        x_unscaled = round(x * sd_x + mean_x, 7),
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
print(paste0("subset done ", Sys.time()))

dt_pred_comp <- rbindlist(for_results_pred) %>% 
  mutate(var_clean = case_when(
    grepl("nitrogen_depo_scaled", var_name) ~ "Nitrogen Deposition",
    grepl("mat_coef_scaled", var_name) ~ "MAT Trend",
    grepl("prec_coef_scaled", var_name) ~ "Precipitation Trend",
    grepl("max_temp_coef_scaled", var_name) ~ "Max. Temperature Trend",
    grepl("hmi_change_scaled", var_name) ~ "HMI Change",
    grepl("fire_frequency_scaled", var_name) ~ "Fire frequency"),
    group_clean = case_when(
           grepl("strict", group) ~  "Protected", 
           grepl("unprotected", group) ~  "Control"))
unique(dt_pred_comp$response_name)
unique(dt_pred_comp$var_name)

fwrite(dt_pred_comp, "builds/model_outputs/pa_dataset_global_change_predictions.csv")

###### Plot - 1 facte_grid for each response, comparing different tiers? 
m <- readRDS(extr_guide[1, ]$model_path)
dat <- m$data

dt_long <- dat %>% pivot_longer(
  cols = c("nitrogen_depo", "mat_coef",
           "prec_coef", "hmi_change", 
           "max_temp_coef", "fire_frequency"), 
  names_to = "var_name", 
  values_to = "var_value") %>% 
  mutate(var_clean = case_when(
    grepl("nitrogen_depo", var_name) ~ "Nitrogen Deposition",
    grepl("mat_coef", var_name) ~ "MAT Trend",
    grepl("prec_coef", var_name) ~ "Precipitation Trend",
    grepl("max_temp_coef", var_name) ~ "Max. Temperature Trend",
    grepl("hmi_change", var_name) ~ "HMI Change",
    grepl("fire_frequency", var_name) ~ "Fire frequency")
  ) %>% 
  left_join(dt_pred_comp %>% 
              dplyr::select(q05_unscaled, q95_unscaled, var_clean) %>% 
              unique())

p_b <- dt_pred_comp %>% 
  #filter(x_unscaled > q05_unscaled,  x_unscaled < q95_unscaled) %>% 
  ggplot() +
  geom_point(data = dt_long, aes(x = var_value, y = evi_coef), alpha = 0.1, size = 0.1, color = "grey25") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_ribbon(aes(x = x_unscaled, ymin = conf.low, ymax = conf.high, 
                  fill = group_clean), alpha = 0.25) +
  geom_line(aes(x = x_unscaled, y = predicted, color = group_clean), linewidth = 1) +
  facet_wrap(~var_clean, scales = "free_x", ncol = 6) +
  labs(y = "Response Value", x = "Predictor Value", color = "", fill = "", title = "B") +
  scale_fill_manual(values = c("#DC9954", "#5B80BC")) +
  scale_color_manual(values = c("#DC9954", "#5B80BC")) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"))

p_b

p_a <- dt_pred_comp %>% 
  #filter(x_unscaled > q05_unscaled,  x_unscaled < q95_unscaled) %>% 
  ggplot() +
 # geom_point(data = dt_long, aes(x = var_value, y = evi_coef), alpha = 0.1, size = 0.1, color = "grey25") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25") +
  geom_ribbon(aes(x = x_unscaled, ymin = conf.low, ymax = conf.high, 
                  fill = group_clean), alpha = 0.25) +
  geom_line(aes(x = x_unscaled, y = predicted, color = group_clean), linewidth = 1) +
  facet_wrap(~var_clean, scales = "free_x", ncol = 6) +
  labs(y = "Response Value", x = "Predictor Value", color = "", fill = "", title = "A") +
  scale_fill_manual(values = c("#DC9954", "#5B80BC")) +
  scale_color_manual(values = c("#DC9954", "#5B80BC")) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        strip.background = element_rect(fill = "snow", color = "snow"))

p_a

p_q = p_a / p_b

ggsave(plot = p_q, "builds/plots/supplement/pa_dataset_predictions_with_points.png",
       dpi = 600, height = 6, width = 12)

