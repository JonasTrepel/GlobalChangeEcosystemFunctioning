source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") %>% 
  as.data.frame()


dt_corr <- dt_raw %>%
  dplyr::select(mean_n_depo_zhu, nitrogen_depo, mean_evi) %>%
  filter(complete.cases(.))

cor.test(dt_corr$mean_n_depo_zhu, dt_corr$nitrogen_depo)
cor.test(dt_corr$mean_evi, dt_corr$nitrogen_depo)
cor.test(dt_corr$mean_n_depo_zhu, dt_corr$mean_evi)

p_corr_1 <- dt_corr %>%
  ggplot() +
  geom_point(aes(x = nitrogen_depo, y = mean_n_depo_zhu*100), alpha = 0.25) +
  labs(x = "N Deposition (Rubin et al. 2023)", y = "N Deposition (Zhu et al. 2025)") +
  geom_smooth(aes(x = nitrogen_depo, y = mean_n_depo_zhu*100), color = "blue", linewidth = 1.1)  +
  geom_abline(color = "red", linewidth = 1.1)  +
  theme_classic()
p_corr_1

p_corr_2 <- dt_corr %>%
  ggplot() +
  geom_point(aes(x = mean_evi/10000, y = nitrogen_depo), alpha = 0.25) +
  labs(x = "EVI", y = "N Deposition (Rubin et al. 2023)") +
  geom_smooth(aes(x = mean_evi/10000, y = nitrogen_depo), color = "blue", linewidth = 1.1)  +
  theme_classic()
p_corr_2

p_corr_3 <- dt_corr %>%
  ggplot() +
  geom_point(aes(x =  mean_evi/10000, y = mean_n_depo_zhu*100), alpha = 0.25) +
  labs(x = "EVI", y = "N Deposition (Zhu et al. 2025)") +
  geom_smooth(aes(x = mean_evi/10000, y = mean_n_depo_zhu*100), color = "blue", linewidth = 1.1)  +
  theme_classic()
p_corr_3

p_corr <- gridExtra::grid.arrange(p_corr_1, p_corr_2, p_corr_3, ncol = 3)
ggsave(plot = p_corr, "builds/plots/supplement/n_depo_corr.png", dpi = 600, height = 3, width = 9)


dt <- dt_raw %>% 
  mutate(
    pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
    nitrogen_depo = scale(nitrogen_depo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(human_modification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA),
    mean_n_depo_zhu = scale(mean_n_depo_zhu)) 

## EVI ------------------------------------------------------------------------------
####### calculate mean_evi trend #######

mean_evi_cols <- grep("mean_evi_", names(dt), value = T)

#subset to complete cases
dt_mean_evi <- dt %>%
  dplyr::select(all_of(mean_evi_cols), functional_biome, X, Y, lon, lat, unique_id, mean_n_depo_zhu,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.))

table(dt_mean_evi$protection_cat_broad)
table(dt_mean_evi$super_biome)


#get y matrix
y_mean_evi <- as.matrix(dt_mean_evi[, mean_evi_cols])

#get coordinate matrix 
coords_mean_evi <- as.matrix(dt_mean_evi[, c("lon", "lat")])

#fit autoregression, accounting for temporal autocorrelation 
ar_mean_evi <- fitAR_map(Y = y_mean_evi, coords = coords_mean_evi)

#extract coefficients and p values 
dt_mean_evi$mean_evi_coef <- ar_mean_evi$coefficients[, "t"]
dt_mean_evi$abs_mean_evi_coef <- abs(ar_mean_evi$coefficients[, "t"])
dt_mean_evi$mean_evi_p_value <- ar_mean_evi$pvals[, 2]


### estimate optimal r parameter (range of spatial autocorrelation)
corfit_mean_evi <- fitCor(resids = residuals(ar_mean_evi), coords = coords_mean_evi, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = 15000)
(range_opt_mean_evi = corfit_mean_evi$spcor) #0.01642254 


############ test hypotheses ############ 

#partition the data 
set.seed(161)
pm_mean_evi <- sample_partitions(npix = nrow(dt_mean_evi), partsize = 1500, npart = NA)
dim(pm_mean_evi) 


## Hypothesis 2: Change depends on climate change, N deposition, human modification -------------

set.seed(161)
gls_mean_evi <- fitGLS_partition(mean_evi_coef ~ 1 +
                                   mean_n_depo_zhu +
                             mat_coef + 
                             max_temp_coef + 
                             map_coef + 
                             human_modification,
                           partmat = pm_mean_evi,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_mean_evi),
                           data = dt_mean_evi,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("lon", "lat"))
gls_mean_evi 

dt_est_mean_evi <- extract_gls_estimates(gls_mean_evi, part = TRUE)

p_mean_evi <- dt_est_mean_evi %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "H2: mean_evi change ~\nglobal change", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_mean_evi


## BURNED AREA ------------------------------------------------------------------------------

####### calculate burned_area trend #######

burned_area_cols <- grep("burned_area_", names(dt), value = T)

#subset to complete cases
dt_burned_area <- dt %>%
  dplyr::select(all_of(burned_area_cols), functional_biome, X, Y, lon, lat, unique_id, mean_n_depo_zhu,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.))

table(dt_burned_area$protection_cat_broad)
table(dt_burned_area$super_biome)


#get y matrix
y_burned_area <- as.matrix(dt_burned_area[, burned_area_cols])

#get coordinate matrix 
coords_burned_area <- as.matrix(dt_burned_area[, c("lon", "lat")])

#fit autoregression, accounting for temporal autocorrelation 
ar_burned_area <- fitAR_map(Y = y_burned_area, coords = coords_burned_area)

#extract coefficients and p values 
dt_burned_area$burned_area_coef <- ar_burned_area$coefficients[, "t"]
dt_burned_area$abs_burned_area_coef <- abs(ar_burned_area$coefficients[, "t"])
dt_burned_area$burned_area_p_value <- ar_burned_area$pvals[, 2]


### estimate optimal r parameter (range of spatial autocorrelation)
corfit_burned_area <- fitCor(resids = residuals(ar_burned_area), coords = coords_burned_area, covar_FUN = "covar_exp", 
                          start = list(range = 0.1), fit.n = 15000)
(range_opt_burned_area = corfit_burned_area$spcor) #0.001645042 


############ test hypotheses ############ 

#partition the data 
set.seed(161)
pm_burned_area <- sample_partitions(npix = nrow(dt_burned_area), partsize = 1500, npart = NA)
dim(pm_burned_area)


## Hypothesis 2: Change depends on climate change, N deposition, human modification -------------

set.seed(161)
gls_burned_area <- fitGLS_partition(burned_area_coef ~ 1 +
                                   mean_n_depo_zhu +
                                   mat_coef + 
                                   max_temp_coef + 
                                   map_coef + 
                                   human_modification,
                                 partmat = pm_burned_area,
                                 covar_FUN = "covar_exp",
                                 covar.pars = list(range = range_opt_burned_area),
                                 data = dt_burned_area,
                                 nugget = NA,
                                 ncores = 25,
                                 progressbar = TRUE, 
                                 parallel = TRUE, 
                                 coord.names = c("lon", "lat"))
gls_burned_area 

dt_est_burned_area <- extract_gls_estimates(gls_burned_area, part = TRUE)

p_burned_area <- dt_est_burned_area %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H2: burned_area change ~\nglobal change", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_burned_area


## GREEN-UP ------------------------------------------------------
####### calculate mean_evi trend #######

greenup_cols <- grep("greenup_", names(dt), value = T)

#subset to complete cases
dt_greenup <- dt %>%
  filter(!grepl("N", functional_biome)) %>% 
  dplyr::select(all_of(greenup_cols), functional_biome, X, Y, lon, lat, unique_id, mean_n_depo_zhu,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.))

table(dt_greenup$protection_cat_broad)
table(dt_greenup$super_biome)


#get y matrix
y_greenup <- as.matrix(dt_greenup[, greenup_cols])

#get coordinate matrix 
coords_greenup <- as.matrix(dt_greenup[, c("lon", "lat")])

#fit autoregression, accounting for temporal autocorrelation 
ar_greenup <- fitAR_map(Y = y_greenup, coords = coords_greenup)

#extract coefficients and p values 
dt_greenup$greenup_coef <- ar_greenup$coefficients[, "t"]
dt_greenup$abs_greenup_coef <- abs(ar_greenup$coefficients[, "t"])
dt_greenup$greenup_p_value <- ar_greenup$pvals[, 2]


### estimate optimal r parameter (range of spatial autocorrelation)
corfit_greenup <- fitCor(resids = residuals(ar_greenup), coords = coords_greenup, covar_FUN = "covar_exp", 
                          start = list(range = 0.1), fit.n = 15000)
(range_opt_greenup = corfit_greenup$spcor)


############ test hypotheses ############ 

#partition the data 
set.seed(161)
pm_greenup <- sample_partitions(npix = nrow(dt_greenup), partsize = 1500, npart = NA)
dim(pm_greenup)


## Hypothesis 2: Change depends on climate change, N deposition, human modification -------------

set.seed(161)
gls_greenup <- fitGLS_partition(greenup_coef ~ 1 +
                                   mean_n_depo_zhu +
                                   mat_coef + 
                                   max_temp_coef + 
                                   map_coef + 
                                   human_modification,
                                 partmat = pm_greenup,
                                 covar_FUN = "covar_exp",
                                 covar.pars = list(range = range_opt_greenup),
                                 data = dt_greenup,
                                 nugget = NA,
                                 ncores = 25,
                                 progressbar = TRUE, 
                                 parallel = TRUE, 
                                 coord.names = c("lon", "lat"))
gls_greenup 

dt_est_greenup <- extract_gls_estimates(gls_greenup, part = TRUE)

p_greenup <- dt_est_greenup %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkcyan", "non-significant" = "grey")) +
  labs(title = "H2: greenup change ~\nglobal change", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_greenup

# Combine -----------------------------------

dt_est <- rbind(
  dt_est_mean_evi %>% mutate(model = "zhu_n_depo", t.stat = NA, response = "mean_evi"),
  dt_est_burned_area %>% mutate(model = "zhu_n_depo", t.stat = NA, response = "burned_area"), 
  dt_est_greenup %>% mutate(model = "zhu_n_depo", t.stat = NA, response = "greenup")
)

fwrite(dt_est, "builds/model_estimates/global_change_models_with_ndep_zhu.csv")
# Plot --------------------------------

dt_plot <- dt_est %>%
  mutate(clean_term = case_when(
    grepl("Intercept", term) ~ "Intercept", 
    grepl("mean_n_depo_zhu", term) ~ "Nitrogen Deposition (Zhu)",
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
      response == "mean_evi" ~ "EVI Trend", 
      response == "burned_area" ~ "Burned Area\nTrend", 
      response == "greenup" ~ "Vegetation Green-up\nTrend"
    ), 
    ci_lb = case_when(
      response == "mean_evi" ~ ci_lb/100,
      response == "burned_area" ~ ci_lb*100,
      response == "greenup" ~ ci_lb),
    ci_ub = case_when(
      response == "mean_evi" ~ ci_ub/100,
      response == "burned_area" ~ ci_ub*100,
      response == "greenup" ~ ci_ub),
    estimate = case_when(
      response == "mean_evi" ~ estimate/100,
      response == "burned_area" ~ estimate*100,
      response == "greenup" ~ estimate)
  ) 

scico::scico(palette = "roma", n = 10)
#"#7E1700" "#995215" "#AF7F2A" "#C7B354" "#CFE3A3" "#A4E5D2" "#5DC0D2" "#3191C1" "#2064AE" "#023198"

p_est <- dt_plot %>%
  mutate(clean_term = factor(clean_term, levels = c(
    "Intercept", 
    "MAP Trend", 
    "Max Temp Trend", 
    "MAT Trend", 
    "Human Modification",
    "Nitrogen Deposition (Zhu)"))) %>%
  mutate(facet_label = factor(facet_label, levels = c(
    "EVI Trend", 
    "Burned Area\nTrend", 
    "Vegetation Green-up\nTrend"))) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~facet_label, scales = "free_x", ncol = 5) +
  scale_color_manual(values = c("Sig. Negative" = "#3191C1",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#AF7F2A"
  )) +
  labs(title = "Nitrogen Deposition by Zhu et al., 2025", y = NULL, x = "Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_est
ggsave(plot = p_est, "builds/plots/supplement/global_models_with_zhu_ndep.png", dpi = 600, height = 3, width = 9)
