source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  as.data.frame()


dt_raw2 <- dt_raw %>% 
  mutate(
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"),
    pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
    nitrogen_depo = scale(nitrogen_depo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(human_modification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA)) 

dt_pa <- dt %>% 
  filter(protection_cat_broad == "strictly_protected")

for(pa_id in unique(dt_pa$unique_id)){

  dt[dt$control_for %in% pa_id, ]$pa_age <- dt[dt$unique_id %in% pa_id, ]$pa_age
  
  print(pa_id)
  }


dt <- dt %>% mutate(area_km2_log = scale(log(area_km2 + 0.0001)),
                     pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA))

### Mean evi
# recalculate trends only on PAs 
mean_evi_cols <- grep("mean_evi_", names(dt), value = T)

#subset to complete cases
dt_mean_evi <- dt %>% filter(protection_cat_broad == "control") %>%
  dplyr::select(all_of(mean_evi_cols),functional_biome, lon, lat, unique_id, 
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad, area_km2_log, pa_age_log) %>% 
  filter(complete.cases(.))

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
                    start = list(range = 0.1), fit.n = nrow(dt_mean_evi))
(range_opt_mean_evi = corfit_mean_evi$spcor)


#partition the data 

pm_mean_evi <- sample_partitions(npix = nrow(dt_pa), partsize = 1200, npart = NA)
dim(pm_mean_evi)

dt_mean_evi$lon

dt_mean_evi2 <- dt_mean_evi %>% slice_sample(n = 1000)

coords_mean_evi <- as.matrix(dt_mean_evi[, c("lon", "lat")])

# get distance matrix 
d_mean_evi <- distm_scaled(coords_mean_evi)

#use r to calculate optimal v parameter 
v_opt_mean_evi <- covar_exp(d_mean_evi, range_opt_mean_evi)

gls_mean_evi <- fitGLS(mean_evi_coef ~ 1 +
                                   area_km2_log + 
                                   pa_age_log,
                                 covar_FUN = "covar_exp",
                                 covar.pars = list(range = range_opt_mean_evi),
                                 data = dt_mean_evi,
                       V = v_opt_mean_evi, 
                       nugget = NA, 
                       no.F = TRUE)
gls_mean_evi 


dt_est_mean_evi <- extract_gls_estimates(gls_mean_evi, part = FALSE)

p_mean_evi <- dt_est_mean_evi %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "Controls: EVI change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_mean_evi

## Hypothesis 5.1 - absolute change PA characteristics ----------------------------
gls_mean_evi.1 <- fitGLS(abs_mean_evi_coef ~ 1 +
                           area_km2_log + 
                           pa_age_log,
                         covar_FUN = "covar_exp",
                         covar.pars = list(range = range_opt_mean_evi),
                         data = dt_mean_evi,
                         V = v_opt_mean_evi, 
                         nugget = NA, 
                         no.F = TRUE)
gls_mean_evi.1 

dt_est_mean_evi.1 <- extract_gls_estimates(gls_mean_evi.1, part = FALSE)

p_mean_evi.1 <- dt_est_mean_evi.1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H5: Abs EVI Change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_mean_evi.1

p_cont <- gridExtra::grid.arrange(p_mean_evi, p_mean_evi.1, ncol = 2)
ggsave(plot = p_cont, "builds/plots/pa_characteristics_in_controls_evi_change.png", height = 3, width = 6, dpi = 600)
