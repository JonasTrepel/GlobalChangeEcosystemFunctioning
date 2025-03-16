source("R/functions/extract_gls_estimates.R")


####### plot USA N depo #######
library(data.table)
library(sf)
library(tidyverse)
library(remotePARTS)

dt_raw <- fread("data/processedData/dataFragments/grid_usa_with_climate_trends.csv") %>% 
  distinct(lon, lat, .keep_all = TRUE)
  
grid_usa <- st_read("data/spatialData/grid_usa.gpkg") %>% 
  dplyr::select(geom, unique_id) %>% 
  left_join(dt_raw) %>% 
  filter(!is.na(nitrogen_depo))



p_mean <- grid_usa %>%
  filter(mean_n_depo_usa > as.numeric(quantile(grid_usa$mean_n_depo_usa, 0.025, na.rm = T)) &
           mean_n_depo_usa < as.numeric(quantile(grid_usa$mean_n_depo_usa, 0.975, na.rm = T))
  ) %>% 
  ggplot() +
  geom_sf(aes(color = mean_n_depo_usa, fill = mean_n_depo_usa)) + 
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  labs(fill = "N depo mean", color = "N depo mean")+
  theme_void()
p_mean

p_trend <- grid_usa %>%
  filter(n_depo_usa_coef > as.numeric(quantile(grid_usa$n_depo_usa_coef, 0.025, na.rm = T)) &
           n_depo_usa_coef < as.numeric(quantile(grid_usa$n_depo_usa_coef, 0.975, na.rm = T))
  ) %>% 
  ggplot() +
  geom_sf(aes(color = n_depo_usa_coef, fill = n_depo_usa_coef)) + 
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  labs(fill = "N depo trend", color = "N depo trend") +
  theme_void()
p_trend



dt_corr <- dt_raw %>%
  dplyr::select(mean_n_depo_usa, n_depo_usa_coef, nitrogen_depo,
                mean_n_depo_zhu, n_depo_zhu_coef,
                n_depo_usa_2010) %>%
  filter(complete.cases(.))

cor.test(dt_corr$mean_n_depo_zhu, dt_corr$nitrogen_depo)
cor.test(dt_corr$mean_n_depo_zhu, dt_corr$mean_n_depo_usa)
cor.test(dt_corr$mean_n_depo_usa, dt_corr$nitrogen_depo)
cor.test(dt_corr$n_depo_usa_2010, dt_corr$nitrogen_depo)
cor.test(dt_corr$mean_n_depo_usa, dt_corr$n_depo_usa_coef)
cor.test(dt_corr$n_depo_zhu_coef, dt_corr$n_depo_usa_coef)



p_r1 <- dt_raw %>%
  ggplot() + 
  geom_point(aes(x = mean_n_depo_usa, y = n_depo_usa_coef), alpha = 0.5, size = 0.7) +
  labs(x = "USA N depo mean", y = "USA N depo trend") +
  theme_classic()
p_r1

p_r2 <- dt_raw %>% ggplot() + 
  geom_point(aes(x = mean_n_depo_usa, y = mean_n_depo_zhu), alpha = 0.5, size = 0.7) +
  labs(x = "USA N depo mean", y = "Zhu et al N depo") +
  geom_abline(color = "red") +
  theme_classic()
p_r2

p_r3 <- dt_raw %>% ggplot() + 
  geom_point(aes(x = mean_n_depo_usa, y = nitrogen_depo/100), alpha = 0.5, size = 0.7) +
  labs(x = "USA N depo mean", y = "Rubin et al N depo") +
  geom_abline(color = "red") +
  theme_classic()
p_r3

p_r4 <- dt_raw %>% ggplot() + 
  geom_point(aes(x = nitrogen_depo/100, y = mean_n_depo_zhu), alpha = 0.5, size = 0.7) +
  geom_abline(color = "red") +
  labs(x = "Rubin et al N depo", y = "Zhu et al N depo") +
  theme_classic()
p_r4

p_r5 <- dt_raw %>% ggplot() + 
  geom_point(aes(x = n_depo_usa_coef, y = n_depo_zhu_coef), alpha = 0.5, size = 0.7) +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red") +
  labs(x = "USA N depo trend", y = "Zhu et al N depo trend") +
  theme_classic()
p_r5

p1 <- gridExtra::grid.arrange(p_mean, p_trend, ncol = 2)
p2 <- gridExtra::grid.arrange(p_r1, p_r2, p_r3, p_r4, p_r5, ncol = 3)

p_c <- gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(1, 2))
ggsave(plot = p_c, "builds/plots/supplement/usa_n_depo_maps_and_xy.png", dpi = 600, height = 8, width = 10)


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
    mean_n_depo_usa = scale(mean_n_depo_usa),
    n_depo_usa_coef = scale(n_depo_usa_coef), 
    mean_n_depo_zhu = scale(mean_n_depo_zhu),
    n_depo_zhu_coef = scale(n_depo_zhu_coef)) 

## EVI ------------------------------------------------------------------------------
####### calculate mean_evi trend #######

mean_evi_cols <- grep("mean_evi_", names(dt), value = T)

#subset to complete cases
dt_mean_evi <- dt %>%
  dplyr::select(all_of(mean_evi_cols), functional_biome, X, Y, lon, lat, unique_id, 
                mean_n_depo_usa, n_depo_usa_coef, mean_n_depo_zhu, n_depo_zhu_coef,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.)) %>% 
  as.data.frame()


#get y matrix
y_mean_evi <- as.matrix(dt_mean_evi[, mean_evi_cols])

#get coordinate matrix 
coords_mean_evi <- as.matrix(dt_mean_evi[, c("lon", "lat")])
any(duplicated(coords_mean_evi))

#fit autoregression, accounting for temporal autocorrelation 
ar_mean_evi <- fitAR_map(Y = y_mean_evi, coords = coords_mean_evi)

#extract coefficients and p values 
dt_mean_evi$mean_evi_coef <- ar_mean_evi$coefficients[, "t"]
dt_mean_evi$abs_mean_evi_coef <- abs(ar_mean_evi$coefficients[, "t"])
dt_mean_evi$mean_evi_p_value <- ar_mean_evi$pvals[, 2]


### estimate optimal r parameter (range of spatial autocorrelation)
corfit_mean_evi <- fitCor(resids = residuals(ar_mean_evi), coords = coords_mean_evi, covar_FUN = "covar_exp", 
                          start = list(range = 0.1), fit.n = 15000)
(range_opt_mean_evi = corfit_mean_evi$spcor) #0.07182121  



############ test hypotheses ############ 

#partition the data 
set.seed(161)
pm_mean_evi <- sample_partitions(npix = nrow(dt_mean_evi), partsize = 1500, npart = NA)
dim(pm_mean_evi) 


## USA N MEAN -------------

set.seed(161)
gls_mean_evi_1 <- fitGLS_partition(mean_evi_coef ~ 1 +
                                   mean_n_depo_usa +
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
                                 coord.names = c("lon", "lat"), 
                                 debug = F)
gls_mean_evi_1 

dt_est_mean_evi_1 <- extract_gls_estimates(gls_mean_evi_1, part = TRUE)

range(dt_mean_evi$lat)

p_mean_evi_1 <- dt_est_mean_evi_1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "USA mean N depo", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_mean_evi_1


## USA N TREND -------------

set.seed(161)
gls_mean_evi_2 <- fitGLS_partition(mean_evi_coef ~ 1 +
                                     n_depo_usa_coef +
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
gls_mean_evi_2 

dt_est_mean_evi_2 <- extract_gls_estimates(gls_mean_evi_2, part = TRUE)

p_mean_evi_2 <- dt_est_mean_evi_2 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "USA N depo trend", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_mean_evi_2

## USA N TREND & MEAN -------------

set.seed(161)
gls_mean_evi_3 <- fitGLS_partition(mean_evi_coef ~ 1 +
                                     n_depo_usa_coef +
                                     mean_n_depo_usa +
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
gls_mean_evi_3 

dt_est_mean_evi_3 <- extract_gls_estimates(gls_mean_evi_3, part = TRUE)

p_mean_evi_3 <- dt_est_mean_evi_3 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "USA mean N depo & trend", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_mean_evi_3

## USA OG N DEPO -------------

set.seed(161)
gls_mean_evi_4 <- fitGLS_partition(mean_evi_coef ~ 1 +
                                     nitrogen_depo +
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
gls_mean_evi_4

dt_est_mean_evi_4 <- extract_gls_estimates(gls_mean_evi_4, part = TRUE)

p_mean_evi_4 <- dt_est_mean_evi_4 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_mean_evi_4

### Zhu et al mean
set.seed(161)
gls_mean_evi_5 <- fitGLS_partition(mean_evi_coef ~ 1 +
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
gls_mean_evi_5

dt_est_mean_evi_5 <- extract_gls_estimates(gls_mean_evi_5, part = TRUE)

p_mean_evi_5 <- dt_est_mean_evi_5 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_mean_evi_5

#### Zhu trend 


set.seed(161)
gls_mean_evi_6 <- fitGLS_partition(mean_evi_coef ~ 1 +
                                     n_depo_zhu_coef +
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
gls_mean_evi_6

dt_est_mean_evi_6 <- extract_gls_estimates(gls_mean_evi_6, part = TRUE)

p_mean_evi_6 <- dt_est_mean_evi_6 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_mean_evi_6

## Zhu mean and trend 


set.seed(161)
gls_mean_evi_7 <- fitGLS_partition(mean_evi_coef ~ 1 +
                                     mean_n_depo_zhu +
                                     n_depo_zhu_coef +
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
gls_mean_evi_7

dt_est_mean_evi_7 <- extract_gls_estimates(gls_mean_evi_7, part = TRUE)

p_mean_evi_7 <- dt_est_mean_evi_7 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_mean_evi_7



p_evi_diff_n <- grid.arrange(p_mean_evi_1, p_mean_evi_2, p_mean_evi_3, p_mean_evi_4, 
                             p_mean_evi_5, p_mean_evi_6, p_mean_evi_7, ncol = 3)
ggsave(plot = p_evi_diff_n, "builds/plots/supplement/mean_evi_vs_various_n_depo.png", dpi = 600)


## BURNED AREA ------------------------------------------------------------------------------

####### calculate burned_area trend #######

burned_area_cols <- grep("burned_area_", names(dt), value = T)

#subset to complete cases
dt_burned_area <- dt %>%
  dplyr::select(all_of(burned_area_cols), functional_biome, X, Y, lon, lat, unique_id, 
                mean_n_depo_usa, n_depo_usa_coef, mean_n_depo_zhu, n_depo_zhu_coef,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.)) %>% as.data.frame()

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


## USA N MEAN -------------


## USA N MEAN -------------

set.seed(161)
gls_burned_area_1 <- fitGLS_partition(burned_area_coef ~ 1 +
                                     mean_n_depo_usa +
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
                                   coord.names = c("lon", "lat"), 
                                   debug = F)
gls_burned_area_1 

dt_est_burned_area_1 <- extract_gls_estimates(gls_burned_area_1, part = TRUE)

range(dt_burned_area$lat)

p_burned_area_1 <- dt_est_burned_area_1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "USA mean N depo", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_burned_area_1


## USA N TREND -------------

set.seed(161)
gls_burned_area_2 <- fitGLS_partition(burned_area_coef ~ 1 +
                                     n_depo_usa_coef +
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
gls_burned_area_2 

dt_est_burned_area_2 <- extract_gls_estimates(gls_burned_area_2, part = TRUE)

p_burned_area_2 <- dt_est_burned_area_2 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "USA N depo trend", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_burned_area_2

## USA N TREND & MEAN -------------

set.seed(161)
gls_burned_area_3 <- fitGLS_partition(burned_area_coef ~ 1 +
                                     n_depo_usa_coef +
                                     mean_n_depo_usa +
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
gls_burned_area_3 

dt_est_burned_area_3 <- extract_gls_estimates(gls_burned_area_3, part = TRUE)

p_burned_area_3 <- dt_est_burned_area_3 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "USA mean N depo & trend", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_burned_area_3

## USA OG N DEPO -------------

set.seed(161)
gls_burned_area_4 <- fitGLS_partition(burned_area_coef ~ 1 +
                                     nitrogen_depo +
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
gls_burned_area_4

dt_est_burned_area_4 <- extract_gls_estimates(gls_burned_area_4, part = TRUE)

p_burned_area_4 <- dt_est_burned_area_4 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_burned_area_4

### Zhu et al mean
set.seed(161)
gls_burned_area_5 <- fitGLS_partition(burned_area_coef ~ 1 +
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
gls_burned_area_5

dt_est_burned_area_5 <- extract_gls_estimates(gls_burned_area_5, part = TRUE)

p_burned_area_5 <- dt_est_burned_area_5 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_burned_area_5

#### Zhu trend 


set.seed(161)
gls_burned_area_6 <- fitGLS_partition(burned_area_coef ~ 1 +
                                     n_depo_zhu_coef +
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
gls_burned_area_6

dt_est_burned_area_6 <- extract_gls_estimates(gls_burned_area_6, part = TRUE)

p_burned_area_6 <- dt_est_burned_area_6 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_burned_area_6

## Zhu mean and trend 


set.seed(161)
gls_burned_area_7 <- fitGLS_partition(burned_area_coef ~ 1 +
                                     mean_n_depo_zhu +
                                     n_depo_zhu_coef +
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
gls_burned_area_7

dt_est_burned_area_7 <- extract_gls_estimates(gls_burned_area_7, part = TRUE)

p_burned_area_7 <- dt_est_burned_area_7 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_burned_area_7



p_evi_diff_n <- grid.arrange(p_burned_area_1, p_burned_area_2, p_burned_area_3, p_burned_area_4, 
                             p_burned_area_5, p_burned_area_6, p_burned_area_7, ncol = 3)
ggsave(plot = p_evi_diff_n, "builds/plots/supplement/burned_area_vs_various_n_depo.png", dpi = 600)
## GREEN-UP ------------------------------------------------------
####### calculate mean_evi trend #######

greenup_cols <- grep("greenup_", names(dt), value = T)

#subset to complete cases
dt_greenup <- dt %>%
  filter(!grepl("N", functional_biome)) %>% 
  dplyr::select(all_of(greenup_cols), functional_biome, X, Y, lon, lat, unique_id, 
                mean_n_depo_usa, n_depo_usa_coef, mean_n_depo_zhu, n_depo_zhu_coef,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.)) %>% as.data.frame()

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



## USA N MEAN -------------

set.seed(161)
gls_greenup_1 <- fitGLS_partition(greenup_coef ~ 1 +
                                     mean_n_depo_usa +
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
                                   coord.names = c("lon", "lat"), 
                                   debug = F)
gls_greenup_1 

dt_est_greenup_1 <- extract_gls_estimates(gls_greenup_1, part = TRUE)

range(dt_greenup$lat)

p_greenup_1 <- dt_est_greenup_1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
  labs(title = "USA mean N depo", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_greenup_1


## USA N TREND -------------

set.seed(161)
gls_greenup_2 <- fitGLS_partition(greenup_coef ~ 1 +
                                     n_depo_usa_coef +
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
gls_greenup_2 

dt_est_greenup_2 <- extract_gls_estimates(gls_greenup_2, part = TRUE)

p_greenup_2 <- dt_est_greenup_2 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
  labs(title = "USA N depo trend", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_greenup_2

## USA N TREND & MEAN -------------

set.seed(161)
gls_greenup_3 <- fitGLS_partition(greenup_coef ~ 1 +
                                     n_depo_usa_coef +
                                     mean_n_depo_usa +
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
gls_greenup_3 

dt_est_greenup_3 <- extract_gls_estimates(gls_greenup_3, part = TRUE)

p_greenup_3 <- dt_est_greenup_3 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
  labs(title = "USA mean N depo & trend", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_greenup_3

## USA OG N DEPO -------------

set.seed(161)
gls_greenup_4 <- fitGLS_partition(greenup_coef ~ 1 +
                                     nitrogen_depo +
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
gls_greenup_4

dt_est_greenup_4 <- extract_gls_estimates(gls_greenup_4, part = TRUE)

p_greenup_4 <- dt_est_greenup_4 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_greenup_4

### Zhu et al mean
set.seed(161)
gls_greenup_5 <- fitGLS_partition(greenup_coef ~ 1 +
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
gls_greenup_5

dt_est_greenup_5 <- extract_gls_estimates(gls_greenup_5, part = TRUE)

p_greenup_5 <- dt_est_greenup_5 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_greenup_5

#### Zhu trend 


set.seed(161)
gls_greenup_6 <- fitGLS_partition(greenup_coef ~ 1 +
                                     n_depo_zhu_coef +
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
gls_greenup_6

dt_est_greenup_6 <- extract_gls_estimates(gls_greenup_6, part = TRUE)

p_greenup_6 <- dt_est_greenup_6 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_greenup_6

## Zhu mean and trend 


set.seed(161)
gls_greenup_7 <- fitGLS_partition(greenup_coef ~ 1 +
                                     mean_n_depo_zhu +
                                     n_depo_zhu_coef +
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
gls_greenup_7

dt_est_greenup_7 <- extract_gls_estimates(gls_greenup_7, part = TRUE)

p_greenup_7 <- dt_est_greenup_7 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
  labs(title = "Rubin et al. N depo", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_greenup_7



p_evi_diff_n <- grid.arrange(p_greenup_1, p_greenup_2, p_greenup_3, p_greenup_4, 
                             p_greenup_5, p_greenup_6, p_greenup_7, ncol = 3)
ggsave(plot = p_evi_diff_n, "builds/plots/supplement/greenup_vs_various_n_depo.png", dpi = 600)

########## STOPPED HERE FOR NOW ############


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
    grepl("mean_n_depo_usa", term) ~ "Nitrogen Deposition (Zhu)",
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
  scale_color_manual(values = c("Sig. Negative" = "#2064AE",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#995215"
  )) +
  labs(title = "Nitrogen Deposition by Zhu et al., 2025", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
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

