#### test burned_area hypothesis ####
source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  as.data.frame()

dt_raw <- dt_raw %>% 
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



dt <- dt_raw


#######################################

####### calculate burned_area trend #######

burned_area_cols <- grep("burned_area_", names(dt), value = T)

#subset to complete cases
dt_burned_area <- dt %>%
  dplyr::select(all_of(burned_area_cols), functional_biome, lon, lat, unique_id, lon_pa, lat_pa, unique_pa_id,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.)) %>% 
  distinct(lon, lat, .keep_all = T)

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

fwrite(dt_burned_area %>% dplyr::select(
  unique_id, burned_area_coef, abs_burned_area_coef, burned_area_p_value), "data/processedData/dataFragments/pas_burned_area_trends.csv")

# get distance matrix 
#d_burned_area <- distm_scaled(coords_burned_area)

### estimate optimal r parameter (range of spatial autocorrelation)
corfit <- fitCor(resids = residuals(ar_burned_area), coords = coords_burned_area, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = 15000)
(range_opt_burned_area = corfit$spcor) ##0.001099433


##### by protected area 

get_mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

dt_sum <- dt_burned_area %>%
  #mutate(unique_pa_id = paste0(lat_pa, lon_pa)) %>% 
  group_by(unique_pa_id) %>% 
  summarize(mat_coef = mean(mat_coef, na.rm = T),
            map_coef = mean(map_coef, na.rm = T),
            max_temp_coef = mean(max_temp_coef, na.rm = T),
            human_modification = mean(human_modification, na.rm = T),
            nitrogen_depo = mean(nitrogen_depo, na.rm = T),
            protection_cat_broad = unique(protection_cat_broad), 
            lon_pa = unique(lon_pa), 
            lat_pa = unique(lat_pa), 
            burned_area_coef = mean(burned_area_coef, na.rm = T),
            abs_burned_area_coef = mean(abs_burned_area_coef, na.rm = T)) 




############ test hypotheses ############ 

#partition the data 
set.seed(161)
pm <- sample_partitions(npix = nrow(dt_sum), partsize = 1000, npart = NA)
dim(pm)

## 
gls_h4.0 <- fitGLS_partition(burned_area_coef ~ 1,
                             partmat = pm,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = range_opt_burned_area),
                             data = dt_sum,
                             nugget = NA,
                             ncores = 10,
                             progressbar = TRUE, 
                             parallel = T, 
                             coord.names = c("lon_pa", "lat_pa")
)
gls_h4.0 # yes. Est: 5.563445 ; SE: 0.4952674 ; pval.t: 2.871102e-29

dt_est_h4.0 <- extract_gls_estimates(gls_h4.0, part = TRUE)


## Hypothesis 4 - different absolute trend in and outside of PAs ----------------------
set.seed(161)
gls_h4 <- fitGLS_partition(burned_area_coef ~ 0 +
                             protection_cat_broad,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_burned_area),
                           data = dt_sum,
                           nugget = NA,
                           ncores = 10,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("lon_pa", "lat_pa"))
gls_h4 # 


dt_est_h4 <- extract_gls_estimates(gls_h4, part = TRUE)
p_h4 <- dt_est_h4 %>% 
  ggplot() + 
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H4: burned_area change ~\nprotection", subtitle = paste0("n = ", nrow(dt_sum)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h4


## Hypothesis 4.1 - different absolute trend in and outside of PAs ------------------
set.seed(161)
gls_h4.1 <- fitGLS_partition(abs_burned_area_coef ~ 0 +
                               protection_cat_broad,
                             partmat = pm,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = range_opt_burned_area),
                             data = dt_sum,
                             nugget = NA,
                             ncores = 10,
                             progressbar = TRUE, 
                             parallel = TRUE, 
                             coord.names = c("lon_pa", "lat_pa"))
gls_h4.1 # 


dt_est_h4.1 <- extract_gls_estimates(gls_h4.1, part = TRUE)

p_h4.1 <- dt_est_h4.1 %>% 
  ggplot() + 
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H4.1: abs burned_area change ~\nprotection", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h4.1


## Hypothesis 5 - PA characteristics ----------------------------

# recalculate trends only on PAs 
burned_area_cols <- grep("burned_area_", names(dt), value = T)

#subset to complete cases
dt_pa <- dt %>% filter(protection_cat_broad == "strictly_protected") %>%
  dplyr::select(all_of(burned_area_cols),functional_biome, lon, lat, unique_id, lon_pa, lat_pa, unique_pa_id,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad, area_km2_log, pa_age_log) %>% 
  filter(complete.cases(.))

#get y matrix
y_pa <- as.matrix(dt_pa[, burned_area_cols])

#get coordinate matrix 
coords_pa <- as.matrix(dt_pa[, c("lon", "lat")])

#fit autoregression, accounting for temporal autocorrelation 
ar_pa <- fitAR_map(Y = y_pa, coords = coords_pa)

#extract coefficients and p values 
dt_pa$burned_area_coef <- ar_pa$coefficients[, "t"]
dt_pa$abs_burned_area_coef <- abs(ar_pa$coefficients[, "t"])
dt_pa$burned_area_p_value <- ar_pa$pvals[, 2]

### estimate optimal r parameter (range of spatial autocorrelation)
corfit_pa <- fitCor(resids = residuals(ar_pa), coords = coords_pa, covar_FUN = "covar_exp", 
                    start = list(range = 0.1), fit.n = 15000)
(range_opt_pa = corfit_pa$spcor)  

#summarize PA 

dt_pa_sum <- dt_pa %>%
  group_by(unique_pa_id) %>% 
  summarize(mat_coef = mean(mat_coef, na.rm = T),
            map_coef = mean(map_coef, na.rm = T),
            max_temp_coef = mean(max_temp_coef, na.rm = T),
            area_km2_log = mean(area_km2_log, na.rm = T),
            pa_age_log = mean(pa_age_log, na.rm = T),
            human_modification = mean(human_modification, na.rm = T),
            nitrogen_depo = mean(nitrogen_depo, na.rm = T),
            protection_cat_broad = unique(protection_cat_broad), 
            lon_pa = unique(lon_pa), 
            lat_pa = unique(lat_pa), 
            burned_area_coef = mean(burned_area_coef, na.rm = T),
            abs_burned_area_coef = mean(abs_burned_area_coef, na.rm = T)) 



coords_sum_pa <- as.matrix(dt_pa_sum[, c("lon_pa", "lat_pa")])

d_m_pa <- distm_scaled(coords_sum_pa)
v_opt_pa <- covar_exp(d_m_pa, range_opt_pa)
#partition the data 

pm_pa <- sample_partitions(npix = nrow(dt_pa_sum), partsize = 1000, npart = NA)
dim(pm_pa)

gls_h5.0 <- fitGLS(burned_area_coef ~ 1,
                   data = dt_pa_sum,
                   nugget = NA,
                   no.F = T,
                   V = v_opt_pa)
gls_h5.0 

dt_est_h5.0 <- extract_gls_estimates(gls_h5.0, part = FALSE)

# change --- 
gls_h5 <- fitGLS(burned_area_coef ~ 1 +
                   area_km2_log + 
                   pa_age_log,
                 data = dt_pa_sum,
                 nugget = NA,
                 no.F = T,
                 V = v_opt_pa)
gls_h5 

dt_est_h5 <- extract_gls_estimates(gls_h5, part = FALSE)

p_h5 <- dt_est_h5 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H5: burned_area change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa_sum)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h5

## Hypothesis 5.1 - absolute change PA characteristics ----------------------------
gls_h5.1 <- fitGLS(abs_burned_area_coef ~ 1 +
                     area_km2_log + 
                     pa_age_log,
                   data = dt_pa_sum,
                   nugget = NA,
                   no.F = T,
                   V = v_opt_pa)
gls_h5.1 

dt_est_h5.1 <- extract_gls_estimates(gls_h5.1, part = FALSE)

p_h5.1 <- dt_est_h5.1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H5.1: abs burned_area change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa_sum)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h5.1

#### Median instead of mean -----------

### PA characteristics ---------

dt_pa_sum_median <- dt_pa %>%
  group_by(unique_pa_id) %>% 
  summarize(mat_coef = median(mat_coef, na.rm = T),
            map_coef = median(map_coef, na.rm = T),
            max_temp_coef = median(max_temp_coef, na.rm = T),
            area_km2_log = median(area_km2_log, na.rm = T),
            pa_age_log = median(pa_age_log, na.rm = T),
            human_modification = median(human_modification, na.rm = T),
            nitrogen_depo = median(nitrogen_depo, na.rm = T),
            protection_cat_broad = unique(protection_cat_broad), 
            lon_pa = unique(lon_pa), 
            lat_pa = unique(lat_pa), 
            burned_area_coef = median(burned_area_coef, na.rm = T),
            abs_burned_area_coef = median(abs_burned_area_coef, na.rm = T)) 



coords_sum_pa <- as.matrix(dt_pa_sum_median[, c("lon_pa", "lat_pa")])

d_m_pa <- distm_scaled(coords_sum_pa)
v_opt_pa <- covar_exp(d_m_pa, range_opt_pa)

gls_h6.0 <- fitGLS(burned_area_coef ~ 1,
                   data = dt_pa_sum_median,
                   nugget = NA,
                   no.F = T,
                   V = v_opt_pa)
gls_h6.0 

dt_est_h6.0 <- extract_gls_estimates(gls_h6.0, part = FALSE)

# change --- 
gls_h6 <- fitGLS(burned_area_coef ~ 1 +
                   area_km2_log + 
                   pa_age_log,
                 data = dt_pa_sum_median,
                 nugget = NA,
                 no.F = T,
                 V = v_opt_pa)
gls_h6 

dt_est_h6 <- extract_gls_estimates(gls_h6, part = FALSE)

p_h6 <- dt_est_h6 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "h6: burned_area change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa_sum_median)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h6

## Hypothesis 5.1 - absolute change PA characteristics ----------------------------
gls_h6.1 <- fitGLS(abs_burned_area_coef ~ 1 +
                     area_km2_log + 
                     pa_age_log,
                   data = dt_pa_sum_median,
                   nugget = NA,
                   no.F = T,
                   V = v_opt_pa)
gls_h6.1 

dt_est_h6.1 <- extract_gls_estimates(gls_h6.1, part = FALSE)

p_h6.1 <- dt_est_h6.1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "h6.1: abs burned_area change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa_sum_median)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h6.1

### Protection --------------------

dt_sum_median <- dt_burned_area %>%
  # mutate(unique_pa_id = paste0(lat_pa, lon_pa)) %>% 
  group_by(unique_pa_id) %>% 
  summarize(mat_coef = median(mat_coef, na.rm = T),
            map_coef = median(map_coef, na.rm = T),
            max_temp_coef = median(max_temp_coef, na.rm = T),
            human_modification = median(human_modification, na.rm = T),
            nitrogen_depo = median(nitrogen_depo, na.rm = T),
            protection_cat_broad = unique(protection_cat_broad), 
            lon_pa = unique(lon_pa), 
            lat_pa = unique(lat_pa), 
            burned_area_coef = median(burned_area_coef, na.rm = T),
            abs_burned_area_coef = median(abs_burned_area_coef, na.rm = T)) 


## 
gls_h7.0 <- fitGLS_partition(burned_area_coef ~ 1,
                             partmat = pm,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = range_opt_burned_area),
                             data = dt_sum_median,
                             nugget = NA,
                             ncores = 10,
                             progressbar = TRUE, 
                             parallel = TRUE, 
                             coord.names = c("lon_pa", "lat_pa"))
gls_h7.0 # 

dt_est_h7.0 <- extract_gls_estimates(gls_h7.0, part = TRUE)


## Hypothesis 4 - different absolute trend in and outside of PAs ----------------------
set.seed(161)
gls_h7 <- fitGLS_partition(burned_area_coef ~ 0 +
                             protection_cat_broad,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_burned_area),
                           data = dt_sum_median,
                           nugget = NA,
                           ncores = 10,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("lon_pa", "lat_pa"))
gls_h7 # 


dt_est_h7 <- extract_gls_estimates(gls_h7, part = TRUE)
p_h7 <- dt_est_h7 %>% 
  ggplot() + 
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "h7: burned_area change ~\nprotection", subtitle = paste0("n = ", nrow(dt_sum_median)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h7


## Hypothesis 4.1 - different absolute trend in and outside of PAs ------------------
set.seed(161)
gls_h7.1 <- fitGLS_partition(abs_burned_area_coef ~ 0 +
                               protection_cat_broad,
                             partmat = pm,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = range_opt_burned_area),
                             data = dt_sum_median,
                             nugget = NA,
                             ncores = 10,
                             progressbar = TRUE, 
                             parallel = TRUE, 
                             coord.names = c("lon_pa", "lat_pa"))
gls_h7.1 # 


dt_est_h7.1 <- extract_gls_estimates(gls_h7.1, part = TRUE)

p_h7.1 <- dt_est_h7.1 %>% 
  ggplot() + 
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "h7.1: abs burned_area change ~\nprotection", subtitle = paste0("n = ", nrow(dt_sum_median)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h7.1


library(gridExtra)
p_burned_area_a <- grid.arrange(p_h4, p_h4.1, p_h5, p_h5.1, p_h6, p_h6.1, p_h7, p_h7.1, ncol = 2)
ggsave(plot = p_burned_area_a, "builds/plots/burned_area_remotePARTS_pas.png", dpi = 600, height = 3, width = 12)



dt_est <- rbind(
  dt_est_h4.0 %>% mutate(model = "H4.0", t.stat = NA, response = "burned_area"),
  dt_est_h4 %>% mutate(model = "H4", t.stat = NA, response = "burned_area"), 
  dt_est_h4.1 %>% mutate(model = "H4.1", t.stat = NA, response = "burned_area"), 
  dt_est_h5 %>% mutate(model = "H5.0", t.stat = NA, response = "burned_area"), 
  dt_est_h5 %>% mutate(model = "H5", t.stat = NA, response = "burned_area"), 
  dt_est_h5.1 %>% mutate(model = "H5.1", t.stat = NA, response = "burned_area"),
  dt_est_h6 %>% mutate(model = "H6.0", t.stat = NA, response = "burned_area"), 
  dt_est_h6 %>% mutate(model = "H6", t.stat = NA, response = "burned_area"), 
  dt_est_h6.1 %>% mutate(model = "H6.1", t.stat = NA, response = "burned_area"),
  dt_est_h7 %>% mutate(model = "H7.0", t.stat = NA, response = "burned_area"), 
  dt_est_h7 %>% mutate(model = "H7", t.stat = NA, response = "burned_area"), 
  dt_est_h7.1 %>% mutate(model = "H7.1", t.stat = NA, response = "burned_area")
)

fwrite(dt_est, "builds/model_estimates/burned_area_pas_estimates.csv")

