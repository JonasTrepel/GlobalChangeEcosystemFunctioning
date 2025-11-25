library(tidyverse)
library(data.table)
library(RColorBrewer)
library(scico)
library(rnaturalearth)
library(sf)
library(gridExtra)
library(tidylog)
library(RColorBrewer)

source("R/functions/extract_gls_estimates.R")


mean_evi_trend <- fread("data/processedData/dataFragments/grid_mean_evi_trends.csv") 
burned_area_trend <- fread("data/processedData/dataFragments/grid_burned_area_trends.csv")
greenup_trend <- fread("data/processedData/dataFragments/grid_greenup_trends.csv")


climate_trends <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") %>% 
  mutate(mean_burned_area = rowMeans(select(., contains("burned_area")), na.rm = TRUE)) %>% 
  dplyr::select(mean_burned_area, unique_id)


dt_raw <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") 


dt <- dt_raw %>% 
  mutate(
    pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
    nitrogen_depo = scale(nitrogen_depo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(human_modification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA)) %>% 
  left_join(climate_trends) %>% 
  left_join(mean_evi_trend) %>% 
  left_join(burned_area_trend) %>% 
  left_join(greenup_trend) %>% 
  mutate(burned_area_coef = scale(burned_area_coef))

mean_evi_cols <- grep("mean_evi_", names(dt), value = T)

#subset to complete cases
dt_mean_evi <- dt %>%
  dplyr::select(all_of(mean_evi_cols), functional_biome, lon, lat, unique_id, 
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad, burned_area_coef) %>% 
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
corfit <- fitCor(resids = residuals(ar_mean_evi), coords = coords_mean_evi, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = nrow(dt_mean_evi))
(range_opt_mean_evi = corfit$spcor)

#range_opt_mean_evi <- 0.01

############ test hypotheses ############ 

#partition the data 
set.seed(161)
pm <- sample_partitions(npix = nrow(dt_mean_evi), partsize = 1500, npart = NA)
dim(pm)

set.seed(161)
gls_h2 <- fitGLS_partition(mean_evi_coef ~ 1 +
                             nitrogen_depo +
                             mat_coef + 
                             max_temp_coef + 
                             map_coef + 
                             human_modification + 
                             burned_area_coef,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_mean_evi),
                           data = dt_mean_evi,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("lon", "lat"))
gls_h2 

dt_est_h2 <- extract_gls_estimates(gls_h2, part = TRUE)

p_h2 <- dt_est_h2 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H2: mean_evi change ~\nglobal change", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h2

######

dt %>% 
  filter(mean_burned_area > 0) %>% 
  ggplot() +
  geom_point(aes(x = burned_area_coef, y = mean_evi_coef))
