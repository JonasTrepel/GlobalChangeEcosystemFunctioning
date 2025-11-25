#wtf is going on with greenup

source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw_raw <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") %>% 
  as.data.frame()

dt_raw <- dt_raw_raw %>% 
  mutate(
    pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
    nitrogen_depo = scale(nitrogen_depo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(human_modification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA)) %>% 
    group_by(functional_biome) %>% 
    mutate(n_per_functional_biome = n()) %>% 
    ungroup() %>% 
    group_by(olson_biome) %>% 
    mutate(n_per_olson_biome = n()) %>% 
    ungroup()  %>% 
  filter(n_per_olson_biome > 2500 &
           !is.na(climatic_region) &
           n_per_functional_biome > 2500 & 
           !climatic_region == "" & !olson_biome == "") %>% 
  #group_by(olson_biome) %>% 
  #slice_sample(n = 1000) %>%
  ungroup()

table(dt_raw$functional_biome)
table(dt_raw$olson_biome)
table(dt_raw$climatic_region)
table(dt_raw$super_biome)

dt <- dt_raw
#######################################

####### calculate greenup trend #######

greenup_cols <- grep("greenup_", names(dt), value = T)

#subset to complete cases
dt_greenup <- dt %>%
  dplyr::select(all_of(greenup_cols),  lon, lat, unique_id,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification,
                super_biome, climatic_region, olson_biome, functional_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.))

table(dt_greenup$functional_biome)
table(dt_greenup$olson_biome)
table(dt_greenup$climatic_region)
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

#fwrite(dt_greenup %>% dplyr::select(
#  unique_id, lon, lat, og_layer, greenup_coef, abs_greenup_coef, greenup_p_value), "data/processedData/dataFragments/pa_greenup_trends.csv")

# get distance matrix 
#d_greenup <- distm_scaled(coords_greenup)

### estimate optimal r parameter (range of spatial autocorrelation)
corfit <- fitCor(resids = residuals(ar_greenup), coords = coords_greenup, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = 15000)
(range_opt_greenup = corfit$spcor)

#use r to calculate optimal v parameter 
#v_opt_greenup <- covar_exp(d_greenup, range_opt_greenup)

############ test hypotheses ############ 

#partition the data 
set.seed(1312)
pm <- sample_partitions(npix = nrow(dt_greenup), partsize = 1500, npart = NA)
dim(pm)

## Hypothesis 1: Ecosystem functioning is overall changing --------------------
set.seed(161)
gls_h1 <- fitGLS_partition(greenup_coef ~ 1,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_greenup),
                           data = dt_greenup,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = T, 
                           coord.names = c("lon", "lat")
)
gls_h1 # yes. Est: -0.05098253; SE: 0.1024441; pval.t: 0.6187328

dt_est_h1 <- extract_gls_estimates(gls_h1, part = TRUE)

# Olson Biomes ---------------------------------
set.seed(161)
gls_biome <- fitGLS_partition(greenup_coef ~ 0 +
                                olson_biome,
                              partmat = pm,
                              covar_FUN = "covar_exp",
                              covar.pars = list(range = range_opt_greenup),
                              data = dt_greenup,
                              nugget = NA,
                              ncores = 25,
                              progressbar = TRUE, 
                              parallel = FALSE, 
                              coord.names = c("lon", "lat"))
gls_biome 

dt_est_biome <- extract_gls_estimates(gls_biome, part = TRUE) %>%
  mutate(term = gsub("olson_biome", "", term), 
         term = gsub("_", " ", term))

p_biome <- dt_est_biome %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkcyan", "non-significant" = "grey")) +
  labs(title = "Greenup Change ~\nOlson Biomes", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_biome

# Functional Biomes  ---------------------------------
set.seed(161)
gls_functional_biome <- fitGLS_partition(greenup_coef ~ 0 +
                                           functional_biome,
                                         partmat = pm,
                                         covar_FUN = "covar_exp",
                                         covar.pars = list(range = range_opt_greenup),
                                         data = dt_greenup,
                                         nugget = NA,
                                         ncores = 25,
                                         progressbar = TRUE, 
                                         parallel = TRUE, 
                                         coord.names = c("lon", "lat"))
gls_functional_biome 

dt_est_functional_biome <- extract_gls_estimates(gls_functional_biome, part = TRUE) %>%
  mutate(term = gsub("functional_biome", "", term))

p_functional_biome <- dt_est_functional_biome %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkcyan", "non-significant" = "grey")) +
  labs(title = "Greenup Change ~\nFunctional Biome", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_functional_biome

# Super Biomes  ---------------------------------
set.seed(161)
gls_super_biome <- fitGLS_partition(greenup_coef ~ 0 +
                                      super_biome,
                                    partmat = pm,
                                    covar_FUN = "covar_exp",
                                    covar.pars = list(range = range_opt_greenup),
                                    data = dt_greenup,
                                    nugget = NA,
                                    ncores = 25,
                                    progressbar = TRUE, 
                                    parallel = TRUE, 
                                    coord.names = c("lon", "lat"))
gls_super_biome 

dt_est_super_biome <- extract_gls_estimates(gls_super_biome, part = TRUE) %>%
  mutate(term = gsub("super_biome", "", term), 
         term = case_when(
           term == "cold_short" ~ "Cold Limited\nShort Vegetation",
           term == "cold_tall" ~ "Cold Limited\nTall Vegetation",
           term == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
           term == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation"))

p_super_biome <- dt_est_super_biome %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkcyan", "non-significant" = "grey")) +
  labs(title = "Greenup Change ~\nSuper Biome", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_super_biome

# Climatic Region  ---------------------------------
set.seed(1312)
gls_climatic_region <- fitGLS_partition(greenup_coef ~ 0 +
                                          climatic_region,
                                        partmat = pm,
                                        covar_FUN = "covar_exp",
                                        covar.pars = list(range = range_opt_greenup),
                                        data = dt_greenup,
                                        nugget = NA,
                                        ncores = 25,
                                        progressbar = TRUE, 
                                        parallel = TRUE, 
                                        coord.names = c("lon", "lat"))
gls_climatic_region 

dt_est_climatic_region <- extract_gls_estimates(gls_climatic_region, part = TRUE) %>%
  mutate(term = gsub("climatic_region", "", term))

p_climatic_region <- dt_est_climatic_region %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkcyan", "non-significant" = "grey")) +
  labs(title = "Greenup Change ~\nClimatic Region", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_climatic_region

p_gr_b <- gridExtra::grid.arrange(p_biome, p_functional_biome, p_super_biome, p_climatic_region,
                                  widths = c(2, 0.9, 1, 1))
ggsave(plot = p_gr_b, "builds/plots/grid_greenup_remotePARTS_different_biome_classifications.png", dpi = 600, 
       height = 3, width = 12)

