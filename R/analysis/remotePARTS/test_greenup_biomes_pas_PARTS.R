#wtf is going on with greenup

source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  group_by(FunctionalBiome) %>% 
  mutate(n_per_functional_biome = n()) %>% 
  ungroup() %>% 
  group_by(Biome) %>% 
  mutate(n_per_biome = n()) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  # filter(n_per_functional_biome > 1000) %>% 
  rename(functional_biome = FunctionalBiome) %>% 
  mutate(
    productivity = case_when(
      grepl("L", functional_biome) ~ "low", 
      grepl("M", functional_biome) ~ "medium", 
      grepl("H", functional_biome) ~ "high"
    ), 
    ndvi_min = case_when(
      grepl("C", functional_biome) ~ "cold", 
      grepl("D", functional_biome) ~ "dry", 
      grepl("B", functional_biome) ~ "cold_and_dry", 
      grepl("N", functional_biome) ~ "non_seasonal"
    ), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"
    ),
    nitrogen_depo = scale(NitrogenDepo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(HumanModification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = scale(log(PaAge + 0.0001))) %>% 
  filter(n_per_biome > 100 &
           !is.na(ClimaticRegion) &
           n_per_functional_biome > 100 & 
           !ClimaticRegion == "") %>% 
  rename(biome = Biome, 
         climatic_region = ClimaticRegion) 

dt <- dt_raw
table(dt$climatic_region)
#######################################

####### calculate greenup trend #######

greenup_cols <- grep("greenup_", names(dt), value = T)

#subset to complete cases
dt_greenup <- dt %>%
  dplyr::select(all_of(greenup_cols),  X, Y, unique_id, productivity, ndvi_min,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification,
                super_biome, climatic_region, biome, functional_biome,
                og_layer) %>% 
  filter(complete.cases(.))

#get y matrix
y_greenup <- as.matrix(dt_greenup[, greenup_cols])

#get coordinate matrix 
coords_greenup <- as.matrix(dt_greenup[, c("X", "Y")])

#fit autoregression, accounting for temporal autocorrelation 
ar_greenup <- fitAR_map(Y = y_greenup, coords = coords_greenup)

#extract coefficients and p values 
dt_greenup$greenup_coef <- ar_greenup$coefficients[, "t"]
dt_greenup$abs_greenup_coef <- abs(ar_greenup$coefficients[, "t"])
dt_greenup$greenup_p_value <- ar_greenup$pvals[, 2]

#fwrite(dt_greenup %>% dplyr::select(
#  unique_id, X, Y, og_layer, greenup_coef, abs_greenup_coef, greenup_p_value), "data/processedData/dataFragments/pa_greenup_trends.csv")

# get distance matrix 
#d_greenup <- distm_scaled(coords_greenup)

### estimate optimal r parameter (range of spatial autocorrelation)
corfit <- fitCor(resids = residuals(ar_greenup), coords = coords_greenup, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = nrow(dt_greenup))
(range_opt_greenup = corfit$spcor)

#use r to calculate optimal v parameter 
#v_opt_greenup <- covar_exp(d_greenup, range_opt_greenup)

############ test hypotheses ############ 

#partition the data 
set.seed(161)
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
                           ncores = 5,
                           progressbar = TRUE, 
                           parallel = T, 
                           coord.names = c("X", "Y")
)
gls_h1 # yes. Est: -0.05098253; SE: 0.1024441; pval.t: 0.6187328

dt_est_h1 <- extract_gls_estimates(gls_h1, part = TRUE)

# Olson Biomes ---------------------------------
set.seed(161)
gls_biome <- fitGLS_partition(greenup_coef ~ 0 +
                             biome,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_greenup),
                           data = dt_greenup,
                           nugget = NA,
                           ncores = 5,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("X", "Y"))
gls_biome 

dt_est_biome <- extract_gls_estimates(gls_biome, part = TRUE) %>%
  mutate(term = gsub("biome", "", term))

p_biome <- dt_est_biome %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
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
                                    ncores = 5,
                                    progressbar = TRUE, 
                                    parallel = TRUE, 
                                    coord.names = c("X", "Y"))
gls_functional_biome 

dt_est_functional_biome <- extract_gls_estimates(gls_functional_biome, part = TRUE) %>%
  mutate(term = gsub("functional_biome", "", term))

p_functional_biome <- dt_est_functional_biome %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
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
                              ncores = 5,
                              progressbar = TRUE, 
                              parallel = TRUE, 
                              coord.names = c("X", "Y"))
gls_super_biome 

dt_est_super_biome <- extract_gls_estimates(gls_super_biome, part = TRUE) %>%
  mutate(term = gsub("super_biome", "", term))

p_super_biome <- dt_est_super_biome %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
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
                                    ncores = 5,
                                    progressbar = TRUE, 
                                    parallel = TRUE, 
                                    coord.names = c("X", "Y"))
gls_climatic_region 

dt_est_climatic_region <- extract_gls_estimates(gls_climatic_region, part = TRUE) %>%
  mutate(term = gsub("climatic_region", "", term))

p_climatic_region <- dt_est_climatic_region %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "darkturquoise", "non-significant" = "grey")) +
  labs(title = "Greenup Change ~\nClimatic Region", subtitle = paste0("n = ", nrow(dt_greenup)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_climatic_region


p_gr_b <- gridExtra::grid.arrange(p_biome, p_functional_biome, p_super_biome, p_climatic_region,
                                  widths = c(2, 0.9, 1, 1))
ggsave(plot = p_gr_b, "builds/plots/greenup_remotePARTS_different_biome_classifications.png", dpi = 600, 
       height = 3, width = 12)
