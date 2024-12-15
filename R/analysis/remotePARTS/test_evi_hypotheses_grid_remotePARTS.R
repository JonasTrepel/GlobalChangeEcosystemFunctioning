#### test EVI hypothesis ####
source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") %>%
  mutate(PaAge = 2023-STATUS_YR, 
         PaAge = ifelse(PaAge > 200, NA, PaAge)) %>%
  group_by(FunctionalBiome) %>% 
  mutate(n_per_functional_biome = n()) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  filter(n_per_functional_biome > 1000) %>% 
  mutate(
    nitrogen_depo = scale(NitrogenDepo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(HumanModification), 
    area_km2_log = scale(log(area_km2)),
    pa_age_log = scale(log(PaAge + 0.0001))) %>% 
  rename(functional_biome = FunctionalBiomeShort)
dt <- dt_raw 

##### Subset for testing only #######

# dt_sub_pa <- dt_raw %>%
#   filter(protection_cat_broad == "protected_areas") %>%
#   filter(complete.cases(.)) %>%
#   sample_n(5000)
# 
# dt_sub_cont <- dt_raw %>%
#   filter(control_for %in% unique(dt_sub_pa$unique_id))
# 
# #dt <- dt_raw
# dt <- rbind(dt_sub_pa, dt_sub_cont) # %>% filter(complete.cases(.))
#######################################

####### calculate EVI trend #######

evi_cols <- grep("evi_", names(dt), value = T)

#subset to complete cases
dt_evi <- dt %>%
  dplyr::select(all_of(evi_cols), functional_biome, X, Y, unique_id,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, 
                protection_cat_broad) %>% 
  filter(complete.cases(.))  %>%
  distinct(X, Y, .keep_all = TRUE)

#get y matrix
y_evi <- as.matrix(dt_evi[, evi_cols])

#get coordinate matrix 
coords_evi <- as.matrix(dt_evi[, c("X", "Y")])

#fit autoregression, accounting for temporal autocorrelation 
ar_evi <- fitAR_map(Y = y_evi, coords = coords_evi)

#extract coefficients and p values 
dt_evi$evi_coef <- ar_evi$coefficients[, "t"]
dt_evi$abs_evi_coef <- abs(ar_evi$coefficients[, "t"])
dt_evi$evi_p_value <- ar_evi$pvals[, 2]

fwrite(dt_evi %>% dplyr::select(
  unique_id, X, Y, protection_cat_broad, evi_coef, abs_evi_coef, evi_p_value
), "data/processedData/dataFragments/grid_evi_trends.csv")

# get distance matrix 
#d_evi <- distm_scaled(coords_evi)

### estimate optimal r parameter (range of spatial autocorrelation)


corfit <- fitCor(resids = residuals(ar_evi), coords = coords_evi, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = 10000)
(range_opt_evi = corfit$spcor)

#use r to calculate optimal v parameter 
#v_opt_evi <- covar_exp(d_evi, range_opt_evi)

############ test hypotheses ############ 

#partition the data 

pm <- sample_partitions(npix = nrow(dt_evi), partsize = 1500, npart = NA)
dim(pm)

## Hypothesis 1: Ecosystem functioning is overall changing --------------------
gls_h1 <- fitGLS_partition(evi_coef ~ 1,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_evi),
                           data = dt_evi,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("X", "Y")
)
gls_h1 # yes. 

## Hypothesis 2: Change depends on climate change, N deposition, human modification -------------

set.seed(161)
gls_h2 <- fitGLS_partition(evi_coef ~ 1 +
                             nitrogen_depo +
                             mat_coef + 
                             max_temp_coef + 
                             map_coef + 
                             human_modification,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_evi),
                           data = dt_evi,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("X", "Y"))
gls_h2 #

dt_est_h2 <- extract_gls_estimates(gls_h2, part = TRUE)
p_h2 <- dt_est_h2 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H2: evi change ~\nglobal change", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none")
p_h2

## Hypothesis 3 - different trends in different biomes ----------------------
set.seed(161)
gls_h3 <- fitGLS_partition(evi_coef ~ 0 +
                             functional_biome,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_evi),
                           data = dt_evi,
                           nugget = NA,
                           ncores = 4,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("X", "Y"))
gls_h3 # 


dt_est_h3 <- extract_gls_estimates(gls_h3, part = TRUE)
p_h3 <- dt_est_h3 %>% 
  ggplot() + 
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H3: evi change ~\nfunctional biome", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none")
p_h3


## Hypothesis 4 - different trend in and outside of PAs ------------------
set.seed(161)
gls_h4 <- fitGLS_partition(evi_coef ~ 0 +
                             protection_cat_broad,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_evi),
                           data = dt_evi,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("X", "Y"))
gls_h4 # 


dt_est_h4 <- extract_gls_estimates(gls_h4, part = TRUE)
p_h4 <- dt_est_h4 %>% 
  ggplot() + 
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H4: evi change ~\nprotection", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none")
p_h4

## get estimates

## Hypothesis 5 - PA size and age ----------------------------

# recalculate trends only on PAs 
evi_cols <- grep("evi_", names(dt), value = T)

#subset to complete cases
dt_pa <- dt %>% filter(protection_cat_broad %in% c("Strict", "Mixed")) %>%
  dplyr::select(all_of(evi_cols), functional_biome, X, Y, unique_id,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, 
                protection_cat_broad, area_km2_log, pa_age_log) %>% 
  filter(complete.cases(.))  %>%
  distinct(X, Y, .keep_all = TRUE)

#get y matrix
y_pa <- as.matrix(dt_pa[, evi_cols])

#get coordinate matrix 
coords_pa <- as.matrix(dt_pa[, c("X", "Y")])

#fit autoregression, accounting for temporal autocorrelation 
ar_pa <- fitAR_map(Y = y_pa, coords = coords_pa)

#extract coefficients and p values 
dt_pa$evi_coef <- ar_pa$coefficients[, "t"]
dt_pa$abs_evi_coef <- abs(ar_pa$coefficients[, "t"])
dt_pa$evi_p_value <- ar_pa$pvals[, 2]

# get distance matrix 
#d_pa <- distm_scaled(coords_pa)

### estimate optimal r parameter (range of spatial autocorrelation)
corfit_pa <- fitCor(resids = residuals(ar_pa), coords = coords_pa, covar_FUN = "covar_exp", 
                    start = list(range = 0.1), fit.n = 10000)
(range_opt_pa = corfit_pa$spcor)

#use r to calculate optimal v parameter 
#v_opt_pa <- covar_exp(d_pa, range_opt_pa)

#partition the data 

pm_pa <- sample_partitions(npix = nrow(dt_pa), partsize = 1500, npart = NA)
dim(pm_pa)

gls_h5 <- fitGLS_partition(evi_coef ~ 1 +
                             area_km2_log + 
                             pa_age_log,
                           partmat = pm_pa,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_evi),
                           data = dt_pa,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("X", "Y"))
gls_h5 

dt_est_h5 <- extract_gls_estimates(gls_h5, part = TRUE)

p_h5 <- dt_est_h5 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H5: evi change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none")
p_h5


library(gridExtra)

p_evi_a <- grid.arrange(p_h2, p_h3, p_h4, p_h5, ncol = 4)
ggsave(plot = p_evi_a, "builds/plots/grid_evi_remotePARTS.png", dpi = 600, height = 3, width = 13)


# Biomes separately -------------------
# test the global change hypothesis again for all biomes separately -

# functional biomes 
## TMN, TMC, TMB, THN, SMD

biome_gls <- function(funbi = NA, col_pattern = NA, dat = NA, part = NA, part_size = NA){
  
  evi_cols <- grep(col_pattern, names(dt), value = T)
  
  #subset to complete cases
  dt_biome <- dat %>% filter(functional_biome == funbi) %>%
    dplyr::select(all_of(evi_cols), functional_biome, X, Y,
                  nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification) %>% 
    filter(complete.cases(.))
  
  #get y matrix
  y_biome <- as.matrix(dt_biome[, evi_cols])
  
  #get coordinate matrix 
  coords_biome <- as.matrix(dt_biome[, c("X", "Y")])
  
  #fit autoregression, accounting for temporal autocorrelation 
  ar_biome <- fitAR_map(Y = y_biome, coords = coords_biome)
  
  #extract coefficients and p values 
  dt_biome$evi_coef <- ar_biome$coefficients[, "t"]
  dt_biome$evi_p_value <- ar_biome$pvals[, 2]
  
  # get distance matrix 
  d_biome <- distm_scaled(coords_biome)
  
  ### estimate optimal r parameter (range of spatial autocorrelation)
  corfit_biome <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp", 
                         start = list(range = 0.1), fit.n = nrow(dt_biome))
  
  (range_opt_biome = corfit_biome$spcor)
  
  print("estimated optimal range")
  
  #use r to calculate optimal v parameter 
  v_opt_biome <- covar_exp(d_biome, range_opt_biome)
  
  cor.test(dt_biome$human_modification, dt_biome$nitrogen_depo)
  
  set.seed(161)
  gls_biome <- fitGLS(evi_coef ~ 1 +
                        nitrogen_depo +
                        mat_coef +
                        max_temp_coef +
                        map_coef +
                        human_modification,
                      data = dt_biome,
                      V = v_opt_biome, 
                      nugget = NA, 
                      no.F = TRUE
  )
  gls_biome 
  
  dt_est_biome <- extract_gls_estimates(gls_biome, part = FALSE)
  
  p_biome <- dt_est_biome %>% 
    ggplot() + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                    alpha = 0.9, linewidth = 1.2) +
    scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
    labs(title = paste0(funbi), subtitle = paste0("n = ", nrow(dt_biome)), y = NULL, x = NULL) +
    theme_classic() +
    theme(legend.position = "none")
  p_biome
  
  return(p_biome)
  
}



p_smd <- biome_gls(funbi = "SMD", col_pattern = "evi_", dat = dt)
p_smd

p_tmn <- biome_gls(funbi = "TMN", col_pattern = "evi_", dat = dt)
p_tmn

p_thn <- biome_gls(funbi = "THN", col_pattern = "evi_", dat = dt)
p_thn

p_tmb <- biome_gls(funbi = "TMB", col_pattern = "evi_", dat = dt)
p_tmb

p_tmc <- biome_gls(funbi = "TMC", col_pattern = "evi_", dat = dt)
p_tmc

p_evi_biome <- grid.arrange(p_smd, p_tmn, 
                            p_thn, 
                            p_tmb, p_tmc, ncol = 5)
ggsave(plot = p_evi_biome, "builds/plots/evi_remotePARTS_biome.png", dpi = 600, height = 3, width = 14)

p_evi <- grid.arrange(p_evi_a, p_evi_biome, ncol = 1)
ggsave(plot = p_evi, "builds/plots/evi_remotePARTS_full.png", dpi = 600, height = 6, width = 14)
