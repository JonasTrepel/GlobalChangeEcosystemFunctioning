#### test EVI hypothesis ####
source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  group_by(FunctionalBiome) %>% 
  mutate(n_per_functional_biome = n()) %>% 
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
    nitrogen_depo = scale(NitrogenDepo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(HumanModification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = scale(log(PaAge + 0.0001))) 

dt <- dt_raw

##### Subset for testing only #######

# dt_sub_pa <- dt_raw %>%
#   filter(og_layer == "protected_areas") %>%
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
  dplyr::select(all_of(evi_cols), functional_biome, X, Y, unique_id, productivity, ndvi_min,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, 
                og_layer) %>% 
  filter(complete.cases(.))

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
  unique_id, X, Y, og_layer, evi_coef, abs_evi_coef, evi_p_value), "data/processedData/dataFragments/pa_evi_trends.csv")

# get distance matrix 
#d_evi <- distm_scaled(coords_evi)

### estimate optimal r parameter (range of spatial autocorrelation)
corfit <- fitCor(resids = residuals(ar_evi), coords = coords_evi, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = nrow(dt_evi))
(range_opt_evi = corfit$spcor)

#use r to calculate optimal v parameter 
#v_opt_evi <- covar_exp(d_evi, range_opt_evi)

############ test hypotheses ############ 

#partition the data 

pm <- sample_partitions(npix = nrow(dt_evi), partsize = 1500, npart = NA)
dim(pm)

## Hypothesis 1: Ecosystem functioning is overall changing --------------------
set.seed(161)
gls_h1 <- fitGLS_partition(evi_coef ~ 1,
               partmat = pm,
               covar_FUN = "covar_exp",
               covar.pars = list(range = range_opt_evi),
               data = dt_evi,
               nugget = NA,
               ncores = 4,
               progressbar = TRUE, 
               parallel = T, 
               coord.names = c("X", "Y")
               )
gls_h1 # yes. Est: 6.489292; SE: 0.665201; pval.t: 2.039496e-22

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
                   ncores = 4,
                   progressbar = TRUE, 
                   parallel = TRUE, 
                   coord.names = c("X", "Y"))
gls_h2 

dt_est_h2 <- extract_gls_estimates(gls_h2, part = TRUE)

p_h2 <- dt_est_h2 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                                        alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H2: evi change ~\nglobal change", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h2
  
## Hypothesis 3 - different trends at different productivity  ----------------------
set.seed(161)
gls_h3 <- fitGLS_partition(evi_coef ~ 0 +
                   productivity,
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
  labs(title = "H3: evi change ~\nproductivity", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h3

## Hypothesis 4 - different trends in different seasonality ----------------------
set.seed(161)
gls_h4 <- fitGLS_partition(evi_coef ~ 0 +
                           ndvi_min,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_evi),
                           data = dt_evi,
                           nugget = NA,
                           ncores = 4,
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
  labs(title = "H4: evi change ~\ngrowth lim.", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h4


## Hypothesis 5 - different absolute trend in and outside of PAs ------------------
set.seed(161)
gls_h5 <- fitGLS_partition(evi_coef ~ 0 +
                 og_layer,
                 partmat = pm,
                 covar_FUN = "covar_exp",
                 covar.pars = list(range = range_opt_evi),
                 data = dt_evi,
                 nugget = NA,
                 ncores = 4,
                 progressbar = TRUE, 
                 parallel = TRUE, 
                 coord.names = c("X", "Y"))
gls_h5 # 


dt_est_h5 <- extract_gls_estimates(gls_h5, part = TRUE)

p_h5 <- dt_est_h5 %>% 
  ggplot() + 
 # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H4: evi change ~\nprotection", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h5

## get estimates

## Hypothesis 6 - PA size and age ----------------------------

# recalculate trends only on PAs 
evi_cols <- grep("evi_", names(dt), value = T)

#subset to complete cases
dt_pa <- dt %>% filter(og_layer == "protected_areas") %>%
  dplyr::select(all_of(evi_cols), functional_biome, X, Y,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, 
                og_layer, area_km2_log, pa_age_log) %>% 
  filter(complete.cases(.))

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
                 start = list(range = 0.1), fit.n = nrow(dt_pa))
(range_opt_pa = corfit_pa$spcor)

#use r to calculate optimal v parameter 
#v_opt_pa <- covar_exp(d_pa, range_opt_pa)

#partition the data 

pm_pa <- sample_partitions(npix = nrow(dt_pa), partsize = 1000, npart = NA)
dim(pm_pa)

gls_h6 <- fitGLS_partition(evi_coef ~ 1 +
                   area_km2_log + 
                   pa_age_log,
                   partmat = pm_pa,
                   covar_FUN = "covar_exp",
                   covar.pars = list(range = range_opt_evi),
                   data = dt_pa,
                   nugget = NA,
                   ncores = 4,
                   progressbar = TRUE, 
                   parallel = TRUE, 
                   coord.names = c("X", "Y"))
gls_h6 

dt_est_h6 <- extract_gls_estimates(gls_h6, part = TRUE)

p_h6 <- dt_est_h6 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H6: evi change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h6


library(gridExtra)

p_evi_a <- grid.arrange(p_h2, p_h3, p_h4, p_h5, p_h6, ncol = 5)
ggsave(plot = p_evi_a, "builds/plots/evi_remotePARTS.png", dpi = 600, height = 3, width = 13)


# Biomes separately -------------------
# test the global change hypothesis again for all biomes separately -
unique(dt$ndvi_min)
unique(dt$productivity)

# unctional biomes 
## TMN, TMC, TMB, THN, SMD

prod <- "medium"
ndvi_m <- NA
col_pattern <- "evi_"
dat = dt

biome_gls <- function(ndvi_m = NA, prod = NA, col_pattern = NA, start = list(range = 0.1),
                      dat = NA, part = TRUE, part_size = NA, fit_n = "row_n"){
  
  evi_cols <- grep(col_pattern, names(dt), value = T)
  
  
  if(!is.na(prod)){
    dat <- dat %>% filter(productivity %in% prod) %>%
      #filter(BurnedAreaMean > 0) %>% 
      dplyr::select(all_of(evi_cols), functional_biome, X, Y, ndvi_min,
                    nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification) %>% 
      filter(complete.cases(.))
  }
  
  if(!is.na(ndvi_m)){
    dat <- dat %>% filter(ndvi_min %in% ndvi_m) %>%
      #filter(BurnedAreaMean > 0) %>% 
      dplyr::select(all_of(evi_cols), functional_biome, X, Y,
                    nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification) %>% 
      filter(complete.cases(.))
  }
  
  dt_biome <- dat 
  #get y matrix
  y_biome <- as.matrix(dt_biome[, evi_cols])
  
  #get coordinate matrix 
  coords_biome <- as.matrix(dt_biome[, c("X", "Y")])
  
  #fit autoregression, accounting for temporal autocorrelation 
  ar_biome <- fitAR_map(Y = y_biome, coords = coords_biome)
  
  #extract coefficients and p values 
  dt_biome$evi_coef <- ar_biome$coefficients[, "t"]
  dt_biome$evi_p_value <- ar_biome$pvals[, 2]
  dt_biome$abs_evi_coef <- abs(ar_biome$coefficients[, "t"])
  
  if(fit_n == "row_n"){
    
    print("start estimating range on full dataset")
    
    fit_n <- nrow(dt_biome)
    
    ### estimate optimal r parameter (range of spatial autocorrelation)
    corfit_biome <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp", 
                               start = start, fit.n = fit_n)
    
    (range_opt_biome = corfit_biome$spcor)
    
    
  }else{
    
    print("start estimating range on three different subsets of the data")
    
    ### estimate optimal r parameter (range of spatial autocorrelation)
    corfit_biome1 <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp", 
                            start = start, fit.n = fit_n)
    
    r1 <-  corfit_biome1$spcor
    
    corfit_biome2 <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp", 
                            start = start, fit.n = fit_n)
    
    r2 <-  corfit_biome2$spcor
    
    corfit_biome3 <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp", 
                            start = start, fit.n = fit_n)
    
    r3 <-  corfit_biome3$spcor
    
    (range_opt_biome <- mean(c(r1, r2, r3)))
    
  }
  
  
  print(paste0("estimated optimal range (", round(range_opt_biome, 4),"), start fitting GLS"))
  
  if(part == FALSE){
    
    # get distance matrix 
    d_biome <- distm_scaled(coords_biome)
    
    #use r to calculate optimal v parameter 
    v_opt_biome <- covar_exp(d_biome, range_opt_biome)
    
    
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
    
  }
  
  if(part == TRUE){
    
    pm_biome <- sample_partitions(npix = nrow(dt_biome), partsize = 1000, npart = NA)
    dim(pm_biome)
    
    
    set.seed(161)
    gls_biome <- fitGLS_partition(evi_coef ~ 1 +
                                    nitrogen_depo +
                                    mat_coef + 
                                    max_temp_coef + 
                                    map_coef + 
                                    human_modification,
                                  partmat = pm_biome,
                                  covar_FUN = "covar_exp",
                                  covar.pars = list(range = range_opt_biome),
                                  data = dt_biome,
                                  nugget = NA,
                                  ncores = 4,
                                  progressbar = TRUE, 
                                  parallel = TRUE, 
                                  coord.names = c("X", "Y"))
    gls_biome 
    
    dt_est_biome <- extract_gls_estimates(gls_biome, part = TRUE)
  }
  
  
  if(!is.na(ndvi_m)){
    label <- paste0("NDVI min:\n", ndvi_m)
  }else{
    label <- paste0("Productivity:\n", prod)
    
  }
  
  p_biome <- dt_est_biome %>% 
    ggplot() + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                    alpha = 0.9, linewidth = 1.2) +
    scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
    labs(title = paste0(label), subtitle = paste0("n = ", nrow(dt_biome)), y = NULL, x = NULL) +
    theme_classic() +
    theme(legend.position = "none", 
          plot.title = element_text(size = 12))
  p_biome
  
  return(p_biome)
  
}

## NDVI min 
nrow(dt[dt$ndvi_min == "cold",])
p_cold <- biome_gls(ndvi_m = "cold", col_pattern = "evi_", dat = dt, part = TRUE, start = list(range = 0.1))
p_cold

nrow(dt[dt$ndvi_min == "dry",])
p_dry <- biome_gls(ndvi_m = "dry", col_pattern = "evi_", dat = dt, part = FALSE)
p_dry

nrow(dt[dt$ndvi_min == "cold_and_dry",])
p_cold_and_dry <- biome_gls(ndvi_m = "cold_and_dry", col_pattern = "evi_", dat = dt, part = FALSE)
p_cold_and_dry

nrow(dt[dt$ndvi_min == "non_seasonal",])
p_non_seasonal <- biome_gls(ndvi_m = "non_seasonal", col_pattern = "evi_", dat = dt, part = TRUE)
p_non_seasonal

## 
nrow(dt[dt$productivity == "low",])
p_low <- biome_gls(prod = "low", col_pattern = "evi_", dat = dt, part = FALSE)
p_low

nrow(dt[dt$productivity == "medium",])
p_medium <- biome_gls(prod = "medium", col_pattern = "evi_", dat = dt, part = TRUE)
p_medium

nrow(dt[dt$productivity == "high",])
p_high <- biome_gls(prod = "high", col_pattern = "evi_", dat = dt, part = FALSE)
p_high


p_evi_ndvi <- gridExtra::grid.arrange(p_cold, p_dry, 
                      p_cold_and_dry, p_non_seasonal, ncol = 5)
ggsave(plot = p_evi_ndvi, "builds/plots/evi_remotePARTS_ndvi.png", dpi = 600, height = 3, width = 14)

p_evi_prod <- gridExtra::grid.arrange(p_low, p_medium, p_high, ncol = 5)
ggsave(plot = p_evi_prod, "builds/plots/evi_remotePARTS_prod.png", dpi = 600, height = 3, width = 14)

p_evi <- gridExtra::grid.arrange(p_evi_a, p_evi_ndvi, p_evi_prod, ncol = 1)
ggsave(plot = p_evi, "builds/plots/evi_remotePARTS_full.png", dpi = 600, height = 8, width = 14)
