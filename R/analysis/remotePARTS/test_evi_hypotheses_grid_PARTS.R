#### test EVI hypothesis ####
source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") %>% 
  as.data.frame()
  
dt_raw <- dt_raw %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "Mixed",
      iucn_cat == "unprotected" ~ "Unprotected"), 
    super_biome = case_when(
      grepl("C", functional_biome) & grepl("T", functional_biome) ~ "cold_tall", 
      grepl("C", functional_biome) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short"),
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

####### calculate EVI trend #######

evi_cols <- grep("evi_", names(dt), value = T)

#subset to complete cases
dt_evi <- dt %>%
  dplyr::select(all_of(evi_cols), functional_biome, X, Y, lon, lat, unique_id, 
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.))

table(dt_evi$protection_cat_broad)
table(dt_evi$super_biome)


#get y matrix
y_evi <- as.matrix(dt_evi[, evi_cols])

#get coordinate matrix 
coords_evi <- as.matrix(dt_evi[, c("lon", "lat")])

#fit autoregression, accounting for temporal autocorrelation 
ar_evi <- fitAR_map(Y = y_evi, coords = coords_evi)

#extract coefficients and p values 
dt_evi$evi_coef <- ar_evi$coefficients[, "t"]
dt_evi$abs_evi_coef <- abs(ar_evi$coefficients[, "t"])
dt_evi$evi_p_value <- ar_evi$pvals[, 2]

fwrite(dt_evi %>% dplyr::select(
  unique_id, evi_coef, abs_evi_coef, evi_p_value), "data/processedData/dataFragments/grid_evi_trends.csv")

# get distance matrix 
#d_evi <- distm_scaled(coords_evi)

### estimate optimal r parameter (range of spatial autocorrelation)
corfit <- fitCor(resids = residuals(ar_evi), coords = coords_evi, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = 15000)
(range_opt_evi = corfit$spcor)

#use r to calculate optimal v parameter 
#v_opt_evi <- covar_exp(d_evi, range_opt_evi)

############ test hypotheses ############ 

#partition the data 
set.seed(161)
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
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = T, 
                           coord.names = c("lon", "lat")
)
gls_h1 # yes. Est: 5.563445 ; SE: 0.4952674 ; pval.t: 2.871102e-29

dt_est_h1 <- extract_gls_estimates(gls_h1, part = TRUE)

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
                           coord.names = c("lon", "lat"))
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

## Hypothesis 3 - different trends in different superbiomes  ----------------------
set.seed(161)
gls_h3 <- fitGLS_partition(evi_coef ~ 0 +
                             super_biome,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_evi),
                           data = dt_evi,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("lon", "lat"))
gls_h3 # 


dt_est_h3 <- extract_gls_estimates(gls_h3, part = TRUE)

p_h3 <- dt_est_h3 %>% 
  ggplot() + 
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H3: evi change ~\nsuper biome", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h3

## Hypothesis 4 - different absolute trend in and outside of PAs ----------------------
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
                           coord.names = c("lon", "lat"))
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
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h4


## Hypothesis 4.1 - different absolute trend in and outside of PAs ------------------
set.seed(161)
gls_h4.1 <- fitGLS_partition(abs_evi_coef ~ 0 +
                               protection_cat_broad,
                             partmat = pm,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = range_opt_evi),
                             data = dt_evi,
                             nugget = NA,
                             ncores = 25,
                             progressbar = TRUE, 
                             parallel = TRUE, 
                             coord.names = c("lon", "lat"))
gls_h4.1 # 


dt_est_h4.1 <- extract_gls_estimates(gls_h4.1, part = TRUE)

p_h4.1 <- dt_est_h4.1 %>% 
  ggplot() + 
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H4.1: abs evi change ~\nprotection", subtitle = paste0("n = ", nrow(dt_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h4.1


## Hypothesis 5 - PA characteristics ----------------------------

# recalculate trends only on PAs 
evi_cols <- grep("evi_", names(dt), value = T)

#subset to complete cases
dt_pa <- dt %>% filter(protection_cat_broad == "Strict") %>%
  dplyr::select(all_of(evi_cols),functional_biome, X, Y, lon, lat, unique_id, 
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad, area_km2_log, pa_age_log) %>% 
  filter(complete.cases(.))

#get y matrix
y_pa <- as.matrix(dt_pa[, evi_cols])

#get coordinate matrix 
coords_pa <- as.matrix(dt_pa[, c("lon", "lat")])

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
                    start = list(range = 0.1), fit.n = 15000)
(range_opt_pa = corfit_pa$spcor)

#use r to calculate optimal v parameter 
#v_opt_pa <- covar_exp(d_pa, range_opt_pa)

#partition the data 

pm_pa <- sample_partitions(npix = nrow(dt_pa), partsize = 1500, npart = NA)
dim(pm_pa)

# change --- 
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
                           coord.names = c("lon", "lat"))
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
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h5

## Hypothesis 5.1 - absolute change PA characteristics ----------------------------
gls_h5.1 <- fitGLS_partition(abs_evi_coef ~ 1 +
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
                             coord.names = c("lon", "lat"))
gls_h5.1 

dt_est_h5.1 <- extract_gls_estimates(gls_h5.1, part = TRUE)

p_h5.1 <- dt_est_h5.1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  labs(title = "H5.1: abs evi change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h5.1

library(gridExtra)
p_evi_a <- grid.arrange(p_h2, p_h3, p_h4, p_h5, ncol = 4)
ggsave(plot = p_evi_a, "builds/plots/evi_remotePARTS_grid.png", dpi = 600, height = 3, width = 12)


# Super Biomes separately -------------------
# test the global change hypothesis again for all biomes separately -
unique(dt$super_biome)


#super_b <- "not_cold_tall"
#col_pattern <- "evi_"
#dat = dt

biome_gls <- function(super_b = NA, col_pattern = NA, start = list(range = 0.1),
                      dat = NA, part = TRUE, part_size = NA, fit_n = "row_n"){
  
  coi <- grep(col_pattern, names(dat), value = T)
  
  dat <- dat %>% filter(super_biome %in% super_b) %>%
    dplyr::select(all_of(coi), functional_biome, X, Y, super_biome, lon, lat, 
                  nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification) %>% 
    filter(complete.cases(.))
  
  dt_biome <- dat 
  
  #get y matrix
  y_biome <- as.matrix(dt_biome[, coi])
  
  #get coordinate matrix 
  coords_biome <- as.matrix(dt_biome[, c("lon", "lat")])
  
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
    
    print(paste0("start estimating range on a subset of ", fit_n, " datapoints"))
    
    ### estimate optimal r parameter (range of spatial autocorrelation)
    corfit_biome1 <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp", 
                            start = start, fit.n = fit_n)
    
    r1 <-  corfit_biome1$spcor
    
   # corfit_biome2 <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp", 
   #                          start = start, fit.n = fit_n)
   # 
   # r2 <-  corfit_biome2$spcor
   # 
   # corfit_biome3 <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp", 
   #                        start = start, fit.n = fit_n)
   #
   # r3 <-  corfit_biome3$spcor
   # 
   # (range_opt_biome <- mean(c(r1, r2, r3)))
    (range_opt_biome <- r1)
      
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
    
    if(is.na(part_size)){part_size <- 1000}
    
    pm_biome <- sample_partitions(npix = nrow(dt_biome), partsize = part_size, npart = NA)
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
                                  ncores = 25,
                                  progressbar = TRUE, 
                                  parallel = TRUE, 
                                  coord.names = c("lon", "lat"))
    gls_biome 
    
    dt_est_biome <- extract_gls_estimates(gls_biome, part = TRUE)
  }
  
  
  label <- paste0("Super Biome:\n", super_b)
  
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

table(dt_raw[dt_raw$super_biome == "cold_tall", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "cold_short", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "not_cold_tall", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "not_cold_short", ]$functional_biome)


nrow(dt[dt$super_biome == "cold_tall",])
p_cold_tall <- biome_gls(super_b = "cold_tall", fit_n = 15000, part_size = 1500,
                         col_pattern = "evi_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_cold_tall <- p_cold_tall$data
p_cold_tall

nrow(dt[dt$super_biome == "cold_short",])
p_cold_short <- biome_gls(super_b = "cold_short", fit_n = 15000, part_size = 1500,
                          col_pattern = "evi_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_cold_short <- p_cold_short$data
p_cold_short

nrow(dt[dt$super_biome == "not_cold_tall",])
p_not_cold_tall <- biome_gls(super_b = "not_cold_tall", fit_n = 15000, part_size = 1500,
                             col_pattern = "evi_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_not_cold_tall <- p_not_cold_tall$data
p_not_cold_tall

nrow(dt[dt$super_biome == "not_cold_short",])
p_not_cold_short <- biome_gls(super_b = "not_cold_short", fit_n = 15000, part_size = 1500,
                              col_pattern = "evi_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_not_cold_short <- p_not_cold_short$data
p_not_cold_short

p_evi_super_b <- gridExtra::grid.arrange(p_cold_tall, p_cold_short, p_not_cold_tall, p_not_cold_short,
                                         ncol = 4)
ggsave(plot = p_evi_super_b, "builds/plots/evi_remotePARTS_super_b_grid.png", dpi = 600, height = 3, width = 12)

p_evi <- gridExtra::grid.arrange(p_evi_a, p_evi_super_b, ncol = 1)
ggsave(plot = p_evi, "builds/plots/evi_remotePARTS_full_grid.png", dpi = 600, height = 5, width = 14)


#### combine and write out the estimates 

dt_est <- rbind(
  dt_est_h1 %>% mutate(model = "H1", t.stat = NA, response = NA),
  dt_est_h2 %>% mutate(model = "H2", t.stat = NA, response = NA), 
  dt_est_h3 %>% mutate(model = "H3", t.stat = NA, response = NA), 
  dt_est_h4 %>% mutate(model = "H4", t.stat = NA, response = NA), 
  dt_est_h4.1 %>% mutate(model = "H4.1", t.stat = NA, response = NA), 
  dt_est_h5 %>% mutate(model = "H5", t.stat = NA, response = NA), 
  dt_est_h5.1 %>% mutate(model = "H5.1", t.stat = NA, response = NA),
  dt_est_cold_tall %>% mutate(model = "cold_tall", t.stat = NA, response = NA), 
  dt_est_cold_short %>% mutate(model = "cold_short", t.stat = NA, response = NA), 
  dt_est_not_cold_tall %>% mutate(model = "not_cold_tall", t.stat = NA, response = NA), 
  dt_est_not_cold_short %>% mutate(model = "not_cold_short", t.stat = NA, response = NA)
)

fwrite(dt_est, "builds/model_estimates/evi_grid_estimates.csv")
