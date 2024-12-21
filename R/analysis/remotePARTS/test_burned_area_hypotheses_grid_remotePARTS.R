#### test burned_area hypothesis ####
source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") %>% 
  group_by(FunctionalBiome) %>% 
  mutate(n_per_functional_biome = n()) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  # filter(n_per_functional_biome > 1000) %>% 
  rename(functional_biome = FunctionalBiomeShort) %>% 
  mutate(
    PaAge = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
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
    pa_age_log = scale(log(PaAge + 0.0001))) 

table(dt_raw[dt_raw$super_biome == "cold_tall", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "cold_short", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "not_cold_tall", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "not_cold_short", ]$functional_biome)

dt <- dt_raw

#######################################

####### calculate burned_area trend #######

burned_area_cols <- grep("burned_area_", names(dt), value = T)

#subset to complete cases
dt_burned_area <- dt %>%
  dplyr::select(all_of(burned_area_cols), functional_biome, X, Y, unique_id, productivity, ndvi_min,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.))

#get y matrix
y_burned_area <- as.matrix(dt_burned_area[, burned_area_cols])

#get coordinate matrix 
coords_burned_area <- as.matrix(dt_burned_area[, c("X", "Y")])

#fit autoregression, accounting for temporal autocorrelation 
ar_burned_area <- fitAR_map(Y = y_burned_area, coords = coords_burned_area)

#extract coefficients and p values 
dt_burned_area$burned_area_coef <- ar_burned_area$coefficients[, "t"]
dt_burned_area$abs_burned_area_coef <- abs(ar_burned_area$coefficients[, "t"])
dt_burned_area$burned_area_p_value <- ar_burned_area$pvals[, 2]

fwrite(dt_burned_area %>% dplyr::select(
  unique_id, X, Y, protection_cat_broad, burned_area_coef, abs_burned_area_coef, burned_area_p_value), "data/processedData/dataFragments/grid_burned_area_trends.csv")

# get distance matrix 
#d_burned_area <- distm_scaled(coords_burned_area)

### estimate optimal r parameter (range of spatial autocorrelation)
corfit <- fitCor(resids = residuals(ar_burned_area), coords = coords_burned_area, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = 15000)
(range_opt_burned_area = corfit$spcor)

#use r to calculate optimal v parameter 
#v_opt_burned_area <- covar_exp(d_burned_area, range_opt_burned_area)

############ test hypotheses ############ 

#partition the data 
set.seed(161)
pm <- sample_partitions(npix = nrow(dt_burned_area), partsize = 1500, npart = NA)
dim(pm)

## Hypothesis 1: Ecosystem functioning is overall changing --------------------
set.seed(161)
gls_h1 <- fitGLS_partition(burned_area_coef ~ 1,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_burned_area),
                           data = dt_burned_area,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = T, 
                           coord.names = c("X", "Y")
)
gls_h1 # yes. Est: -0.0002312834; SE: 3.678995e-05; pval.t: 3.253132e-10

dt_est_h1 <- extract_gls_estimates(gls_h1, part = TRUE)

## Hypothesis 2: Change depends on climate change, N deposition, human modification -------------

set.seed(161)
gls_h2 <- fitGLS_partition(burned_area_coef ~ 1 +
                             nitrogen_depo +
                             mat_coef + 
                             max_temp_coef + 
                             map_coef + 
                             human_modification,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_burned_area),
                           data = dt_burned_area,
                           nugget = NA,
                           ncores = 25,
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
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H2: burned_area change ~\nglobal change", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h2

## Hypothesis 3 - different trends in different superbiomes  ----------------------
set.seed(161)
gls_h3 <- fitGLS_partition(burned_area_coef ~ 0 +
                             super_biome,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_burned_area),
                           data = dt_burned_area,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("X", "Y"))
gls_h3 # 


dt_est_h3 <- extract_gls_estimates(gls_h3, part = TRUE)

p_h3 <- dt_est_h3 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H3: burned_area change ~\nsuper biome", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h3

## Hypothesis 4 - different absolute trend in and outside of PAs ----------------------
set.seed(161)
gls_h4 <- fitGLS_partition(burned_area_coef ~ 0 +
                             protection_cat_broad,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_burned_area),
                           data = dt_burned_area,
                           nugget = NA,
                           ncores = 25,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("X", "Y"))
gls_h4 # 


dt_est_h4 <- extract_gls_estimates(gls_h4, part = TRUE)
p_h4 <- dt_est_h4 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H4: burned_area change ~\nprotection", subtitle = paste0("n = ", nrow(dt_burned_area)), y = NULL, x = NULL) +
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
                             data = dt_burned_area,
                             nugget = NA,
                             ncores = 25,
                             progressbar = TRUE, 
                             parallel = TRUE, 
                             coord.names = c("X", "Y"))
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
dt_pa <- dt %>% filter(protection_cat_broad == "Strict") %>%
  dplyr::select(all_of(burned_area_cols), functional_biome, X, Y, super_biome,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, 
                protection_cat_broad, area_km2_log, pa_age_log) %>% 
  filter(complete.cases(.))

#get y matrix
y_pa <- as.matrix(dt_pa[, burned_area_cols])

#get coordinate matrix 
coords_pa <- as.matrix(dt_pa[, c("X", "Y")])

#fit autoregression, accounting for temporal autocorrelation 
ar_pa <- fitAR_map(Y = y_pa, coords = coords_pa)

#extract coefficients and p values 
dt_pa$burned_area_coef <- ar_pa$coefficients[, "t"]
dt_pa$abs_burned_area_coef <- abs(ar_pa$coefficients[, "t"])
dt_pa$burned_area_p_value <- ar_pa$pvals[, 2]

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

# change --- 
gls_h5 <- fitGLS_partition(burned_area_coef ~ 1 +
                             area_km2_log + 
                             pa_age_log,
                           partmat = pm_pa,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_burned_area),
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
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H5: burned_area change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h5

## Hypothesis 5.1 - absolute change PA characteristics ----------------------------
gls_h5.1 <- fitGLS_partition(abs_burned_area_coef ~ 1 +
                               area_km2_log + 
                               pa_age_log,
                             partmat = pm_pa,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = range_opt_burned_area),
                             data = dt_pa,
                             nugget = NA,
                             ncores = 25,
                             progressbar = TRUE, 
                             parallel = TRUE, 
                             coord.names = c("X", "Y"))
gls_h5.1 

dt_est_h5.1 <- extract_gls_estimates(gls_h5.1, part = TRUE)

p_h5.1 <- dt_est_h5.1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
  labs(title = "H5.1: abs burned_area change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h5.1

library(gridExtra)
p_burned_area_a <- grid.arrange(p_h2, p_h3, p_h4, p_h5, ncol = 4)
ggsave(plot = p_burned_area_a, "builds/plots/burned_area_remotePARTS_grid.png", dpi = 600, height = 3, width = 12)


# Super Biomes separately -------------------
# test the global change hypothesis again for all biomes separately -
unique(dt$super_biome)


super_b <- "not_cold_tall"
col_pattern <- "burned_area_"
dat = dt

biome_gls <- function(super_b = NA, col_pattern = NA, start = list(range = 0.1),
                      dat = NA, part = TRUE, part_size = NA, fit_n = "row_n"){
  
  coi <- grep(col_pattern, names(dat), value = T)
  
  dat <- dat %>% filter(super_biome %in% super_b) %>%
    dplyr::select(all_of(coi), functional_biome, X, Y, ndvi_min,
                  nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification) %>% 
    filter(complete.cases(.))
  
  dt_biome <- dat 
  
  #get y matrix
  y_biome <- as.matrix(dt_biome[, coi])
  
  #get coordinate matrix 
  coords_biome <- as.matrix(dt_biome[, c("X", "Y")])
  
  #fit autoregression, accounting for temporal autocorrelation 
  ar_biome <- fitAR_map(Y = y_biome, coords = coords_biome)
  
  #extract coefficients and p values 
  dt_biome$burned_area_coef <- ar_biome$coefficients[, "t"]
  dt_biome$burned_area_p_value <- ar_biome$pvals[, 2]
  dt_biome$abs_burned_area_coef <- abs(ar_biome$coefficients[, "t"])
  
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
    gls_biome <- fitGLS(burned_area_coef ~ 1 +
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
    gls_biome <- fitGLS_partition(burned_area_coef ~ 1 +
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
                                  coord.names = c("X", "Y"))
    gls_biome 
    
    dt_est_biome <- extract_gls_estimates(gls_biome, part = TRUE)
  }
  
  
  label <- paste0("Super Biome:\n", super_b)
  
  p_biome <- dt_est_biome %>% 
    ggplot() + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                    alpha = 0.9, linewidth = 1.2) +
    scale_color_manual(values = c("significant" = "firebrick", "non-significant" = "grey")) +
    labs(title = paste0(label), subtitle = paste0("n = ", nrow(dt_biome)), y = NULL, x = NULL) +
    theme_classic() +
    theme(legend.position = "none", 
          plot.title = element_text(size = 12))
  p_biome
  
  return(p_biome)
  
}

## NDVI min 
table(dt_raw[dt_raw$super_biome == "cold_tall", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "cold_short", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "not_cold_tall", ]$functional_biome)
table(dt_raw[dt_raw$super_biome == "not_cold_short", ]$functional_biome)


nrow(dt[dt$super_biome == "cold_tall",])
p_cold_tall <- biome_gls(super_b = "cold_tall", fit_n = 15000, part_size = 1500, 
                         col_pattern = "burned_area_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_cold_tall <- p_cold_tall$data
p_cold_tall

nrow(dt[dt$super_biome == "cold_short",])
p_cold_short <- biome_gls(super_b = "cold_short", fit_n = 15000, part_size = 1500,
                          col_pattern = "burned_area_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_cold_short <- p_cold_short$data
p_cold_short

nrow(dt[dt$super_biome == "not_cold_tall",])
p_not_cold_tall <- biome_gls(super_b = "not_cold_tall", fit_n = 15000, part_size = 1500,
                             col_pattern = "burned_area_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_not_cold_tall <- p_not_cold_tall$data
p_not_cold_tall

nrow(dt[dt$super_biome == "not_cold_short",])
p_not_cold_short <- biome_gls(super_b = "not_cold_short", fit_n = 15000, part_size = 1500,
                              col_pattern = "burned_area_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_not_cold_short <- p_not_cold_short$data
p_not_cold_short

p_burned_area_super_b <- gridExtra::grid.arrange(p_cold_tall, p_cold_short, p_not_cold_tall, p_not_cold_short,
                                         ncol = 4)
ggsave(plot = p_burned_area_super_b, "builds/plots/burned_area_remotePARTS_super_b_grid.png", dpi = 600, height = 3, width = 12)

p_burned_area <- gridExtra::grid.arrange(p_burned_area_a, p_burned_area_super_b, ncol = 1)
ggsave(plot = p_burned_area, "builds/plots/burned_area_remotePARTS_full_grid.png", dpi = 600, height = 5, width = 14)


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

fwrite(dt_est, "builds/model_estimates/burned_area_grid_estimates.csv")
