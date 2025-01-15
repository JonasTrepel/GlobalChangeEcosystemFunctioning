#### test mean_evi hypothesis ####
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

####### calculate mean_evi trend #######

mean_evi_cols <- grep("mean_evi_", names(dt), value = T)

#subset to complete cases
dt_mean_evi <- dt %>%
  dplyr::select(all_of(mean_evi_cols), functional_biome, lon, lat, unique_id, lon_pa, lat_pa,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad) %>% 
  filter(complete.cases(.)) %>% 
  distinct(lon, lat, .keep_all = T)

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

fwrite(dt_mean_evi %>% dplyr::select(
  unique_id, mean_evi_coef, abs_mean_evi_coef, mean_evi_p_value), "data/processedData/dataFragments/pas_mean_evi_trends.csv")

# get distance matrix 
#d_mean_evi <- distm_scaled(coords_mean_evi)

### estimate optimal r parameter (range of spatial autocorrelation)
corfit <- fitCor(resids = residuals(ar_mean_evi), coords = coords_mean_evi, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = 15000)
(range_opt_mean_evi = corfit$spcor)


##### by protected area 

get_mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

dt_sum <- dt_mean_evi %>%
  group_by(unique_pa_id) %>% 
  summarize(mat_coef = mean(mat_coef, na.rm = T),
            map_coef = mean(map_coef, na.rm = T),
            max_temp_coef = mean(max_temp_coef, na.rm = T),
            human_modification = mean(human_modification, na.rm = T),
            nitrogen_depo = mean(nitrogen_depo, na.rm = T),
            protection_cat_broad = unique(protection_cat_broad), 
            lon_pa = unique(lon_pa), 
            lat_pa = unique(lat_pa), 
            mean_evi_coef = mean(mean_evi_coef, na.rm = T),
            abs_mean_evi_coef = mean(abs_mean_evi_coef, na.rm = T)) 
  



############ test hypotheses ############ 

#partition the data 
set.seed(161)
pm <- sample_partitions(npix = nrow(dt_sum), partsize = 1500, npart = NA)
dim(pm)


## Hypothesis 4 - different absolute trend in and outside of PAs ----------------------
set.seed(161)
gls_h4 <- fitGLS_partition(mean_evi_coef ~ 0 +
                           protection_cat_broad,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_mean_evi),
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
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "H4: mean_evi change ~\nprotection", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h4


## Hypothesis 4.1 - different absolute trend in and outside of PAs ------------------
set.seed(161)
gls_h4.1 <- fitGLS_partition(abs_mean_evi_coef ~ 0 +
                              protection_cat_broad,
                             partmat = pm,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = range_opt_mean_evi),
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
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "H4.1: abs mean_evi change ~\nprotection", subtitle = paste0("n = ", nrow(dt_mean_evi)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))
p_h4.1


## Hypothesis 5 - PA characteristics ----------------------------

# recalculate trends only on PAs 
mean_evi_cols <- grep("mean_evi_", names(dt), value = T)

#subset to complete cases
dt_pa <- dt %>% filter(protection_cat_broad == "strictly_protected") %>%
  dplyr::select(all_of(mean_evi_cols),functional_biome, lon, lat, unique_id, lon_pa, lat_pa,
                nitrogen_depo, mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
                protection_cat_broad, area_km2_log, pa_age_log) %>% 
  filter(complete.cases(.))

#get y matrix
y_pa <- as.matrix(dt_pa[, mean_evi_cols])

#get coordinate matrix 
coords_pa <- as.matrix(dt_pa[, c("lon", "lat")])

#fit autoregression, accounting for temporal autocorrelation 
ar_pa <- fitAR_map(Y = y_pa, coords = coords_pa)

#extract coefficients and p values 
dt_pa$mean_evi_coef <- ar_pa$coefficients[, "t"]
dt_pa$abs_mean_evi_coef <- abs(ar_pa$coefficients[, "t"])
dt_pa$mean_evi_p_value <- ar_pa$pvals[, 2]

### estimate optimal r parameter (range of spatial autocorrelation)
corfit_pa <- fitCor(resids = residuals(ar_pa), coords = coords_pa, covar_FUN = "covar_exp", 
                    start = list(range = 0.1), fit.n = nrow(dt_pa))
(range_opt_pa = corfit_pa$spcor)

#summarize PA 

dt_pa_sum <- dt_mean_evi %>%
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
            mean_evi_coef = mean(mean_evi_coef, na.rm = T),
            abs_mean_evi_coef = mean(abs_mean_evi_coef, na.rm = T)) 


#partition the data 

pm_pa <- sample_partitions(npix = nrow(dt_pa_sum), partsize = 1000, npart = NA)
dim(pm_pa)

# change --- 
gls_h5 <- fitGLS_partition(mean_evi_coef ~ 1 +
                             area_km2_log + 
                             pa_age_log,
                           partmat = pm_pa,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = range_opt_mean_evi),
                           data = dt_pa_sum,
                           nugget = NA,
                           ncores = 10,
                           progressbar = TRUE, 
                           parallel = TRUE, 
                           coord.names = c("lon_pa", "lat_pa"))
gls_h5 

dt_est_h5 <- extract_gls_estimates(gls_h5, part = TRUE)

p_h5 <- dt_est_h5 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "H5: mean_evi change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h5

## Hypothesis 5.1 - absolute change PA characteristics ----------------------------
gls_h5.1 <- fitGLS_partition(abs_mean_evi_coef ~ 1 +
                               area_km2_log + 
                               pa_age_log,
                             partmat = pm_pa,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = range_opt_mean_evi),
                             data = dt_pa_sum,
                             nugget = NA,
                             ncores = 10,
                             progressbar = TRUE, 
                             parallel = TRUE, 
                             coord.names = c("lon_pa", "lat_pa"))
gls_h5.1 

dt_est_h5.1 <- extract_gls_estimates(gls_h5.1, part = TRUE)

p_h5.1 <- dt_est_h5.1 %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
  labs(title = "H5.1: abs mean_evi change ~\npa age & size", subtitle = paste0("n = ", nrow(dt_pa)), y = NULL, x = NULL) +
  theme_classic() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12)) 
p_h5.1

library(gridExtra)
p_mean_evi_a <- grid.arrange(p_h2, p_h3, p_h4, p_h5, ncol = 4)
ggsave(plot = p_mean_evi_a, "builds/plots/mean_evi_remotePARTS_pas.png", dpi = 600, height = 3, width = 12)


# Super Biomes separately -------------------
# test the global change hypothesis again for all biomes separately -
unique(dt$super_biome)


#super_b <- "not_cold_tall"
#col_pattern <- "mean_evi_"
#dat = dt

biome_gls <- function(super_b = NA, col_pattern = NA, start = list(range = 0.1),
                      dat = NA, part = TRUE, part_size = NA, fit_n = "row_n"){
  
  coi <- grep(col_pattern, names(dat), value = T)
  
  dat <- dat %>% filter(super_biome %in% super_b) %>%
    dplyr::select(all_of(coi), functional_biome, super_biome, lon, lat, 
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
  dt_biome$mean_evi_coef <- ar_biome$coefficients[, "t"]
  dt_biome$mean_evi_p_value <- ar_biome$pvals[, 2]
  dt_biome$abs_mean_evi_coef <- abs(ar_biome$coefficients[, "t"])
  
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
    gls_biome <- fitGLS(mean_evi_coef ~ 1 +
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
    gls_biome <- fitGLS_partition(mean_evi_coef ~ 1 +
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
                                  ncores = 10,
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
    scale_color_manual(values = c("significant" = "olivedrab", "non-significant" = "grey")) +
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
p_cold_tall <- biome_gls(super_b = "cold_tall", fit_n = "row_n", part_size = 1000,
                         col_pattern = "mean_evi_", dat = dt, part = TRUE, start = list(range = 0.1))
dt_est_cold_tall <- p_cold_tall$data
p_cold_tall

nrow(dt[dt$super_biome == "cold_short",])
p_cold_short <- biome_gls(super_b = "cold_short", fit_n = "row_n", part_size = 1500,
                          col_pattern = "mean_evi_", dat = dt, part = FALSE, start = list(range = 0.1))
dt_est_cold_short <- p_cold_short$data
p_cold_short

nrow(dt[dt$super_biome == "not_cold_tall",])
p_not_cold_tall <- biome_gls(super_b = "not_cold_tall", fit_n = "row_n", part_size = 1000,
                             col_pattern = "mean_evi_", dat = dt, part = FALSE, start = list(range = 0.01))
dt_est_not_cold_tall <- p_not_cold_tall$data
p_not_cold_tall

nrow(dt[dt$super_biome == "not_cold_short",])
p_not_cold_short <- biome_gls(super_b = "not_cold_short", fit_n = "row_n", part_size = 1500,
                              col_pattern = "mean_evi_", dat = dt, part = FALSE, start = list(range = 0.1))
dt_est_not_cold_short <- p_not_cold_short$data
p_not_cold_short

p_mean_evi_super_b <- gridExtra::grid.arrange(p_cold_tall, p_cold_short, p_not_cold_tall, p_not_cold_short,
                                                ncol = 4)
ggsave(plot = p_mean_evi_super_b, "builds/plots/mean_evi_remotePARTS_super_b_pas.png", dpi = 600, height = 3, width = 12)

p_mean_evi <- gridExtra::grid.arrange(p_mean_evi_a, p_mean_evi_super_b, ncol = 1)
ggsave(plot = p_mean_evi, "builds/plots/mean_evi_remotePARTS_full_pas.png", dpi = 600, height = 5, width = 14)


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

fwrite(dt_est, "builds/model_estimates/mean_evi_pas_estimates.csv")

