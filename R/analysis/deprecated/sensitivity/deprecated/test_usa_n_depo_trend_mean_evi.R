### test for different biomes separately 
source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/dataFragments/grid_usa_with_climate_trends.csv") %>% 
  distinct(lon, lat, .keep_all = TRUE)


dt <- dt_raw %>% 
  mutate(
    pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(human_modification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA),
    n_depo_usa_coef = scale(n_depo_usa_coef),
    mean_n_depo_zhu = scale(mean_n_depo_zhu),
    n_depo_zhu_coef = scale(n_depo_zhu_coef))  %>% 
  group_by(functional_biome) %>% 
  mutate(n_fb = n()) %>% 
  ungroup() %>% 
  filter(n_fb > 10000)

table(dt$functional_biome)
#write function 

biome_gls <- function(f_biome = NA, col_pattern = NA, start = list(range = 0.1),
                      dat = NA, part = TRUE, part_size = NA, fit_n = "row_n"){
  
  coi <- grep(col_pattern, names(dat), value = T)
  
  dat <- dat %>% filter(functional_biome %in% f_biome) %>%
    dplyr::select(all_of(coi), olson_biome, X, Y, functional_biome, lon, lat,  n_depo_usa_coef,
                  mat_coef, map_coef, max_temp_coef, human_modification) %>% 
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
    
    range_opt_biome <- NULL
    
    # First attempt with the original start value
    if (is.null(range_opt_biome)) {
      tryCatch({
        corfit_biome <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp",
                               start = start, fit.n = fit_n
        )
        range_opt_biome <- corfit_biome$spcor
        
        message("Fit successful on first attempt with start = ", start)
      }, error = function(e) {
        message("Error on first attempt with start = ", start, ": ", conditionMessage(e))
      })
    }
    
    # Second attempt with start = 0.01
    if (is.null(range_opt_biome)) {
      tryCatch({
        start <- list(range = 0.01)
        corfit_biome <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp",
                               start = start, fit.n = fit_n
        )
        range_opt_biome <- corfit_biome$spcor
        
        message("Fit successful on second attempt with start = 0.01")
      }, error = function(e) {
        message("Error on second attempt with start = 0.01: ", conditionMessage(e))
      })
    }
    
    # Third attempt with start = 0.5
    if (is.null(range_opt_biome)) {
      tryCatch({
        start <- list(range = 0.5)
        corfit_biome <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp",
                               start = start, fit.n = fit_n
        )
        range_opt_biome <- corfit_biome$spcor
        
        message("Fit successful on third attempt with start = 0.5")
      }, error = function(e) {
        stop("All attempts failed. Last error with start = 0.5: ", conditionMessage(e))
      })
    }
    
    print(paste0("Range opt: ", range_opt_biome))
  }
  
  
  print(paste0("estimated optimal range (", round(range_opt_biome, 4),"), start fitting GLS"))
  
  if(part == FALSE){
    
    # get distance matrix 
    d_biome <- distm_scaled(coords_biome)
    
    #use r to calculate optimal v parameter 
    v_opt_biome <- covar_exp(d_biome, range_opt_biome)
    
    
    set.seed(161)
    gls_biome <- fitGLS(mean_evi_coef ~ 1 +
                          n_depo_usa_coef +
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
                                    n_depo_usa_coef +
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
  
  
  label <- paste0("Functional Biome:\n", f_biome)
  
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

unique(dt$functional_biome)

nrow(dt[dt$functional_biome == "SLN",])

f_biomes <- unique(dt$functional_biome)
dt_est <- data.table()

for(biome in unique(f_biomes)){
  
  p_biome <- biome_gls(f_biome = biome, fit_n = 15000, part_size = 1500,
                       col_pattern = "mean_evi_", dat = dt, part = TRUE, start = list(range = 0.1))
  tmp_est <- p_biome$data
  
  tmp_est <- tmp_est %>% mutate(
    functional_biome = paste0(biome))
  
  dt_est <- rbind(dt_est, tmp_est)
  
  print(p_biome)
  print(paste0(biome, " done"))
}




dt_plot <- dt_est %>% 
  as.data.table() %>%
  mutate(clean_term = case_when(
    grepl("Intercept", term) ~ "Intercept", 
    grepl("n_depo_usa_coef", term) ~ "N Deposition Trend",
    grepl("mat_coef", term) ~ "MAT Trend", 
    grepl("map_coef", term) ~ "MAP Trend",
    grepl("max_temp_coef", term) ~ "Max Temp Trend",
    grepl("human_modification", term) ~ "Human Modification"),
    estimate = estimate/100, 
    ci_lb = ci_lb/100, 
    ci_ub = ci_ub/100, 
    sdt_error = std_error/100, 
    sig_pn = case_when(
      estimate < 0 & p_value < 0.05 ~ "Sig. Negative", 
      estimate > 0 & p_value < 0.05 ~ "Sig. Positive", 
      p_value >= 0.05 ~ "Non-Significant"
    ), 
    response = "mean_evi"
  ) 

fwrite(dt_plot, "builds/model_estimates/usa_mean_evi_functional_biomes_sep_n_depo_trend.csv")
dt_plot <- fread("builds/model_estimates/usa_mean_evi_functional_biomes_sep_n_depo_trend.csv")

p_est_functional_biome <- dt_plot %>%
  mutate(clean_term = factor(clean_term, levels = c(
    "Intercept", 
    "MAP Trend", 
    "Max Temp Trend", 
    "MAT Trend", 
    "Human Modification",
    "N Deposition Trend"))) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~functional_biome, scales = "free_x", ncol = 5) +
  scale_color_manual(values = c("Sig. Negative" = "#9E3C85",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_est_functional_biome

ggsave(plot = p_est_functional_biome, "builds/plots/supplement/usa_mean_evi_in_functional_biomes_n_depo_trend.png", dpi = 600, height = 10, width = 10)

############################################### SUPER BIOME #########################################
dt <- dt_raw %>% 
  mutate(
    pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(human_modification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA),
    n_depo_usa_coef = scale(n_depo_usa_coef), 
    mean_n_depo_zhu = scale(mean_n_depo_zhu),
    n_depo_zhu_coef = scale(n_depo_zhu_coef))  %>% 
  group_by(super_biome) %>% 
  mutate(n_fb = n()) %>% 
  ungroup() %>% 
  filter(n_fb > 10000)

table(dt$super_biome)
#write function 

biome_gls <- function(f_biome = NA, col_pattern = NA, start = list(range = 0.1),
                      dat = NA, part = TRUE, part_size = NA, fit_n = "row_n"){
  
  coi <- grep(col_pattern, names(dat), value = T)
  
  dat <- dat %>% filter(super_biome %in% f_biome) %>%
    dplyr::select(all_of(coi), olson_biome, X, Y, super_biome, lon, lat, n_depo_usa_coef,
                  mat_coef, map_coef, max_temp_coef, human_modification) %>% 
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
    
    range_opt_biome <- NULL
    
    # First attempt with the original start value
    if (is.null(range_opt_biome)) {
      tryCatch({
        corfit_biome <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp",
                               start = start, fit.n = fit_n
        )
        range_opt_biome <- corfit_biome$spcor
        
        message("Fit successful on first attempt with start = ", start)
      }, error = function(e) {
        message("Error on first attempt with start = ", start, ": ", conditionMessage(e))
      })
    }
    
    # Second attempt with start = 0.01
    if (is.null(range_opt_biome)) {
      tryCatch({
        start <- list(range = 0.01)
        corfit_biome <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp",
                               start = start, fit.n = fit_n
        )
        range_opt_biome <- corfit_biome$spcor
        
        message("Fit successful on second attempt with start = 0.01")
      }, error = function(e) {
        message("Error on second attempt with start = 0.01: ", conditionMessage(e))
      })
    }
    
    # Third attempt with start = 0.5
    if (is.null(range_opt_biome)) {
      tryCatch({
        start <- list(range = 0.5)
        corfit_biome <- fitCor(resids = residuals(ar_biome), coords = coords_biome, covar_FUN = "covar_exp",
                               start = start, fit.n = fit_n
        )
        range_opt_biome <- corfit_biome$spcor
        
        message("Fit successful on third attempt with start = 0.5")
      }, error = function(e) {
        stop("All attempts failed. Last error with start = 0.5: ", conditionMessage(e))
      })
    }
    
    print(paste0("Range opt: ", range_opt_biome))
  }
  
  
  print(paste0("estimated optimal range (", round(range_opt_biome, 4),"), start fitting GLS"))
  
  if(part == FALSE){
    
    # get distance matrix 
    d_biome <- distm_scaled(coords_biome)
    
    #use r to calculate optimal v parameter 
    v_opt_biome <- covar_exp(d_biome, range_opt_biome)
    
    
    set.seed(161)
    gls_biome <- fitGLS(mean_evi_coef ~ 1 +
                          n_depo_usa_coef +
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
                                    n_depo_usa_coef +
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
  
  
  label <- paste0("Functional Biome:\n", f_biome)
  
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

unique(dt$super_biome)

nrow(dt[dt$super_biome == "SLN",])

f_biomes <- unique(dt$super_biome)
dt_est <- data.table()

for(biome in unique(f_biomes)){
  
  p_biome <- biome_gls(f_biome = biome, fit_n = 15000, part_size = 1500,
                       col_pattern = "mean_evi_", dat = dt, part = TRUE, start = list(range = 0.1))
  tmp_est <- p_biome$data
  
  tmp_est <- tmp_est %>% mutate(
    super_biome = paste0(biome))
  
  dt_est <- rbind(dt_est, tmp_est)
  
  print(p_biome)
  print(paste0(biome, " done"))
}




dt_plot <- dt_est %>% 
  as.data.table() %>%
  mutate(clean_term = case_when(
    grepl("Intercept", term) ~ "Intercept", 
    grepl("n_depo_usa_coef", term) ~ "N Deposition Trend",
    grepl("mat_coef", term) ~ "MAT Trend", 
    grepl("map_coef", term) ~ "MAP Trend",
    grepl("max_temp_coef", term) ~ "Max Temp Trend",
    grepl("human_modification", term) ~ "Human Modification"),
    estimate = estimate/100, 
    ci_lb = ci_lb/100, 
    ci_ub = ci_ub/100, 
    sdt_error = std_error/100, 
    sig_pn = case_when(
      estimate < 0 & p_value < 0.05 ~ "Sig. Negative", 
      estimate > 0 & p_value < 0.05 ~ "Sig. Positive", 
      p_value >= 0.05 ~ "Non-Significant"
    ), 
    response = "mean_evi"
  ) 

fwrite(dt_plot, "builds/model_estimates/usa_mean_evi_super_biomes_sep_n_depo_trend.csv")
dt_plot <- fread("builds/model_estimates/usa_mean_evi_super_biomes_sep_n_depo_trend.csv")


p_est_super_biome <- dt_plot %>%
  mutate(clean_term = factor(clean_term, levels = c(
    "Intercept", 
    "MAP Trend", 
    "Max Temp Trend", 
    "MAT Trend", 
    "Human Modification",
    "N Deposition Trend"))) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~super_biome, scales = "free_x", ncol = 4) +
  scale_color_manual(values = c("Sig. Negative" = "#9E3C85",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_est_super_biome

ggsave(plot = p_est_super_biome, "builds/plots/supplement/usa_mean_evi_in_super_biomes_n_depo_trend.png", dpi = 600, height = 10, width = 10)


##### Global usa #####

library(data.table)
library(sf)
library(tidyverse)
library(remotePARTS)

dt_raw <- fread("data/processedData/dataFragments/grid_usa_with_climate_trends.csv") %>% 
  distinct(lon, lat, .keep_all = TRUE)

dt <- dt_raw %>% 
  mutate(
    pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(human_modification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA),
    n_depo_usa_coef = scale(n_depo_usa_coef), 
    mean_n_depo_zhu = scale(mean_n_depo_zhu),
    n_depo_zhu_coef = scale(n_depo_zhu_coef)) 

## EVI ------------------------------------------------------------------------------
####### calculate mean_evi trend #######

mean_evi_cols <- grep("mean_evi_", names(dt), value = T)

#subset to complete cases
dt_mean_evi <- dt %>%
  dplyr::select(all_of(mean_evi_cols), functional_biome, X, Y, lon, lat, unique_id, 
                n_depo_usa_coef,
                mat_coef, map_coef, max_temp_coef, human_modification, super_biome,
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


## usa N MEAN -------------

set.seed(161)
gls_mean_evi_1 <- fitGLS_partition(mean_evi_coef ~ 1 +
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
                                   coord.names = c("lon", "lat"), 
                                   debug = F)
gls_mean_evi_1 

dt_est_mean_evi_1 <- extract_gls_estimates(gls_mean_evi_1, part = TRUE)

dt_plot <- dt_est_mean_evi_1 %>% 
  as.data.table() %>%
  mutate(clean_term = case_when(
    grepl("Intercept", term) ~ "Intercept", 
    grepl("n_depo_usa_coef", term) ~ "N Deposition Trend",
    grepl("mat_coef", term) ~ "MAT Trend", 
    grepl("map_coef", term) ~ "MAP Trend",
    grepl("max_temp_coef", term) ~ "Max Temp Trend",
    grepl("human_modification", term) ~ "Human Modification"),
    estimate = estimate/100, 
    ci_lb = ci_lb/100, 
    ci_ub = ci_ub/100, 
    sdt_error = std_error/100, 
    sig_pn = case_when(
      estimate < 0 & p_value < 0.05 ~ "Sig. Negative", 
      estimate > 0 & p_value < 0.05 ~ "Sig. Positive", 
      p_value >= 0.05 ~ "Non-Significant"
    ), 
    response = "mean_evi"
  ) 


p_est_global <- dt_plot %>%
  mutate(clean_term = factor(clean_term, levels = c(
    "Intercept", 
    "MAP Trend", 
    "Max Temp Trend", 
    "MAT Trend", 
    "Human Modification",
    "N Deposition Trend"))) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("Sig. Negative" = "#9E3C85",
                                "Non-Significant" = "grey60", 
                                "Sig. Positive" = "#457B2A"
  )) +
  labs(title = "usa", y = NULL, x = "EVI Trend Estimate", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none", 
        # panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11.5),   
        strip.background = element_rect(fill = "snow2", color = "snow2"),
        panel.background = element_rect(fill = "snow1", color = "snow1"),
        panel.border = element_blank(),
        plot.title = element_text(size = 12))

p_est_global

ggsave(plot = p_est_global, "builds/plots/supplement/usa_mean_evi_n_depo_n_depo_trend.png", dpi = 600)

p_comb <- gridExtra::grid.arrange(p_est_functional_biome, p_est_super_biome, p_est_global)
ggsave(plot = p_comb, "builds/plots/supplement/usa_mean_evi_n_depo_estimates_n_depo_trend.png", dpi = 600, height = 6, width = 7)

