### test for different biomes separately 
source("R/functions/extract_gls_estimates.R")


library(tidyverse)
library(data.table)
library(remotePARTS)


dt_raw <- fread("data/processedData/dataFragments/grid_sample_with_climate_trends.csv") %>% 
  as.data.frame()

dt_raw <- dt_raw %>% 
  mutate(
    pa_age = ifelse(STATUS_YR > 1800, 2023-STATUS_YR, NA), 
    nitrogen_depo = scale(nitrogen_depo),
    mat_coef = scale(mat_coef),
    max_temp_coef = scale(max_temp_coef),
    map_coef = scale(map_coef),
    human_modification = scale(human_modification), 
    area_km2_log = scale(log(area_km2 + 0.0001)),
    pa_age_log = ifelse(!pa_age == 0, scale(log(pa_age + 0.0001)), NA)) %>% 
  group_by(olson_biome) %>% 
  mutate(n_ob = n()) %>% 
  ungroup() %>% 
  filter(n_ob > 10000)

table(dt_raw$olson_biome)


dt <- dt_raw

#write function 

biome_gls <- function(o_biome = NA, col_pattern = NA, start = list(range = 0.1),
                      dat = NA, part = TRUE, part_size = NA, fit_n = "row_n"){
  
  coi <- grep(col_pattern, names(dat), value = T)
  
  dat <- dat %>% filter(olson_biome %in% o_biome) %>%
    dplyr::select(all_of(coi), olson_biome, X, Y, functional_biome, lon, lat, 
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
                                  ncores = 25,
                                  progressbar = TRUE, 
                                  parallel = TRUE, 
                                  coord.names = c("lon", "lat"))
    gls_biome 
    
    dt_est_biome <- extract_gls_estimates(gls_biome, part = TRUE)
  }
  
  
  label <- paste0("Olson Biome:\n", o_biome)
  
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

unique(dt$olson_biome)


o_biomes <- unique(dt$olson_biome)
dt_est <- data.table()

for(biome in unique(o_biomes)){
  
  p_biome <- biome_gls(o_biome = biome, fit_n = 15000, part_size = 1500,
                       col_pattern = "mean_evi_", dat = dt, part = TRUE, start = list(range = 0.1))
  tmp_est <- p_biome$data
  
  tmp_est <- tmp_est %>% mutate(
    olson_biome = paste0(biome))
  
  dt_est <- rbind(dt_est, tmp_est)
  
  print(p_biome)
  print(paste0(biome, " done"))
}




dt_plot <- dt_est %>% 
  as.data.table() %>%
  mutate(clean_term = case_when(
    grepl("Intercept", term) ~ "Intercept", 
    grepl("nitrogen_depo", term) ~ "Nitrogen Deposition",
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

fwrite(dt_plot, "builds/model_estimates/mean_evi_olson_biomes_sep.csv")


p_est <- dt_plot %>%
  mutate(clean_term = factor(clean_term, levels = c(
    "Intercept", 
    "MAP Trend", 
    "Max Temp Trend", 
    "MAT Trend", 
    "Human Modification",
    "Nitrogen Deposition")),
    olson_biome = gsub("&", "&\n", olson_biome), 
    olson_biome = gsub(",", ",\n", olson_biome),
    olson_biome = gsub("Moist", "Moist\n", olson_biome)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(x = estimate, xmin = ci_lb, xmax = ci_ub, y = clean_term, color = sig_pn),
                  alpha = 0.9, linewidth = 1.2) +
  facet_wrap(~olson_biome, scales = "free", ncol = 3) +
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

p_est

ggsave(plot = p_est, "builds/plots/supplement/mean_evi_in_olson_biomes.png", dpi = 600, height = 10, width = 10)

