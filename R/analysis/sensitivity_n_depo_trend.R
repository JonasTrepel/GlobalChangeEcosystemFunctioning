library(tidyverse)
library(data.table)
library(tidylog)
library(ggcorrplot)
library(broom)
library("sdmTMB")
library(future)
library(furrr)
library(groupdata2)
library(GGally)
library(glmmTMB)

#### Leave one site out CV


#1 HOUSEKEEPING -------------------------------------

#load data 
#sf_parks <- st_read("data/spatial_data/protected_areas/park_boundaries.gpkg") 

dt <- fread("data/processed_data/analysis_ready_n_depo_trend_grid.csv")

setDT(dt)


# get dataframe with comlete and clean data fro modeling 

acceptable_numbers = seq(1, 10000000, 5)
# mutate(park_row_nr = 1:n()) %>% 
#  filter(park_row_nr %in% acceptable_numbers) %>% 
dt_mod <- dt %>% 
  mutate(
    x_moll_km = x_moll/1000, 
    y_moll_km = y_moll/1000)


subsets <- c("Europe", "USA") 

responses <- c("evi_coef")

mesh_guide_raw <- CJ(subset = subsets, 
                     response = responses) %>% 
  filter(subset != "") %>% 
  mutate(subset_col = "region",
  tier = subset, 
  filter_call = paste0(subset_col, " == '", subset,"'"), 
  divisor = as.numeric(1))

#get divisor to keem sample size on check 

for(i in 1:nrow(mesh_guide_raw)){
  
  (filter_call = mesh_guide_raw[i, ]$filter_call)
  dt_sub <- dt_mod %>% 
    dplyr::filter(eval(parse(text = filter_call)))
  
  if(nrow(dt_sub) < 10000){
    
    (divisor = NA)
    
  } else if(nrow(dt_sub) <= 100000) {
    
    (divisor = 1)
    
  }else if(nrow(dt_sub) > 100000) {
    
    x = nrow(dt_sub)/100000
    (divisor = floor(x))
    
  }
  
  mesh_guide_raw[i, ]$divisor <- divisor
  
}

mesh_guide_raw <- mesh_guide_raw %>% 
  mutate(divisor = ifelse(subset == "yes", 10, divisor)) %>% 
  filter(!is.na(divisor))

### 2 - Choose Mesh ------------------
#https://www.biorxiv.org/content/10.1101/2022.03.24.485545v4.full.pdf

#https://becarioprecario.bitbucket.io/spde-gitbook/ch-intro.html on how to construct meshs , chapter 2.6 and 2.7

#the following may take a good couple of hours to finish 
mesh_grid <- expand.grid(max_inner_edge = 1000, cutoff = c(20, 40, 80, 160, 320), 
                         tier = unique(mesh_guide_raw$tier)) %>% 
  group_by(tier) %>% 
  mutate(mesh_id = paste0("mesh_", 1:n())) %>% 
  ungroup() %>% 
  unique()


mesh_guide = mesh_guide_raw %>% 
  unique() %>% 
  left_join(mesh_grid)

plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GiB
start_time <- Sys.time()


mesh_res_list <- list()
dt_mesh_res <- data.frame()



all_mesh_results <- future_map(
  1:nrow(mesh_guide),
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE),
  function(j) {
    
    tier  <- mesh_guide[j, ]$tier
    resp <- mesh_guide[j, ]$response
    filter_call <- mesh_guide[j, ]$filter_call
    divisor <- mesh_guide[j, ]$divisor
    cutoff <- mesh_guide[j, ]$cutoff
    max_inner <- mesh_guide[j, ]$max_inner_edge
    mesh_id <- mesh_guide[j, ]$mesh_id
    
    message("Starting tier: ", tier, " at ", Sys.time())
    
    acceptable_numbers <- seq(1, 10000000, by = divisor)
    
    dt_sub <- dt_mod %>%
      filter(eval(parse(text = filter_call))) %>%
      mutate(row_nr = 1:n()) %>%
      filter(row_nr %in% acceptable_numbers) %>%
      mutate(
        n_depo_coef_scaled = as.numeric(scale(n_depo_coef)),
        fire_frequency_scaled = as.numeric(scale(fire_frequency)),
        mat_coef_scaled = as.numeric(scale(mat_coef)),
        prec_coef_scaled = as.numeric(scale(prec_coef)),
        hmi_change_scaled = as.numeric(scale(hmi_change)),
        max_temp_coef_scaled = as.numeric(scale(max_temp_coef)), 
        mean_prec_scaled = as.numeric(scale(mean_prec)), 
        mean_mat_scaled = as.numeric(scale(mean_mat))
      )
    
    inla_mesh <- fmesher::fm_mesh_2d_inla(
      loc = cbind(dt_sub$x_moll_km, dt_sub$y_moll_km),
      cutoff = cutoff,
      max.edge = c(max_inner, 100000)
    )
    
    mesh <- make_mesh(
      data = dt_sub,
      xy_cols = c("x_moll_km", "y_moll_km"),
      mesh = inla_mesh
    )
    
    formula <- as.formula(
      paste0(resp, " ~ n_depo_coef_scaled +
                     mat_coef_scaled +
                     prec_coef_scaled +
                     hmi_change_scaled +
                     max_temp_coef_scaled +
                     fire_frequency_scaled +
                     mean_mat_scaled +
                     mean_prec_scaled"))
    
    fit_cv <- tryCatch({
      sdmTMB::sdmTMB_cv(
        formula,
        data = dt_sub,
        mesh = mesh,
        # family = sdmTMB::student(),
        spatial = "on",
        reml = T
      )
    }, error = function(e) {
      message("Skipping CV model: ", e$message)
      return(NULL)
    })
    
    if (is.null(fit_cv)) return(NULL)
    
    tier_ = gsub(" ", "_", tier)
    tier_ = gsub("/", "_", tier_)
    
    cv_model_id <- paste0("cv_", tier_, "_", mesh_id, "_n_depo_trend")
    
    result_row <- data.frame(
      tier = tier,
      filter_call = filter_call,
      cutoff = cutoff,
      max_inner_edge = max_inner,
      mesh_id  = mesh_id,
      n_vertices  = nrow(mesh$mesh$loc),
      all_converged = fit_cv$converged,
      p_d_hessian = sum(fit_cv$pdHess),
      response = resp,
      divisor = divisor,
      sum_loglik = fit_cv$sum_loglik,
      n  = nrow(dt_sub),
      log_cpo_approx = fit_cv$sum_loglik / nrow(dt_sub),
      model_id = cv_model_id,
      model_path = paste0("builds/cv_models/", cv_model_id, ".Rds")
    )
    
    saveRDS(fit_cv, file = result_row$model_path)
    
    rm(fit_cv)
    gc()
    
    message("Finished tier: ", tier, " at ", Sys.time())
    
    return(result_row)
  }
) 

dt_mesh_res <- rbindlist(all_mesh_results)

plan(sequential)

print(paste0("Started loop at: ", start_time, " and finished at: ", Sys.time()))
print(paste0("Estimate time for CV: ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), " mins"))

## bind results 
unique(responses)
dt_mesh_res_fin <- dt_mesh_res  %>% 
  mutate(clean_response = case_when(
    .default = response,
    response == "evi_coef" ~ "EVI Trend"))
summary(dt_mesh_res_fin)

fwrite(dt_mesh_res_fin, "builds/model_outputs/cv_mesh_selection_sdmtmb_n_depo_trend.csv")
n_distinct(dt_mesh_res_fin$model_id)


p_loglik <- dt_mesh_res_fin %>% 
  filter(tier == "full_dataset_yes") %>% 
  ggplot() +
  geom_point(aes(x = cutoff, y = sum_loglik), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "B", direction = - 1, begin = 0.2, end = 0.8) +
  labs(y = "Log Likelihood (sum)", x = "Cutoff (km)", color = "Max\nInner\nEdge\n(km)") +
  # facet_wrap(~clean_response, scales = "free") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = .5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "snow"), 
        strip.background = element_rect(fill = "linen", color = "linen"))
p_loglik
ggsave(plot = p_loglik, "builds/plots/supplement/sum_loglik_different_meshs_n_depo_trend.png", dpi = 600, height = 4, width = 8)


##### run models --------------------------------------------

dt_mesh_res_fin <- fread("builds/model_outputs/cv_mesh_selection_sdmtmb_n_depo_trend.csv")

dt_best_mesh <- dt_mesh_res_fin %>% 
  group_by(tier) %>% 
  slice_max(sum_loglik) %>% 
  ungroup()



### - Use best  Mesh ------------------

plan(multisession, workers = 42)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GiB
start_time <- Sys.time()

best_mesh_res_list <- future_map(1:nrow(dt_best_mesh),
                                 .progress = TRUE,
                                 .options = furrr_options(seed = TRUE),
                                 function(i) {
                                   
                                   tier = dt_best_mesh[i,]$tier
                                   resp = dt_best_mesh[i,]$response 
                                   filter_call = dt_best_mesh[i,]$filter_call
                                   divisor = dt_best_mesh[i, ]$divisor
                                   
                                   
                                   acceptable_numbers = seq(1, 10000000, divisor)
                                   
                                   dt_sub <- dt_mod %>% 
                                     dplyr::filter(eval(parse(text = filter_call))) %>% 
                                     mutate(row_nr = 1:n()) %>% 
                                     filter(row_nr %in% acceptable_numbers) %>% 
                                     mutate(
                                       n_depo_coef_scaled = as.numeric(scale(n_depo_coef)),
                                       fire_frequency_scaled = as.numeric(scale(fire_frequency)),
                                       mat_coef_scaled = as.numeric(scale(mat_coef)),
                                       prec_coef_scaled = as.numeric(scale(prec_coef)),
                                       hmi_change_scaled = as.numeric(scale(hmi_change)), 
                                       max_temp_coef_scaled = as.numeric(scale(max_temp_coef)), 
                                       mean_mat_scaled = as.numeric(scale(mean_mat)), 
                                       mean_prec_scaled = as.numeric(scale(mean_prec)))
                                   
                                   co <- as.numeric(dt_best_mesh[i, ]$cutoff)
                                   i_e <- as.numeric(dt_best_mesh[i, ]$max_inner_edge)
                                   
                                   inla_mesh <- fmesher::fm_mesh_2d_inla(
                                     loc = cbind(dt_sub$x_moll_km, dt_sub$y_moll_km),
                                     cutoff = co,
                                     max.edge = c(i_e, 100000)
                                   )
                                   
                                   mesh <- make_mesh(
                                     data = dt_sub,
                                     xy_cols = c("x_moll_km", "y_moll_km"),
                                     mesh = inla_mesh
                                   )
                                   
                                   
                                   int_formula <- as.formula(paste0(resp, "~ 1"))
                                   
                                   n_dep_formula <- as.formula(paste0(resp, "~ n_depo_coef_scaled"))
                                   
                                   no_dep_formula <- as.formula(paste0(resp, "~
                                   prec_coef_scaled +
                                   hmi_change_scaled + 
                                   mat_coef_scaled +
                                   max_temp_coef_scaled +
                                   fire_frequency_scaled +
                                   mean_mat_scaled +
                                   mean_prec_scaled"))
                                   
                                   
                                   fixed_formula <- as.formula(paste0(resp, " ~
                                   n_depo_coef_scaled +
                                   prec_coef_scaled +
                                   mat_coef_scaled +
                                   hmi_change_scaled + 
                                   max_temp_coef_scaled +
                                   fire_frequency_scaled +
                                   mean_mat_scaled +
                                   mean_prec_scaled"))
                                   
                                   full_formula <- as.formula(paste0(resp, " ~
                                   n_depo_coef_scaled +
                                   prec_coef_scaled +
                                   hmi_change_scaled +
                                   mat_coef_scaled +
                                   max_temp_coef_scaled +
                                   fire_frequency_scaled +
                                   mean_mat_scaled +
                                   mean_prec_scaled"))
                                   
                                   #https://github.com/pbs-assess/sdmTMB/issues/466#issuecomment-3119589818
                                   
                                   # intercept only
                                   fit_int <- tryCatch(
                                     sdmTMB(int_formula,
                                            spatial = "off",
                                            data = dt_sub,
                                            mesh = mesh,
                                            reml = TRUE),
                                     error = function(e) NULL
                                   )
                                   if (is.null(fit_int)){return(NULL)}
                                   
                                   # intercept + spatial
                                   fit_sp <- tryCatch(
                                     update(fit_int,
                                            spatial = "on",
                                            reml = TRUE),
                                     error = function(e) NULL
                                   )
                                   if (is.null(fit_sp)){return(NULL)}
                                   
                                   # full model
                                   fit_full <- tryCatch(
                                     sdmTMB(formula = full_formula,
                                            mesh = mesh, 
                                            data = dt_sub,
                                            spatial = "on",
                                            reml = TRUE),
                                     error = function(e) NULL
                                   )
                                   if (is.null(fit_full)){return(NULL)}
                                   
                                   # fixed effects only
                                   fit_fixed <- tryCatch(
                                     sdmTMB(formula = fixed_formula,
                                            mesh = mesh, 
                                            data = dt_sub,
                                            spatial = "off",
                                            reml = TRUE),
                                     error = function(e) NULL
                                   )
                                   if (is.null(fit_fixed)){return(NULL)}
                                   
                                   # n deposition only
                                   fit_n_dep <- tryCatch(
                                     sdmTMB(formula = n_dep_formula,
                                            mesh = mesh, 
                                            data = dt_sub,
                                            spatial = "off",
                                            reml = TRUE),
                                     error = function(e) NULL
                                   )
                                   if (is.null(fit_n_dep)){return(NULL)}
                                   
                                   # no deposition
                                   fit_no_dep <- tryCatch(
                                     sdmTMB(formula = no_dep_formula,
                                            mesh = mesh, 
                                            data = dt_sub,
                                            spatial = "off",
                                            reml = TRUE),
                                     error = function(e) NULL
                                   )
                                   if (is.null(fit_no_dep)){return(NULL)}
                                   
                                   # total proportion deviance explained by our full model:
                                   (dev_explained_full <- 1 - deviance(fit_full) / deviance(fit_int))
                                   
                                   # proportion deviance explained by the mesh:
                                   (dev_explained_spatial <- 1 - deviance(fit_sp) / deviance(fit_int))
                                   
                                   # proportion deviance explained by the n deposition only:
                                   (dev_explained_fixed <- 1 - deviance(fit_fixed) / deviance(fit_int))
                                   
                                   # proportion deviance explained by the covariate:
                                   (dev_explained_n_dep <- 1 - deviance(fit_n_dep) / deviance(fit_int))
                                   
                                   # proportion deviance explained by covariates but no n dep:
                                   (dev_explained_no_dep <- 1 - deviance(fit_no_dep) / deviance(fit_int))
                                   
                                   
                                   san <- sdmTMB::sanity(fit_full)
                                   
                                   tier_ = gsub(" ", "_", tier)
                                   tier_ = gsub("/", "_", tier_)
                                   model_id = paste0("subset_",tier_, "_n_depo_trend")
                                   
                                   tmp_tidy <- broom::tidy(fit_full, conf.int = TRUE) %>%
                                     #dplyr::filter(!grepl("(Intercept)", term)) %>%
                                     dplyr::mutate(sig = case_when(
                                       .default = "non-significant",
                                       conf.low > 0 ~ "positive",
                                       conf.high < 0 ~ "negative"
                                     )) %>%
                                     dplyr::mutate(
                                       response = resp, 
                                       tier = tier,
                                       dev_explained_fixed = dev_explained_fixed, 
                                       dev_explained_full = dev_explained_full,
                                       dev_explained_n_dep = dev_explained_n_dep, 
                                       dev_explained_no_dep = dev_explained_no_dep,
                                       dev_explained_spatial = dev_explained_spatial,
                                       cutoff = co,
                                       max_inner_edge = i_e,
                                       n_vertices = nrow(mesh$mesh$loc),
                                       aic = AIC(fit_full),
                                       n = nrow(dt_sub), 
                                       model_id = model_id, 
                                       model_path = paste0("builds/models/", model_id, ".Rds"),
                                       hessian_ok = san$hessian_ok, 
                                       eigen_values_ok = san$eigen_values_ok,
                                       nlminb_ok = san$nlminb_ok,
                                       range_ok = san$range_ok,
                                       gradients_ok = san$gradients_ok,
                                       se_magnitude_ok = san$se_magnitude_ok,
                                       se_na_ok = san$se_na_ok,
                                       sigmas_ok = san$sigmas_ok,
                                       all_ok = san$all_ok
                                       
                                     )
                                   
                                   saveRDS(fit_full, file = unique(tmp_tidy$model_path))
                                   
                                   print(paste0(i, " done"))
                                   
                                   rm(fit_int, fit_no_dep, fit_sp, fit_n_dep, fit_fixed, fit_full)
                                   gc()
                                   
                                   return(tmp_tidy)
                                   
                                 })
plan(sequential)

print(paste0("Started loop at: ", start_time, " and finished at: ", Sys.time()))

## bind results 
unique(responses)
dt_res <- rbindlist(best_mesh_res_list) %>% 
  mutate(clean_response = case_when(
    .default = response,
    response == "evi_coef" ~ "EVI Trend",
  ), 
  clean_term = case_when(
    .default = term,
    term == "n_depo_coef_scaled" ~ "N Deposition Trend",
    term == "mat_coef_scaled" ~ "MAT Trend",
    term == "prec_coef_scaled" ~ "Precipitation Trend",
    term == "max_temp_coef_scaled" ~ "Max. Temperature Trend",
    term == "hmi_change_scaled" ~ "HMI Change",
    term == "mean_mat_scaled" ~ "MAT", 
    term == "mean_prec_scaled" ~ "MAP", 
    term == "fire_frequency_scaled" ~ "Fire frequency"))
unique(dt_res$clean_term)
summary(dt_res)
fwrite(dt_res, "builds/model_outputs/sdmtmb_results_n_depo_trend.csv")

############## PLOTS ######################################
library(data.table)
library(tidyverse)
library(sf)
library(scico)
dt_res <- fread("builds/model_outputs/sdmtmb_results_n_depo_trend.csv")

p_est <- dt_res %>% 
  mutate(sig = (conf.low > 0 | conf.high < 0)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = clean_term, 
                      x = estimate, xmin = conf.low, xmax = conf.high, 
                      color = sig), 
                  shape = 18, linewidth = 1, size = 1.3, alpha = .85) +
  scale_color_manual(values = c("#D7E7C0", "#0C4C00")) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~tier) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = "white"), 
    strip.background = element_rect(fill = "snow", color = "snow")
  ) +
  labs(x = "Estimate",y = NULL)

p_est
ggsave(plot = p_est, "builds/plots/supplement/usa_europe_n_depo_trend_estimates.png", 
       dpi = 900, height = 3, width = 7)


m_us = readRDS(unique(dt_res[dt_res$tier == "USA", ]$model_path))
dt_us <- m_us$data

m_eu = readRDS(unique(dt_res[dt_res$tier == "Europe", ]$model_path))
dt_eu <- m_eu$data


p_us_mean = dt_mod %>% 
  filter(region == "USA") %>%
  ggplot() +
  geom_tile(aes(x = x_moll, y = y_moll, fill = nitrogen_depo)) +
  scale_fill_scico(palette = "batlow", begin = .1, end = .9) +
  theme_void() +
  coord_fixed() +
  theme(legend.position = "bottom") +
  labs(fill = "Mean N deposition")
p_us_mean


p_us_trend = dt_mod %>% 
  filter(region == "USA") %>%
  mutate(n_depo_coef = case_when(
    .default = n_depo_coef, 
    n_depo_coef < quantile(n_depo_coef, .025) ~ quantile(n_depo_coef, .025), 
    n_depo_coef > quantile(n_depo_coef, .975) ~ quantile(n_depo_coef, .975)
  )) %>% 
  ggplot() +
  geom_tile(aes(x = x_moll, y = y_moll, fill = n_depo_coef)) +
  scale_fill_scico(palette = "roma",midpoint = 0, direction = -1) +
  theme_void() +
  coord_fixed() +
  theme(legend.position = "bottom") +
  labs(fill = "N deposition Trend")
p_us_trend


p_eu_mean = dt_mod %>% 
  filter(region == "Europe") %>%
  ggplot() +
  geom_tile(aes(x = x_moll, y = y_moll, fill = nitrogen_depo)) +
  scale_fill_scico(palette = "batlow", begin = .1, end = .9, 
                   breaks = scales::pretty_breaks(n = 3)) +
  theme_void() +
  coord_fixed() +
  theme(legend.position = "bottom") +
  labs(fill = "Mean N deposition")
p_eu_mean


p_eu_trend = dt_mod %>% 
  filter(region == "Europe") %>%
  mutate(n_depo_coef = case_when(
    .default = n_depo_coef, 
    n_depo_coef < quantile(n_depo_coef, .025) ~ quantile(n_depo_coef, .025), 
    n_depo_coef > quantile(n_depo_coef, .975) ~ quantile(n_depo_coef, .975)
  )) %>% 
  ggplot() +
  geom_tile(aes(x = x_moll, y = y_moll, fill = n_depo_coef)) +
  scale_fill_scico(palette = "roma",midpoint = 0, direction = -1) +
  theme_void() +
  coord_fixed() +
  theme(legend.position = "bottom") +
  labs(fill = "N deposition Trend")
p_eu_trend


library(patchwork)

p_maps = ((p_eu_mean | p_eu_trend) / (p_us_mean | p_us_trend))
p_all = (p_maps / p_est) + plot_layout(ncol = 1, heights = c(1.25, 1.25, 1)) + plot_annotation(tag_levels = "A")
ggsave(plot = p_all, "builds/plots/supplement/n_depo_trend_full.png", dpi = 900, 
       height = 10, width = 10)
