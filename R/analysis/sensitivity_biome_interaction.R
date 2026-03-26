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



#1 HOUSEKEEPING -------------------------------------

#load data 
#sf_parks <- st_read("data/spatial_data/protected_areas/park_boundaries.gpkg") 

dt <- fread("data/processed_data/analysis_ready_world_grid.csv")

setDT(dt)


# get dataframe with complete and clean data for modeling 

acceptable_numbers = seq(1, 10000000, 3)
# mutate(park_row_nr = 1:n()) %>% 
#  filter(park_row_nr %in% acceptable_numbers) %>% 
dt_mod <- dt %>% 
  mutate(
    x_moll_km = x_moll/1000, 
    y_moll_km = y_moll/1000, 
    full_dataset = "yes")

#check corr

dt_corr <- dt_mod %>%
  dplyr::select(nitrogen_depo, hmi_change,
                mat_coef, prec_coef, max_temp_coef, 
                fire_frequency, evi_coef, 
                mean_mat, mean_prec)
corr <- round(cor(dt_corr), 1)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE, colors = c("#6D9EC1", "white", "#E46726"))


tiers <- c("functional_biome", 
             "climatic_region", 
             "olson_biome",
             "super_biome") 

responses <- c("evi_coef")

mesh_guide_raw <- CJ(tier = tiers, 
                     response = responses) %>% 
  filter(tier != "") %>% 
  mutate(divisor = as.numeric(10))



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

plan(multisession, workers = 20)
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
    divisor <- mesh_guide[j, ]$divisor
    cutoff <- mesh_guide[j, ]$cutoff
    max_inner <- mesh_guide[j, ]$max_inner_edge
    mesh_id <- mesh_guide[j, ]$mesh_id
    
    message("Starting tier: ", tier, " at ", Sys.time())
    
    acceptable_numbers <- seq(1, 10000000, by = divisor)
    
    dt_sub <- dt_mod %>%
      mutate(row_nr = 1:n()) %>%
      filter(row_nr %in% acceptable_numbers) %>%
      mutate(
        nitrogen_depo_scaled = as.numeric(scale(nitrogen_depo)),
        fire_frequency_scaled = as.numeric(scale(fire_frequency)),
        mat_coef_scaled = as.numeric(scale(mat_coef)),
        prec_coef_scaled = as.numeric(scale(prec_coef)),
        hmi_change_scaled = as.numeric(scale(hmi_change)),
        max_temp_coef_scaled = as.numeric(scale(max_temp_coef)),
        mean_mat_scaled = as.numeric(scale(mean_mat)),
        mean_prec_scaled = as.numeric(scale(mean_prec))
      )
    
    
    if(tier == "functional_biome"){
      dt_sub <- dt_sub %>% 
        mutate(biome_col = functional_biome)
    }else if(tier == "climatic_region"){
      dt_sub <- dt_sub %>% 
        mutate(biome_col = climatic_region)
    }else if(tier == "olson_biome"){
      dt_sub <- dt_sub %>% 
        mutate(biome_col = olson_biome)
    }else if(tier == "super_biome"){
      dt_sub <- dt_sub %>% 
        mutate(biome_col = super_biome)
    }
    
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
      paste0(resp, " ~ biome_col*nitrogen_depo_scaled +
                     biome_col*mat_coef_scaled +
                     biome_col*prec_coef_scaled +
                     biome_col*hmi_change_scaled +
                     biome_col*max_temp_coef_scaled +
                     biome_col*fire_frequency_scaled +
                     biome_col*mean_mat_scaled +
                     biome_col*mean_prec_scaled"))
    
    fit_cv <- tryCatch({
      sdmTMB::sdmTMB_cv(
        formula,
        data = dt_sub,
        mesh = mesh,
        # family = sdmTMB::student(),
        spatial = "on",
        k_folds = 5,
        reml = T
      )
    }, error = function(e) {
      message("Skipping CV model: ", e$message)
      return(NULL)
    })
    
    if (is.null(fit_cv)) return(NULL)
    
    tier_ = gsub(" ", "_", tier)
    tier_ = gsub("/", "_", tier_)
    
    cv_model_id <- paste0("cv_", tier_, "_", mesh_id, "_biome_interaction_sensitivity")
    
    result_row <- data.frame(
      tier = tier,
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

fwrite(dt_mesh_res_fin, "builds/model_outputs/cv_mesh_selection_sdmtmb_biome_interaction_sensitivity.csv")
n_distinct(dt_mesh_res_fin$model_id)

##### run models --------------------------------------------

dt_mesh_res_fin <- fread("builds/model_outputs/cv_mesh_selection_sdmtmb_biome_interaction_sensitivity.csv")

dt_best_mesh <- dt_mesh_res_fin %>% 
  group_by(tier) %>% 
  slice_max(sum_loglik) %>% 
  ungroup()



### - Use best  Mesh ------------------

plan(multisession, workers = 4)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GiB
start_time <- Sys.time()

best_mesh_res_list <- future_map(1:nrow(dt_best_mesh),
                                 .progress = TRUE,
                                 .options = furrr_options(seed = TRUE),
                                 function(i) {
                                   
                                   tier = dt_best_mesh[i,]$tier
                                   resp = dt_best_mesh[i,]$response 
                                   divisor = dt_best_mesh[i, ]$divisor
                                   
                                   
                                   acceptable_numbers = seq(1, 10000000, divisor)
                                   
                                   dt_sub <- dt_mod %>% 
                                     mutate(row_nr = 1:n()) %>% 
                                     filter(row_nr %in% acceptable_numbers) %>% 
                                     mutate(
                                       nitrogen_depo_scaled = as.numeric(scale(nitrogen_depo)),
                                       fire_frequency_scaled = as.numeric(scale(fire_frequency)),
                                       mat_coef_scaled = as.numeric(scale(mat_coef)),
                                       prec_coef_scaled = as.numeric(scale(prec_coef)),
                                       hmi_change_scaled = as.numeric(scale(hmi_change)), 
                                       max_temp_coef_scaled = as.numeric(scale(max_temp_coef)),
                                       mean_mat_scaled = as.numeric(scale(mean_mat)),
                                       mean_prec_scaled = as.numeric(scale(mean_prec)))
                                   
                                   if(tier == "functional_biome"){
                                     dt_sub <- dt_sub %>% 
                                       mutate(biome_col = functional_biome)
                                   }else if(tier == "climatic_region"){
                                     dt_sub <- dt_sub %>% 
                                       mutate(biome_col = climatic_region)
                                   }else if(tier == "olson_biome"){
                                     dt_sub <- dt_sub %>% 
                                       mutate(biome_col = olson_biome)
                                   }else if(tier == "super_biome"){
                                     dt_sub <- dt_sub %>% 
                                       mutate(biome_col = super_biome)
                                   }
                                   
                                   
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
                                   
                                   n_dep_formula <- as.formula(paste0(resp, "~ biome_col*nitrogen_depo_scaled"))
                                   
                                   no_dep_formula <- as.formula(paste0(resp, "~
                                   biome_col*prec_coef_scaled +
                                   biome_col*hmi_change_scaled + 
                                   biome_col*mat_coef_scaled +
                                   biome_col*max_temp_coef_scaled +
                                   biome_col*fire_frequency_scaled +
                                   biome_col*mean_mat_scaled +
                                   biome_col*mean_prec_scaled"))
                                   
                                   
                                   fixed_formula <- as.formula(paste0(resp, " ~
                                   biome_col*nitrogen_depo_scaled +
                                   biome_col*prec_coef_scaled +
                                   biome_col*mat_coef_scaled +
                                   biome_col*hmi_change_scaled + 
                                   biome_col*max_temp_coef_scaled +
                                   biome_col*fire_frequency_scaled +
                                   biome_col*mean_mat_scaled +
                                   biome_col*mean_prec_scaled"))
                                   
                                   full_formula <- as.formula(paste0(resp, " ~
                                   biome_col*nitrogen_depo_scaled +
                                   biome_col*prec_coef_scaled +
                                   biome_col*hmi_change_scaled +
                                   biome_col*mat_coef_scaled +
                                   biome_col*max_temp_coef_scaled +
                                   biome_col*fire_frequency_scaled +
                                   biome_col*mean_mat_scaled +
                                   biome_col*mean_prec_scaled"))
                                   
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
                                   
                                   model_id = paste0("subset_",tier_, "biome_interaction_sensitivity")
                                   
                                   
                                   ### Test spatial autocorrelation -----
                                   library(spdep)
                                   
                                   dt_mod_sf <- dt_sub %>% 
                                     st_as_sf(.,
                                              coords = c("x_moll", "y_moll"), 
                                              crs = "ESRI:54009") %>% 
                                     mutate(resids_full = residuals(fit_full), 
                                            resids_fixed = residuals(fit_fixed)) %>% 
                                     filter(!is.infinite(resids_full), !is.infinite(resids_fixed))
                                   
                                   coords <- st_coordinates(dt_mod_sf)
                                   # knn <- knearneigh(coords, k = 250)
                                   # nb_knn <- knn2nb(knn)
                                   
                                   nb_knn <- dnearneigh(coords, 0, 100000)
                                   
                                   lw <- nb2listw(nb_knn, style = "W", zero.policy = T)
                                   
                                   mi_fixed <- moran.test(dt_mod_sf$resids_fixed, lw)
                                   mi_fixed
                                   
                                   mi_full <- moran.test(dt_mod_sf$resids_full, lw)
                                   mi_full
                                   
                                   
                                   
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
                                       morans_i_fixed = as.numeric(mi_fixed$estimate[1]), 
                                       morans_i_full = as.numeric(mi_full$estimate[1]), 
                                       morans_i_p_val_fixed = as.numeric(mi_fixed$p.value), 
                                       morans_i_p_val_full = as.numeric(mi_full$p.value),
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


#####################################################################################

## bind results 
unique(responses)
dt_res <- rbindlist(best_mesh_res_list) %>% 
  mutate(clean_response = case_when(
    .default = response,
    response == "evi_coef" ~ "EVI Trend",
  ), 
  clean_term = case_when(
    .default = term,
    term == "nitrogen_depo_scaled" ~ "Nitrogen Deposition",
    term == "mat_coef_scaled" ~ "MAT Trend",
    term == "prec_coef_scaled" ~ "Precipitation Trend",
    term == "max_temp_coef_scaled" ~ "Max. Temperature Trend",
    term == "hmi_change_scaled" ~ "HMI Change",
    term == "fire_frequency_scaled" ~ "Fire frequency",
    term == "mean_mat_scaled" ~ "MAT", 
    term == "mean_prec_scaled" ~ "MAP"))

unique(dt_res$clean_term)

summary(dt_res)

fwrite(dt_res, "builds/model_outputs/sdmtmb_results_world_grid.csv")

### Plot 

dt_res <- fread("builds/model_outputs/sdmtmb_results_world_grid.csv")
glimpse(dt_res)

dt_cr <- dt_res %>% 
  filter(tier == "climatic_region") %>% 
  mutate(
    biome = case_when(
      str_detect(term, "biome_colCold") ~ "Cold",
      str_detect(term, "biome_colPolar") ~ "Polar",
      str_detect(term, "biome_colTemperate") ~ "Temperate",
      str_detect(term, "biome_colTropical") ~ "Tropical",
      TRUE ~ "Reference (Arid)"),  #reference biome
    cont_var = case_when(
      str_detect(term, "nitrogen_depo_scaled") ~ "N Deposition",
      str_detect(term, "prec_coef_scaled") ~ "Precipitation trend",
      str_detect(term, "hmi_change_scaled") ~ "HMI Change",
      str_detect(term, "mat_coef_scaled") ~ "Mat Trend",
      str_detect(term, "max_temp_coef_scaled") ~ "Max. Temperature Trend",
      str_detect(term, "fire_frequency_scaled") ~ "Fire Frequency",
      str_detect(term, "mean_mat_scaled") ~ "MAT",
      str_detect(term, "mean_prec_scaled") ~ "MAP",
      TRUE ~ "Intercept"), 
    estimate_ci = paste0(round(estimate, 2), " (", round(conf.low, 2), "; ", round(conf.high, 2), ")")
  ) %>% 
  dplyr::select(Predictor = cont_var, Biome = biome, Estimate = estimate_ci) %>% 
  arrange(Predictor)
dt_cr

fwrite(dt_cr, "builds/model_outputs/clean_res_biome_interaction_climatic_regions.csv")

dt_sb <- dt_res %>% 
  filter(tier == "super_biome") %>% 
  mutate(
    biome = case_when(
      str_detect(term, "biome_colnot_cold_tall") ~ "Tall Vegetation, Not Cold Limited",
      str_detect(term, "biome_colnot_cold_short") ~ "Short Vegetation, Not Cold Limited",
      str_detect(term, "biome_colcold_tall") ~ "Tall Vegetation, Cold Limited",
      TRUE ~ "Reference (Short Vegetation, Cold Limited)"),  #reference biome
    cont_var = case_when(
      str_detect(term, "nitrogen_depo_scaled") ~ "N Deposition",
      str_detect(term, "prec_coef_scaled") ~ "Precipitation trend",
      str_detect(term, "hmi_change_scaled") ~ "HMI Change",
      str_detect(term, "mat_coef_scaled") ~ "Mat Trend",
      str_detect(term, "max_temp_coef_scaled") ~ "Max. Temperature Trend",
      str_detect(term, "fire_frequency_scaled") ~ "Fire Frequency",
      str_detect(term, "mean_mat_scaled") ~ "MAT",
      str_detect(term, "mean_prec_scaled") ~ "MAP",
      TRUE ~ "Intercept"), 
    estimate_ci = paste0(round(estimate, 2), " (", round(conf.low, 2), "; ", round(conf.high, 2), ")")
  ) %>% 
  dplyr::select(Predictor = cont_var, Biome = biome, Estimate = estimate_ci) %>% 
  arrange(Predictor)
dt_sb
fwrite(dt_sb, "builds/model_outputs/clean_res_biome_interaction_super_biome.csv")

