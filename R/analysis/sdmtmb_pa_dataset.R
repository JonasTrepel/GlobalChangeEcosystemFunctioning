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

dt <- fread("data/processed_data/analysis_ready_pa_dataset.csv")

setDT(dt)


# get dataframe with comlete and clean data fro modeling 

acceptable_numbers = seq(1, 10000000, 5)
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
                fire_frequency, evi_coef)
corr <- round(cor(dt_corr), 2)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE, colors = c("#6D9EC1", "white", "#E46726"))


tiers <- c("prot_importance", 
             "prot_global_change") 

responses <- c("evi_coef")

mesh_guide_raw <- CJ(tier = tiers, 
                     response = responses) %>% 
  mutate(formula =  case_when(
    tier %in% unique("prot_importance") ~ "0 + protection_cat_broad + (1 | pair_id)", 
    tier %in% unique("prot_global_change") ~ "nitrogen_depo_scaled:protection_cat_broad +
    prec_coef_scaled:protection_cat_broad +
    hmi_change_scaled:protection_cat_broad +
    max_temp_coef_scaled:protection_cat_broad +
    mat_coef_scaled:protection_cat_broad +
    fire_frequency_scaled:protection_cat_broad +
    (1 | pair_id)"))



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
    cutoff <- mesh_guide[j, ]$cutoff
    max_inner <- mesh_guide[j, ]$max_inner_edge
    mesh_id <- mesh_guide[j, ]$mesh_id
    form <- mesh_guide[j, ]$formula
    
    
    message("Starting tier: ", tier, " at ", Sys.time())
    
    #acceptable_numbers <- seq(1, 10000000, by = divisor)
    
    dt_sub <- dt_mod %>%
      # filter(eval(parse(text = filter_call))) %>%
      # mutate(row_nr = 1:n()) %>%
      # filter(row_nr %in% acceptable_numbers) %>%
      mutate(
        nitrogen_depo_scaled = as.numeric(scale(nitrogen_depo)),
        fire_frequency_scaled = as.numeric(scale(fire_frequency)),
        mat_coef_scaled = as.numeric(scale(mat_coef)),
        prec_coef_scaled = as.numeric(scale(prec_coef)),
        hmi_change_scaled = as.numeric(scale(hmi_change)),
        max_temp_coef_scaled = as.numeric(scale(max_temp_coef)),
        pair_id = as.factor(pair_id)
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
      paste0(resp, " ~ ", form))
    
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
    
    cv_model_id <- paste0("cv_", tier_, "_", mesh_id, "_pa_dataset")
    
    result_row <- data.frame(
      tier = tier,
      formula = form,
      cutoff = cutoff,
      max_inner_edge = max_inner,
      mesh_id  = mesh_id,
      n_vertices  = nrow(mesh$mesh$loc),
      all_converged = fit_cv$converged,
      p_d_hessian = sum(fit_cv$pdHess),
      response = resp,
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

fwrite(dt_mesh_res_fin, "builds/model_outputs/cv_mesh_selection_sdmtmb_pa_dataset.csv")
n_distinct(dt_mesh_res_fin$model_id)


p_loglik <- dt_mesh_res_fin %>% 
 # filter(tier == "full_dataset_yes") %>% 
  ggplot() +
  geom_point(aes(x = cutoff, y = sum_loglik), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "B", direction = - 1, begin = 0.2, end = 0.8) +
  labs(y = "Log Likelihood (sum)", x = "Cutoff (km)", color = "Max\nInner\nEdge\n(km)") +
  facet_wrap(~tier, scales = "free") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = .5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "snow"), 
        strip.background = element_rect(fill = "linen", color = "linen"))
p_loglik
ggsave(plot = p_loglik, "builds/plots/supplement/sum_loglik_different_meshs_pa_dataset.png", dpi = 600, height = 4, width = 8)



##### run models --------------------------------------------

dt_mesh_res_fin <- fread("builds/model_outputs/cv_mesh_selection_sdmtmb_pa_dataset.csv")

dt_best_mesh <- dt_mesh_res_fin %>% 
  group_by(tier) %>% 
  slice_max(sum_loglik) %>% 
  ungroup()



### - Use best  Mesh ------------------

plan(multisession, workers = 2)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GiB
start_time <- Sys.time()

best_mesh_res_list <- future_map(1:nrow(dt_best_mesh),
                                 .progress = TRUE,
                                 .options = furrr_options(seed = TRUE),
                                 function(i) {
                                   
                                   tier = dt_best_mesh[i,]$tier
                                   resp = dt_best_mesh[i,]$response 
                                   form = dt_best_mesh[i,]$formula


                                   

                                   dt_sub <- dt_mod %>% 
                                     mutate(
                                       nitrogen_depo_scaled = as.numeric(scale(nitrogen_depo)),
                                       fire_frequency_scaled = as.numeric(scale(fire_frequency)),
                                       mat_coef_scaled = as.numeric(scale(mat_coef)),
                                       prec_coef_scaled = as.numeric(scale(prec_coef)),
                                       hmi_change_scaled = as.numeric(scale(hmi_change)), 
                                       max_temp_coef_scaled = as.numeric(scale(max_temp_coef)), 
                                       pair_id = as.factor(pair_id))
                                   
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
                                   
                                   re_formula <- as.formula(paste0(resp, "~ 1 + ( 1 | pair_id)"))
                                   
                                   form_fixed = gsub("+ ( 1 | pair_id)", "", form)
                                   fixed_formula <- as.formula(paste0(resp, " ~ ",form_fixed))
                                   
                                   full_formula <- as.formula(paste0(resp, " ~ ",form))
                                   
                                   #https://github.com/pbs-assess/sdmTMB/issues/466#issuecomment-3119589818
                                   
                                   #intercept only 
                                   fit_int <- sdmTMB(int_formula,
                                                     spatial = "off",
                                                     data = dt_sub,
                                                     mesh = mesh, 
                                                     reml = T)
                                   
                                   #intercept and random effect  
                                   fit_re <- sdmTMB(data = dt_sub,
                                                    mesh = mesh, 
                                                  spatial = "off",
                                                  formula = re_formula,
                                                  reml = T)

                                   #intercept and spatial  
                                   fit_sp <- sdmTMB(formula = int_formula,
                                                    data = dt_sub,
                                                    mesh = mesh, 
                                                    spatial = "on", 
                                                    reml = T)
                                   
                                   
                                   #full model 
                                   fit_full <- sdmTMB(spatial = "on", 
                                                      formula = full_formula, 
                                                      reml = T, 
                                                      data = dt_sub,
                                                      mesh = mesh)
                                   
                                   #fixed effects only 
                                   fit_fixed <- sdmTMB(spatial = "off", 
                                                       formula = fixed_formula, 
                                                       reml = T, 
                                                       data = dt_sub,
                                                       mesh = mesh)
                                   
                                 
                                   # total proportion deviance explained by our full model:
                                   (dev_explained_full <- 1 - deviance(fit_full) / deviance(fit_int))
                                   
                                   # proportion deviance explained by the mesh:
                                   (dev_explained_spatial <- 1 - deviance(fit_sp) / deviance(fit_int))
                                   
                                   # proportion deviance explained by the n deposition only:
                                   (dev_explained_fixed <- 1 - deviance(fit_fixed) / deviance(fit_int))
                                   
                                   # proportion deviance explained by random effect:
                                   (dev_explained_re <- 1 - deviance(fit_re) / deviance(fit_int))

                                   
                                   san <- sdmTMB::sanity(fit_full)
                                   
                                   model_id = paste0("subset_",tier, "_pa_dataset")
                                   
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
                                       dev_explained_re = dev_explained_re, 
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
    grepl("nitrogen_depo_scaled", term) ~ "Nitrogen Deposition",
    grepl("mat_coef_scaled", term) ~ "MAT Trend",
    grepl("prec_coef_scaled", term) ~ "Precipitation Trend",
    grepl("max_temp_coef_scaled", term) ~ "Max. Temperature Trend",
    grepl("hmi_change_scaled", term) ~ "HMI Change",
    grepl("fire_frequency_scaled", term) ~ "Fire frequency", 
    grepl("protection_cat_broadunprotected", term) ~ "Control", 
    grepl("protection_cat_broadstrict", term) ~ "Protected",
    grepl("(Intercept)", term) ~ "Intercept"),
  protection = case_when(
    grepl("strict", term) ~ "Protected", 
    grepl("unprotected", term) ~ "Control"
  ))
unique(dt_res$clean_term)
summary(dt_res)
fwrite(dt_res, "builds/model_outputs/sdmtmb_results_pa_dataset.csv")


dt_res %>% 
  filter(tier == "prot_importance") %>% 
  ggplot() +
  geom_pointrange(aes(x = estimate, y = clean_term,
                      xmin =  conf.low, xmax= conf.high, 
                      color = protection), 
                  linewidth = 1) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Estimate",
    y = NULL
  )


dt_res %>% 
  filter(tier != "prot_importance" & 
           !clean_term == "Intercept") %>% 
  ggplot(aes(x = estimate, y = clean_term, xmin = conf.low,  xmax = conf.high, color = protection
  )) +
  geom_pointrange(position = position_dodge(width = 0.5), 
                  linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Estimate",
    y = NULL
  )
