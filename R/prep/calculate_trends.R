library(remotePARTS)
library(tidyverse)
library(data.table)
library(future)
library(furrr)
library(tictoc)

param = "pas"
#param = "pa_controls"
#param = "grid"
#param = "usa"
#param = "europe"

if(param == "grid"){
  dt <- fread("data/processed_data/grid_with_timeseries.csv") %>% 
    as.data.frame() 
} else if(param == "pas"){
  dt <- fread("data/processed_data/pas_with_timeseries.csv") %>% 
    as.data.frame()
} else if(param == "pa_controls"){
  dt <- fread("data/processed_data/pa_controls_with_timeseries.csv") %>% 
    as.data.frame()
} else if(param == "usa") {
  dt <- fread("data/processed_data/grid_usa_with_timeseries.csv") %>% 
    as.data.frame()
} else if(param == "europe") {
  dt <- fread("data/processed_data/grid_europe_with_timeseries.csv") %>% 
    as.data.frame()
}



#dt <- dt %>% sample_n(1000)

# Define a helper function to process trends
process_trend <- function(cols_pattern, trend_name, dt) {
  
  cols <- grep(cols_pattern, names(dt), value = TRUE)
  
  dt_subset <- dt %>% dplyr::select(all_of(cols), lon, lat, unique_id) %>% 
    filter(complete.cases(.)) %>% as.data.frame()
  
  Y <- as.matrix(dt_subset[, cols])
  coords <- as.matrix(dt_subset[, c("lon", "lat")])
  
  ar_results <- fitAR_map(Y = Y, coords = coords)
  
  dt_subset[[paste0(trend_name, "_coef")]] <- coefficients(ar_results)[, "t"] 
  dt_subset[[paste0(trend_name, "_p_value")]] <- ar_results$pvals[, 2]
  
  
  dt_subset <- dt_subset %>% 
    dplyr::select(paste0(trend_name, "_coef"),
                  paste0(trend_name, "_p_value"),
                  unique_id)
  
  return(dt_subset)
}

# List of trends
trend_configs <- data.frame(
  pattern = c("mat_", "max_temp_", "prec_", "n_depo_zhu_", "usa_n_depo_", "n_depo_europe_", "evi_"),
  name = c("mat", "max_temp", "prec", "n_depo_zhu", "n_depo_usa", "n_depo_europe", "evi"),
  stringsAsFactors = FALSE
)

if(param != "usa"){
  trend_configs <- trend_configs %>% 
    filter(!name == "n_depo_usa")
}
if(param != "europe"){
  trend_configs <- trend_configs %>% 
    filter(!name == "n_depo_europe")
}

#define chunks and remove usa values for the main grids
if(param %in% c("pas", "pa_controls", "usa", "europe")){
  dt <- dt %>% 
    mutate(chunk_id = 1) 
} else if(param %in% c("grid")){
  chunk_size <- 250000
  dt$chunk_id <- ceiling(seq_len(nrow(dt)) / chunk_size)
  table(dt$chunk_id)
} 



############### create cluster and run loop in parallel ####################
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
plan(multisession, workers = 7)
tic()


all_trends_list <- list()
dt_trend <- data.table()

tic()

for(chunk in unique(dt$chunk)){
  
  print(paste0("Starting with chunk ", chunk, " of ", max(dt$chunk_id)))
  
  dt_chunk <- dt[dt$chunk_id == chunk, ]
  
  
  dt_trend_chunk_list <- future_map(1:nrow(trend_configs),
                                    .progress = TRUE,
                                    .options = furrr_options(seed = TRUE),
                                    function(i) {
                                      
                                      #for(i in 1:nrow(trend_configs)){
                                      config <- trend_configs[i, ]
                                      
                                      dt_sub <- process_trend(config$pattern, config$name, dt_chunk)
                                      
                                      return(dt_sub)
                                      
                                      print(paste0(config$name, " done! ", Sys.time()))
                                      
                                    }
  )
  
  dt_trend_chunk <- dt_trend_chunk_list %>%
    reduce(~ left_join(.x, .y, by = "unique_id"))
  
  dt_trend <- rbind(dt_trend, dt_trend_chunk)
  all_trends_list[[as.character(chunk)]] <- dt_trend_chunk
  
  
  print(paste0(chunk, " done! ", Sys.time()))
  
  rm(dt_trend_chunk)
  gc()
  
  
}

plan(sequential)
Sys.time()
dt_trend_from_list <- rbindlist(all_trends_list)


ctk <- dt %>% dplyr::select(unique_id,
                            mean_burned_area,
                            mean_prec, mean_mat, mean_max_temp,
                            mean_evi)

dt_res <- dt %>% 
  as.data.table() %>% 
  left_join(dt_trend) %>%
  dplyr::select(-all_of(grep("mat_", names(dt), value = T)),
                -all_of(grep("max_temp_", names(dt), value = T)),
                -all_of(grep("prec_", names(dt), value = T)),
                -all_of(grep("burned_area_", names(dt), value = T)),
                -all_of(grep("n_depo_zhu_", names(dt), value = T))) %>% 
  left_join(ctk)


summary(dt_res)

if(param == "grid"){
  fwrite(dt_res, "data/processed_data/grid_with_trends.csv")
} else if(param == "pas"){
  fwrite(dt_res, "data/processed_data/pas_with_trends.csv")
} else if(param == "pa_controls"){
  fwrite(dt_res, "data/processed_data/pa_controls_with_trends.csv")
} else if(param == "usa"){
  fwrite(dt_res, "data/processed_data/grid_usa_with_trends.csv")
} else if(param == "europe"){
  fwrite(dt_res, "data/processed_data/grid_europe_with_trends.csv")
}
