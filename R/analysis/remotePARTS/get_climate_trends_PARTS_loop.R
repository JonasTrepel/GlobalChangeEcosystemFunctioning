library(remotePARTS)
library(tidyverse)
library(data.table)

#param = "pas"
#param = "grid"
#param = "usa"
param = "europe"

if(param == "grid"){
  dt <- fread("data/processedData/dataFragments/grid_sample_with_raw_timeseries.csv") %>% 
    as.data.frame() 
} else if(param == "pas"){
  dt <- fread("data/processedData/dataFragments/pa_and_controls_with_raw_timeseries.csv") %>% 
      as.data.frame()
} else if(param == "usa") {
  dt <- fread("data/processedData/dataFragments/grid_usa_with_raw_timeseries.csv") %>% 
    as.data.frame()
} else if(param == "europe") {
  dt <- fread("data/processedData/dataFragments/grid_europe_with_raw_timeseries.csv") %>% 
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
  pattern = c("mat_", "max_temp_", "map_", "n_depo_zhu_", "n_depo_usa_", "n_depo_europe_"),
  name = c("mat", "max_temp", "map", "n_depo_zhu", "n_depo_usa", "n_depo_europe"),
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
if(param %in% c("grid", "usa", "europe")){
  dt <- dt %>% 
    mutate(chunk = "chunk") 
} else if(param == "pas"){
  dt <- dt %>% 
    mutate(chunk = Continent)
} 


chunks <- unique(dt$chunk)

############### create cluster and run loop in parallel ####################

library(doSNOW)
library(foreach)
library(tictoc)

# Create and register a cluster
clust <- makeCluster(5)
registerDoSNOW(clust)

## progress bar 
iterations <- nrow(trend_configs)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

dt_trend <- data.table()

for(chu in unique(chunks)){
  
print(paste0("starting with ", chu))
  
dt_chunk <- dt %>% filter(chunk %in% c(chu))

dt_trend_chunk <- foreach(i = 1:nrow(trend_configs),
                   .packages = c('tidyverse', 'remotePARTS', 'data.table', 'sf'),
                   .options.snow = opts,
                   .inorder = TRUE,
                   .verbose = TRUE,
                   .combine = left_join) %dopar% {

  config <- trend_configs[i, ]
 
  dt_sub <- process_trend(config$pattern, config$name, dt_chunk)
  
  return(dt_sub)
  
  print(paste0(config$name, " done! ", Sys.time()))
  
}


dt_trend <- rbind(dt_trend, dt_trend_chunk)

print(paste0(chu, " done! ", Sys.time()))

}

stopCluster(clust)
print(paste0("done! ", Sys.time()))



ctk <- dt %>% dplyr::select(unique_id,
                            mean_burned_area, max_burned_area, 
                            map_era, mat_era, mat_era,
                            mean_evi, mean_greenup)

dt_res <- dt %>% 
  as.data.table() %>% 
  left_join(dt_trend) %>%
  dplyr::select(-all_of(grep("mat_", names(dt), value = T)),
                -all_of(grep("max_temp_", names(dt), value = T)),
                -all_of(grep("map_", names(dt), value = T))) %>% 
  left_join(ctk)


summary(dt_res)

if(param == "grid"){
  fwrite(dt_res, "data/processedData/dataFragments/grid_sample_with_climate_trends.csv")
} else if(param == "pas"){
  fwrite(dt_res, "data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv")
} else if(param == "usa"){
  fwrite(dt_res, "data/processedData/dataFragments/grid_usa_with_climate_trends.csv")
} else if(param == "europe"){
  fwrite(dt_res, "data/processedData/dataFragments/grid_europe_with_climate_trends.csv")
}
