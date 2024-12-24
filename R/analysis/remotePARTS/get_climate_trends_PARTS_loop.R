library(remotePARTS)
library(tidyverse)
library(data.table)


# Read the data
dt <- fread("data/processedData/dataFragments/grid_sample_with_raw_timeseries.csv") %>% 
  as.data.frame() 

# dt <- fread("data/processedData/dataFragments/pa_and_controls_with_raw_timeseries.csv") %>% 
#  as.data.frame() %>% 
#  rename(X = Longitude, 
#         Y = Latitude)

#dt <- dt %>% sample_n(1000)

# Define a helper function to process trends
process_trend <- function(cols_pattern, trend_name, dt) {
  
  cols <- grep(cols_pattern, names(dt), value = TRUE)
  
  dt_subset <- dt %>% dplyr::select(all_of(cols), X, Y, unique_id) %>% 
    filter(complete.cases(.)) %>% as.data.frame()
  
  Y <- as.matrix(dt_subset[, cols])
  coords <- as.matrix(dt_subset[, c("X", "Y")])
  
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
  pattern = c("mat_", "max_temp_", "map_"#, "evi_", "burned_area_", "doy_greenup_1_"
              ),
  name = c("mat", "max_temp", "map"#, "evi", "burned_area", "doy_greenup_1"
           ),
  stringsAsFactors = FALSE
)

# 
for (i in 1:nrow(trend_configs)) {
  
  config <- trend_configs[i, ]
 
  dt_sub <- process_trend(config$pattern, config$name, dt)
  
  if(i == 1){
    dt_trend <- dt_sub
  }else{
    dt_trend <- left_join(dt_trend, dt_sub)
  }
  
  print(paste0(config$name, " done! ", Sys.time()))
  
}

print(paste0("done! ", Sys.time()))



dt_res <- dt %>% 
  as.data.table() %>% 
  left_join(dt_trend) %>%
  dplyr::select(-all_of(grep("mat_", names(dt), value = T)),
                -all_of(grep("max_temp_", names(dt), value = T)),
        #       -all_of(grep("evi_", names(dt), value = T)),
        #       -all_of(grep("burned_area_", names(dt), value = T)),
        #       -all_of(grep("doy_greenup_1_", names(dt), value = T)),
                -all_of(grep("map_", names(dt), value = T))
) 

fwrite(dt_res, "data/processedData/dataFragments/grid_sample_with_climate_trends.csv")
#fwrite(dt_res, "data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv")

