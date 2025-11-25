#install.packages("remotePARTS")
library(remotePARTS)
library(tidyverse)
library(data.table)


dt <- fread("data/processedData/dataFragments/grid_sample_with_raw_timeseries.csv") %>% 
  as.data.frame()

#reshape2::melt(dt, measure = c("evi_2001", "evi_2011", "evi_2023")) %>% 
#  ggplot(aes(x = X, y = Y, col = value )) + 
#  geom_point(size = .1) +
#  labs(col = "evi") +
#  facet_wrap(~ gsub("evi_", "", variable), ncol = 3) +
#  scale_color_viridis_c(option = "magma") +
#  labs(x = "Longitude", y = "Latitude")


#### MEAN ANNUAL TEMPERATURE ####
# identify columns of interest 
mat_cols <- grep("mat_", names(dt), value = TRUE)

dt_mat <- dt %>% dplyr::select(all_of(mat_cols), X, Y, gridID) %>% 
  filter(complete.cases(.)) %>% as.data.frame()

str(dt_mat)
# make them a matrix 
Y_mat <- as.matrix(dt_mat[, mat_cols])

# get coordinate matrix 
coords_mat <- as.matrix(dt_mat[, c("X", "Y")])

# get time trends (accounts for temporal autocorrelation)
ar_mat <- fitAR_map(Y = Y_mat, coords = coords_mat)

# get coefficients 
dt_mat$mat_coef <- coefficients(ar_mat)[, "t"] 

# get p values 
dt_mat$mat_p_value <- ar_mat$pvals[,2]

dt_mat_trend <- dt_mat %>% 
  dplyr::select(mat_coef, mat_p_value, gridID)

print(paste0("mat done! ", Sys.time()))

#### MAX ANNUAL TEMPERATURE ####
#define columns of interest
max_temp_cols <- grep("max_temp_", names(dt), value = TRUE)

dt_max_temp <- dt %>% dplyr::select(all_of(max_temp_cols), X, Y, gridID) %>% 
  filter(complete.cases(.)) %>% as.data.frame()

# make them a matrix 
Y_max_temp <- as.matrix(dt_max_temp[, max_temp_cols])

# get coordinate matrix 
coords_max_temp <- as.matrix(dt_max_temp[, c("X", "Y")])

# get time trends (accounts for temporal autocorrelation)
ar_max_temp <- fitAR_map(Y = Y_max_temp, coords = coords_max_temp)

# get coefficients 
dt_max_temp$max_temp_coef <- coefficients(ar_max_temp)[, "t"] 

# get p values 
dt_max_temp$max_temp_p_value <- ar_max_temp$pvals[,2]

dt_max_temp_trend <- dt_max_temp %>% 
  dplyr::select(max_temp_coef, max_temp_p_value, gridID)

print(paste0("max_temp done! ", Sys.time()))


#### ANNUAL RAINFALL ####
#define columns of interest
map_cols <- grep("map_", names(dt), value = TRUE)

dt_map <- dt %>% dplyr::select(all_of(map_cols), X, Y, gridID) %>% 
  filter(complete.cases(.)) %>% as.data.frame()

# make them a matrix 
Y_map <- as.matrix(dt_map[, map_cols])

# get coordinate matrix 
coords_map <- as.matrix(dt_map[, c("X", "Y")])

# get time trends (accounts for temporal autocorrelation)
ar_map <- fitAR_map(Y = Y_map, coords = coords_map)

# get coefficients 
dt_map$map_coef <- coefficients(ar_map)[, "t"] 

# get p values 
dt_map$map_p_value <- ar_map$pvals[,2]

dt_map_trend <- dt_map %>% 
  dplyr::select(map_coef, map_p_value, gridID)

print(paste0("map done! ", Sys.time()))


#### EVI ####
#define columns of interest
evi_cols <- grep("evi_", names(dt), value = TRUE)

dt_evi <- dt %>% dplyr::select(all_of(evi_cols), X, Y, gridID) %>% 
  filter(complete.cases(.)) %>% as.data.frame()

# make them a matrix 
Y_evi <- as.matrix(dt_evi[, evi_cols])

# get coordinate matrix 
coords_evi <- as.matrix(dt_evi[, c("X", "Y")])

# get time trends (accounts for temporal autocorrelation)
ar_evi <- fitAR_map(Y = Y_evi, coords = coords_evi)

# get coefficients 
dt_evi$evi_coef <- coefficients(ar_evi)[, "t"] 

# get p values 
dt_evi$evi_p_value <- ar_evi$pvals[,2]

dt_evi_trend <- dt_evi %>% 
  dplyr::select(evi_coef, evi_p_value, gridID)

print(paste0("evi done! ", Sys.time()))

#### BURNED AREA ####

#define columns of interest
burned_area_cols <- grep("burned_area_", names(dt), value = TRUE)

dt_burned_area <- dt %>% dplyr::select(all_of(burned_area_cols), X, Y, gridID) %>% 
  filter(complete.cases(.)) %>% as.data.frame()

# make them a matrix 
Y_burned_area <- as.matrix(dt_burned_area[, burned_area_cols])

# get coordinate matrix 
coords_burned_area <- as.matrix(dt_burned_area[, c("X", "Y")])

# get time trends (accounts for temporal autocorrelation)
ar_burned_area <- fitAR_map(Y = Y_burned_area, coords = coords_burned_area)

# get coefficients 
dt_burned_area$burned_area_coef <- coefficients(ar_burned_area)[, "t"] 

# get p values 
dt_burned_area$burned_area_p_value <- ar_burned_area$pvals[,2]

dt_burned_area_trend <- dt_burned_area %>% 
  dplyr::select(burned_area_coef, burned_area_p_value, gridID)

print(paste0("burned_area done! ", Sys.time()))


#### START OF GROWING SEASON ####
#define columns of interest
doy_greenup_1_cols <- grep("doy_greenup_1_", names(dt), value = TRUE)

dt_doy_greenup_1 <- dt %>% dplyr::select(all_of(doy_greenup_1_cols), X, Y, gridID) %>% 
  filter(complete.cases(.)) %>% as.data.frame()

# make them a matrix 
Y_doy_greenup_1 <- as.matrix(dt_doy_greenup_1[, doy_greenup_1_cols])

# get coordinate matrix 
coords_doy_greenup_1 <- as.matrix(dt_doy_greenup_1[, c("X", "Y")])

# get time trends (accounts for temporal autocorrelation)
ar_doy_greenup_1 <- fitAR_map(Y = Y_doy_greenup_1, coords = coords_doy_greenup_1)

# get coefficients 
dt_doy_greenup_1$doy_greenup_1_coef <- coefficients(ar_doy_greenup_1)[, "t"] 

# get p values 
dt_doy_greenup_1$doy_greenup_1_p_value <- ar_doy_greenup_1$pvals[,2]

dt_doy_greenup_1_trend <- dt_doy_greenup_1 %>% 
  dplyr::select(doy_greenup_1_coef, doy_greenup_1_p_value, gridID)

print(paste0("doy_greenup_1 done! ", Sys.time()))



dt_trend <- dt %>% 
  as.data.table() %>% 
  dplyr::select(-all_of(mat_cols),
                -all_of(max_temp_cols),
                -all_of(map_cols),
                -all_of(evi_cols),
                -all_of(burned_area_cols),
                -all_of(doy_greenup_1_cols)) %>% 
  left_join(dt_mat_trend) %>% 
  left_join(dt_max_temp_trend) %>% 
  left_join(dt_map_trend) %>% 
  left_join(dt_evi_trend) %>% 
  left_join(dt_burned_area_trend) %>% 
  left_join(dt_doy_greenup_1_trend)

fwrite(dt_trend, "data/processedData/dataFragments/grid_sample_with_trends_seq.csv")



#plot significant trends 

dt_trend %>% 
  filter(evi_p_value < 0.05) %>% 
  ggplot(aes(x = X, y = Y, col = evi_coef)) + 
  geom_point(size = .1) + 
  scale_color_gradient2(high = "red", low = "blue", 
                        mid = "grey90", midpoint = 0) + 
  guides(fill = "none") + 
  labs(y = "Latitude", x = "Longitude", col = expression(beta[1])) +
  theme_void()
