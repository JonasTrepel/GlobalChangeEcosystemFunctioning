#install.packages("remotePARTS")
library(remotePARTS)
library(tidyverse)
library(data.table)
library(terra)
library(sf)

dt_raw <- fread("data/rawData/raw_time_series/data_frames/era5_mat_df.csv")

dt <- dt_raw %>% sample_n(10000) %>% 
  filter(complete.cases(.)) %>% 
  as.data.frame()
str(dt)


reshape2::melt(dt, measure = c("mat_1950", "mat_1990", "mat_2023")) %>% 
  ggplot(aes(x = x, y = y, col = value )) + 
  geom_point(size = .1) +
  labs(col = "temp") +
  facet_wrap(~ gsub("mat_", "", variable), ncol = 3) +
  scale_color_viridis_c(option = "magma") +
  labs(x = "Longitude", y = "Latitude")


### identify columns of interest 
mat_cols <- grep("mat", names(dt), value = TRUE)

# make them a matrix 
Y <- as.matrix(dt[, mat_cols])

### get coordinate matrix 
coords <- as.matrix(dt[, c("x", "y")])

## get time trends (accounts for temporal autocorrelation)
ARfit <- fitAR_map(Y = Y, coords = coords)

# get coefficients 
dt$AR_coef <- coefficients(ARfit)[, "t"] 

# get p values 
dt$AR_p_value <- ARfit$pvals[,2]

#plot significant trends 

dt %>% 
  filter(AR_p_value < 0.05) %>% 
  ggplot(aes(x = x, y = y, col = AR_coef)) + 
  geom_point(size = .1) + 
  scale_color_gradient2(high = "red", low = "blue", 
                        mid = "grey90", midpoint = 0) + 
  guides(fill = "none") + 
  labs(y = "Latitude", x = "Longitude", col = expression(beta[1])) +
  theme_void()

## now dealing with spatial autocorrelation 

D <- distm_scaled(coords)

corfit <- fitCor(resids = residuals(ARfit), coords = coords, covar_FUN = "covar_exp", 
                 start = list(range = 0.1), fit.n = 3000)

(range.opt = corfit$spcor)

V.opt <- covar_exp(D, range.opt)


(GLS.int <- fitGLS(AR_coef ~ 1, data = dt, V = V.opt, nugget = NA, no.F = TRUE))

GLS.int$
