library(terra)
library(tidyverse)

years <- c(2001:2023)

#https://thredds.met.no/thredds/catalog/data/EMEP/2024_Reporting/catalog.html
files <- data.frame(filepath = list.files("data/rawData/raw_time_series/europe_n_depo/emep_nc/", full.names = T))

for(year in years){
  
  year_cr <- as.character(year)
  
  path <- files %>% filter(grepl(year_cr, filepath)) %>% dplyr::select(filepath) %>% pull()
  
  r <- rast(path)
  
  if(year %in% c(2022, 2023)){
    
    
  r <- subset(r, c(1:59))
    
  } 
  
  n_layers <- subset(r, c("DDEP_OXN_m2Grid", "WDEP_OXN", "DDEP_RDN_m2Grid", "WDEP_RDN"))
    

  total_n_dep <- sum(n_layers)
  
  plot(total_n_dep, main = paste0(year))
  
  writeRaster(total_n_dep, paste0("data/rawData/raw_time_series/europe_n_depo/europe_n_depo_", year, ".tif"), overwrite=TRUE)
  
}

#unit: mgN/m2
