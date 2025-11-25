library(terra)
library(tidyverse)

years <- c(2001:2024)

#https://thredds.met.no/thredds/catalog/data/EMEP/2025_Reporting/catalog.html
(files <- data.frame(filepath = list.files("data/raw_data/emep_nc/", full.names = T)))

for(year in years){
  
  year_cr <- as.character(year)
  
  path <- files %>% filter(grepl(year_cr, filepath)) %>% dplyr::select(filepath) %>% pull()
  
  r <- rast(path)
  
  if(year %in% c(2023, 2024)){
    
    
  r <- subset(r, c(1:58))
    
  } 
  
  n_layers <- subset(r, c("DDEP_OXN_m2Grid", "WDEP_OXN", "DDEP_RDN_m2Grid", "WDEP_RDN"))
    

  total_n_dep <- sum(n_layers)
  
  plot(total_n_dep, main = paste0(year))
  
  writeRaster(total_n_dep, paste0("data/raw_data/time_series/europe_n_depo_", year, ".tif"), overwrite=TRUE)
  
}

#unit: mgN/m2
