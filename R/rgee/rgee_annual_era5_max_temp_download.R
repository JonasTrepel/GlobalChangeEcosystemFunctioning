##### Load all rasters ##### 



library(rgee)
library(data.table)
library(tidyverse)
library(googledrive)

source("R/functions/monitor_gee_task.R")
#  
rgee_env_dir <- c("C:\\Users\\au713983\\.conda\\envs\\rgee_env")
reticulate::use_python(rgee_env_dir, required=T)
#ee_clean_user_credentials()
#ee$Authenticate(auth_mode='notebook')
ee$Initialize(project = "jonas-trepel")
drive_auth(email = "jonas.trepel@gmail.com")
ee$String('Hello from the Earth Engine servers!')$getInfo()
years <- c(1951:2024)


for(year in years){
  
  print(paste0("Starting with: ", year))
  
  start_date <- paste0(year, "-01-01")
  end_date <- paste0(year, "-12-31")
  
  annual_img <- ee$
    ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
    select('temperature_2m_max')$
    filterDate(start_date, end_date)$
    max()$subtract(273.15)
  
  Map$addLayer(annual_img)
  
  # ee_print(annual_img)
  
  world_ext <- ee$Geometry$Rectangle(
    coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
    #coords = c(20.7, 41.3, 21, 41.5), 
    proj = "EPSG:4326",
    geodesic = FALSE
  )
  
  #Map$addLayer(world_ext)
  
  export_task <- ee_image_to_drive(image = annual_img,
                                   region = world_ext,
                                   folder = "rgee_backup_max_temp",
                                   description = "annual_max_temp",
                                   scale = 11132, 
                                   timePrefix = FALSE, 
                                   maxPixels = 1e13
  )
  export_task$start()
  
  Sys.sleep(300)
  monitor_gee_task(pattern = "annual_max_temp", path = "rgee_backup_max_temp",
                   last_sleep_time = 600, mail = "jonas.trepel@gmail.com")
  
  Sys.sleep(600)  
  
  
  drive_files <- drive_ls(path = "rgee_backup_max_temp", pattern = "annual_max_temp") %>%
    dplyr::select(name) %>% 
    unique()
  
  
  for(filename in unique(drive_files$name)){
    
    path_name = paste0("data/raw_data/raw_tiles/", filename)
    
    drive_download(file = filename, path = path_name, overwrite = TRUE)
    
  }
  
  googledrive::drive_rm(unique(drive_files$name))
  #googledrive::drive_empty_trash()
  
  
  files <- list.files("data/raw_data/raw_tiles/", full.names = T, pattern = "annual_max_temp")
  
  raster_list <- lapply(files, rast)
  
  file_name_merge <- paste0("data/raw_data/time_series/era5_max_temp_", year, ".tif")
  
  
  if(length(raster_list) > 1){
    global_max_temp<- do.call(merge, c(raster_list,
                                   list(filename = file_name_merge, overwrite = TRUE)))
  } else {
    global_max_temp <- raster_list[[1]]
    writeRaster(global_max_temp, filename = file_name_merge, overwrite = TRUE)
  }
  
  plot(global_max_temp, main = paste0(year))
  
  file.remove(files)
  
  print(paste0(year, " max_temp done"))
  
}


