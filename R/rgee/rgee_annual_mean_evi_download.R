##### Load all rasters 😵‍💫 ##### 

library(rgee)
library(data.table)
library(tidyverse)
library(googledrive)
library(terra)

ee_Initialize(project = "ee-jonastrepel", drive = TRUE)
drive_auth(email = "jonas.trepel@bio.au.dk")


years <- c(2001:2023)

for(year in years){
  
  print(paste0("Starting with: ", year))
  
  start_date <- paste0(year, "-01-01")
  end_date <- paste0(year, "-12-31")
  
  annual_img <- ee$
    ImageCollection('MODIS/061/MOD13A1')$
    map(function(img) {
      qa <- img$select("SummaryQA")
      img$updateMask(qa$eq(0))})$ ## select only high quality data 
    select('EVI')$
    filterDate(start_date, end_date)$
    mean()
  
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
                                   folder = "rgee_backup",
                                   description = "annual_envi",
                                   scale = 500, 
                                   timePrefix = FALSE, 
                                   maxPixels = 1e13
  )
  export_task$start()
  
  #monitor task
  
  # Monitor the task status
  
  drive_auth(email = "jonas.trepel@bio.au.dk")
  
  for (i in 1:1000) {
    
    drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_envi") %>% 
      dplyr::select(name)
    
    # Check if the folder is empty
    if (n_distinct(drive_files) == 0) {
      Sys.sleep(30)  
      print(paste0("Attempt ", i, ": Drive still empty"))
    } else {
      print("Files found:")
      print(drive_files)
      
      if(n_distinct(drive_files) < 8){
        Sys.sleep(150) #to make sure all tiles are there
        drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_envi") %>% 
          dplyr::select(name)
      }
      #check again
      if(n_distinct(drive_files) < 8){
        Sys.sleep(150) #to make sure all tiles are there
      }
      drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_envi") %>%  dplyr::select(name)
      print(drive_files)
      
      break  #
    }
  }
  
  drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_envi") %>%
    dplyr::select(name) %>% 
    unique()
  
  
  for(filename in unique(drive_files$name)){
    
    path_name = paste0("data/rawData/raw_time_series/mean_evi/mean_evi_tmp_tiles/", filename)
    
    drive_download(file = filename, path = path_name, overwrite = TRUE)
    
  }
  
  googledrive::drive_rm(unique(drive_files$name))
  #googledrive::drive_empty_trash()
  
  
  files <- list.files("data/rawData/raw_time_series/mean_evi/mean_evi_tmp_tiles/", full.names = T)
  
  raster_list <- lapply(files, rast)
  
  file_name_merge <- paste0("data/rawData/raw_time_series/mean_evi/envi_mean_500m_", year, ".tif")
  
  global_evi <- do.call(merge, c(raster_list, list(filename = file_name_merge, overwrite = TRUE)))
  
  plot(global_evi, main = paste0(year))
  
  file.remove(files)
  
  print(paste0(year, " EVI mean done"))
  
}


