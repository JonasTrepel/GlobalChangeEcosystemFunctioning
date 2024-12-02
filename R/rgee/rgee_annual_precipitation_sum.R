##### Load all rasters 😵‍💫 ##### 

library(rgee)
library(data.table)
library(tidyverse)
library(googledrive)

ee_Initialize(project = "ee-jonastrepel", drive = TRUE)
drive_auth(email = "jonas.trepel@bio.au.dk")


years <- c(1950:2023)


for(year in years){
  
  print(paste0("Starting with: ", year))
  
  start_date <- paste0(year, "-01-01")
  end_date <- paste0(year, "-12-31")
  
  annual_img <- ee$
    ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
    select('total_precipitation_sum')$
    filterDate(start_date, end_date)$
    sum()$multiply(1000)
  
  Map$addLayer(annual_img, vis_params)
  
  world_ext <- ee$Geometry$Rectangle(
    coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
    #coords = c(20.7, 41.3, 21, 41.5), 
    proj = "EPSG:4326",
    geodesic = FALSE
  )
  
  Map$addLayer(world_ext)
  
  export_task <- ee_image_to_drive(image = annual_img,
                                   region = world_ext,
                                   folder = "rgee_backup",
                                   description = "annual_prec",
                                   scale = 11132, 
                                   timePrefix = FALSE, 
                                   maxPixels = 1e13
  )
  export_task$start()
  
  #monitor task
  
  # Monitor the task status
  task_monitor <- function(task) {
    status <- task$status()$state
    print(paste("Initial status:", status))  # Debugging output
    
    while (status %in% c("READY", "RUNNING")) {
      cat("Task is still running...\n")
      Sys.sleep(30)  # Wait 30 seconds before checking again
      status <- task$status()$state
      print(paste("Updated status:", status))  # Debugging output
    }
    
    # Final status update
    if (status == "COMPLETED") {
      cat("Task completed successfully!\n")
    } else {
      cat("Task failed or cancelled.\n")
    }
  }
  
  task_monitor(export_task)
  
  drive_auth(email = "jonas.trepel@bio.au.dk")
  drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_prec") %>% dplyr::select(name)
  
  
  for(filename in unique(drive_files$name)){
    
    path_name = paste0("data/rawData/raw_time_series/era5/era5_map/map_tmp_images/", filename)
    
    drive_download(file = filename, path = path_name, overwrite = TRUE)
    
  }
  
  googledrive::drive_rm(unique(drive_files$name))
  #googledrive::drive_empty_trash()
  
  
  files <- list.files("data/rawData/raw_time_series/era5/era5_map/map_tmp_images/", full.names = T)
  
  r1 <- rast(files[1])
  # r2 <- rast(files[2])
  # r3 <- rast(files[3])
  # r4 <- rast(files[4])
  # r5 <- rast(files[5])
  # r6 <- rast(files[6])
  # r7 <- rast(files[7])
  # r8 <- rast(files[8])
  
  file_name_merge <- paste0("data/rawData/raw_time_series/era5/era5_map/era_5_mat_", year, ".tif")
  
  writeRaster(r1, # r2, r3, r4, r5, r6, r7, r8,
              filename = file_name_merge, 
              overwrite = TRUE)
  
  plot(r1)
  
  print(paste0(year, " done"))
  
}
