##### Load all rasters 😵‍💫 ##### 

### get phenology

library(rgee)
library(data.table)
library(tidyverse)
library(googledrive)
library(terra)
#ee_install_upgrade()
ee_Initialize(project = "ee-jonastrepel", drive = TRUE)
drive_auth(email = "jonas.trepel@bio.au.dk")


years <- c(2001:2023)

for(year in years){
  
  print(paste0("Starting with: ", year))
  
  start_date <- paste0(year, "-01-01")
  end_date <- paste0(year, "-12-31")
  
  # Define the annual MODIS image with quality mask
  annual_img <- ee$ImageCollection('MODIS/061/MCD12Q2')$
    map(function(img) {
      qa <- img$select("QA_Overall_1")
      img$updateMask(qa$lte(2)) # Select only high-quality data
    })$
    select('Greenup_1')$
    filterDate(start_date, end_date)$
    median() # Result is a single ee$Image
  
  
  ref_date <- ee$Date('1970-01-01') # Reference date for DOY calculation (all units are days since 1970)
  
  
  # Function to compute DOY
  doy_fun <- function(img) {
    # Calculate the difference in days from the reference date
    year_start_date <- ee$Date$fromYMD(year, 1, 1)
    day_diff <- year_start_date$difference(ref_date, 'day')
    
    # Subtract day_diff to calculate DOY
    doy <- img$subtract(day_diff)$byte()$
      set('year', year)$
      set('system:time_start', ee$Date$fromYMD(year, 1, 1))
    
    return(doy)
  }
  

  sos_doy <- doy_fun(annual_img)
 
  Map$addLayer(
    sos_doy,
    list(min = 0, max = 250, palette = c('yellow', 'green', 'blue')),
    'Day of Year'
  )
  

  world_ext <- ee$Geometry$Rectangle(
    coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
    #coords = c(20.7, 41.3, 21, 41.5), 
    proj = "EPSG:4326",
    geodesic = FALSE
  )
  

  export_task <- ee_image_to_drive(image = sos_doy,
                                   region = world_ext,
                                   folder = "rgee_backup",
                                   description = "annual_sos_doy_1",
                                   scale = 500, 
                                   timePrefix = FALSE, 
                                   maxPixels = 1e13
  )
  export_task$start()
  
  #monitor task
  
  # Monitor the task status
  
  drive_auth(email = "jonas.trepel@bio.au.dk")
  
  for (i in 1:1000) {
    
    drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_sos_doy_1") %>% 
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
        drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_sos_doy_1") %>% 
          dplyr::select(name)
      }
      #check again
      if(n_distinct(drive_files) < 8){
        Sys.sleep(150) #to make sure all tiles are there
      }
      drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_sos_doy_1") %>%  dplyr::select(name)
      print(drive_files)
      
      break  #
    }
  }
  
  drive_files <- drive_ls(path = "rgee_backup", pattern = "annual_sos_doy_1") %>%
    dplyr::select(name) %>% 
    unique()
  
  
  for(filename in unique(drive_files$name)){
    
    path_name = paste0("data/rawData/raw_time_series/phenology/phenology_tmp_tiles/", filename)
    
    drive_download(file = filename, path = path_name, overwrite = TRUE)
    
  }
  
  googledrive::drive_rm(unique(drive_files$name))
  #googledrive::drive_empty_trash()
  
  
  files <- list.files("data/rawData/raw_time_series/phenology/phenology_tmp_tiles/",
                      full.names = T, 
                      pattern = "annual_sos_doy_1")
  
  raster_list <- lapply(files, rast)
  
  file_name_merge <- paste0("data/rawData/raw_time_series/phenology/doy_greenup_1_modis_500m_", year, ".tif")
  
  global_doy <- do.call(merge, c(raster_list, list(filename = file_name_merge, overwrite = TRUE)))
  
  plot(global_doy, main = paste0(year))
  
  file.remove(files)
  
  print(paste0(year, "SOS DOY done"))
  
}


