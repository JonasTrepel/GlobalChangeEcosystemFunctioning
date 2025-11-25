
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

years <- c(2001:2024)


for(year in years){
  
  print(paste0("Starting with: ", year))
  
  start_date <- paste0(year, "-01-01")
  end_date <- paste0(year, "-12-31")
  
  ## get burn date band 
  annual_img_col <- ee$ImageCollection("MODIS/061/MCD64A1")$
    filterDate(start_date, end_date)$
    select("BurnDate")
  
  Map$addLayer(annual_img_col$mean())
  
  # Define the image collection and select QA band
  QA <- ee$ImageCollection("MODIS/061/MCD64A1")$
    filterDate(start_date, end_date)$
    select("QA")
  
  #print(QA$getInfo())
  
  # Function to modify the BurnDate band
  set_non_zero_to_one <- function(image) {
    band <- image$select("BurnDate")
    modifiedBand <- band$where(band$neq(0), 1)
    image$addBands(modifiedBand, NULL, TRUE)$rename("BurnDateModified")
  }
  
  # Apply the function to each image in the collection

  annual_01 <- annual_img_col$map(set_non_zero_to_one)
  Map$addLayer(annual_01$max())
  
 bitwise_extract = function(input, fromBit, toBit) {
    maskSize = ee$Number(1)$add(toBit)$subtract(fromBit)
    mask = ee$Number(1)$leftShift(maskSize)$subtract(1)
    return(input$rightShift(fromBit)$bitwiseAnd(mask))
  }
  
 # Function to count the number of fires in a pixel per year
yearly_fires <- function(year) {
    fire_sub <- annual_01$
      filter(ee$Filter$calendarRange(year, year, "year"))$
      max()
    
    qa_sub <-  QA$filter(ee$Filter$calendarRange(year, year, 'year'))$max();
    
    mask <- bitwise_extract(qa_sub, 0, 0)$eq(1)
    
    fire_sub2 <- fire_sub$
      unmask(0)$
      updateMask(mask)$
      set("system:time_start", ee$Date$fromYMD(year, 1, 1))
    
    return(fire_sub2)
  }
  
  years_ls <- ee$List$sequence(year, year)  # 
  
  # Create an ImageCollection of annual burned area
  annual_col <- ee$ImageCollection$fromImages(
    years_ls$map(ee_utils_pyfunc(function(anno) {
      yearly_fires(anno)
    }))
  )
  
  # Function to reduce resolution of images and calculate fraction burned area
  burned_area_fun <- function(year) {

    annual_burned_sub <- annual_col$
      filter(ee$Filter$calendarRange(year, year, "year"))$
      max()$
      setDefaultProjection(
        crs = "SR-ORG:6974",
        scale = 500
      )$
      reduceResolution(
        reducer = ee$Reducer$mean(), # Fractional cover 
        maxPixels = 255
      )$
      reproject(
        scale = 5000,
        crs = "SR-ORG:6974"
      )$
      rename("BurnedFraction")
    
    annual_burned_sub$
      set("year", year)$
      set("system:time_start", ee$Date$fromYMD(year, 1, 1))
  }
  
  # Create an ImageCollection from the function applied to each year
  burned_area <- ee$ImageCollection$fromImages(
    years_ls$map(ee_utils_pyfunc(function(anno) {
      burned_area_fun(anno)
    }))
  )
  
  # Compute mean burned fraction and display
  mean_burned <- burned_area$reduce(ee$Reducer$mean())
  
  #get world extent
  world_ext <- ee$Geometry$Rectangle(
    coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
    #coords = c(20.7, 41.3, 21, 41.5), 
    proj = "EPSG:4326",
    geodesic = FALSE
  )
  
  export_task <- ee_image_to_drive(image = mean_burned,
                                   region = world_ext,
                                   folder = "rgee_backup_burned_area",
                                   description = "annual_burned_area",
                                   scale = 5000, 
                                   timePrefix = FALSE, 
                                   maxPixels = 1e13
  )
  export_task$start()
  
  Sys.sleep(300)
  monitor_gee_task(pattern = "annual_burned_area", path = "rgee_backup_burned_area",
                   last_sleep_time = 300, mail = "jonas.trepel@gmail.com")
  
  Sys.sleep(60)  
  
  
  drive_files <- drive_ls(path = "rgee_backup_burned_area", pattern = "annual_burned_area") %>%
    dplyr::select(name) %>% 
    unique()
  
  
  for(filename in unique(drive_files$name)){
    
    path_name = paste0("data/raw_data/raw_tiles/", filename)
    
    drive_download(file = filename, path = path_name, overwrite = TRUE)
    
  }
  
  googledrive::drive_rm(unique(drive_files$name))
  #googledrive::drive_empty_trash()
  
  
  files <- list.files("data/raw_data/raw_tiles/", full.names = T, pattern = "annual_burned_area")
  
  raster_list <- lapply(files, rast)
  
  file_name_merge <- paste0("data/raw_data/time_series/burned_area_modis_5000m_", year, ".tif")
  
  
  if(length(raster_list) > 1){
    global_burned_area<- do.call(merge, c(raster_list,
                                  list(filename = file_name_merge, overwrite = TRUE)))
  } else {
    global_burned_area <- raster_list[[1]]
    writeRaster(global_burned_area, filename = file_name_merge, overwrite = TRUE)
  }
  
  
  plot(global_burned_area, main = paste0(year))
  
  file.remove(files)
  
  print(paste0(year, " Burned area done"))
  
}
