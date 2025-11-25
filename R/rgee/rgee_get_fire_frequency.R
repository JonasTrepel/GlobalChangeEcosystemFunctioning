rgee_env_dir <- c("C:\\Users\\au713983\\.conda\\envs\\rgee_env")
reticulate::use_python(rgee_env_dir, required=T)

library(rgee)
library(data.table)
library(tidyverse)
library(googledrive)
library(terra)
library(exactextractr)

############################# HOUSEKEEPING #############################

source("R/functions/monitor_gee_task.R")

#ee_clean_user_credentials()
#ee$Authenticate(auth_mode='notebook')
#when on GIS04
#ee$Initialize(project = "ee-jonastrepel")
#drive_auth(email = "jonas.trepel@bio.au.dk")
#mail <- "jonas.trepel@bio.au.dk"
#when GIS07
ee$Initialize(project = "jonas-trepel")
drive_auth(email = "jonas.trepel@gmail.com")
mail <- "jonas.trepel@gmail.com"

ee$String('Hello from the Earth Engine servers!')$getInfo()


### Get sub-saharan Africa extent
aoi <- ee$Geometry$Rectangle(
  coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
  #coords = c(20.7, 41.3, 21, 41.5), 
  proj = "EPSG:4326",
  geodesic = FALSE
)


year_list <- ee$List$sequence(2001, 2024)
n_years <- 2024-2001

modis_burndate <- ee$ImageCollection("MODIS/061/MCD64A1")$
  select("BurnDate")

modis_qa <- ee$ImageCollection("MODIS/061/MCD64A1")$
  select("QA")

# Bitwise extraction function (bit 0: valid burn)
bitwise_extract <- function(input, from_bit, to_bit) {
  mask_size <- ee$Number(1)$add(to_bit)$subtract(from_bit)
  mask <- ee$Number(1)$leftShift(mask_size)$subtract(1)
  input$rightShift(from_bit)$bitwiseAnd(mask)
}

# Function to get annual binary burn map (1 = burned, 0 = unburned)
get_annual_binary <- function(year) {
  year <- ee$Number(year)

  # BurnDate and QA filtered by year
  start_date <- ee$Date$fromYMD(year, 1, 1)
  end_date <- ee$Date$fromYMD(year, 12, 31)
  
  burn_img <- modis_burndate$
    filterDate(start_date, end_date)$
    select("BurnDate")$
    max()
  
  qa_img <- modis_qa$
    filterDate(start_date, end_date)$
    select("QA")$
    max()
  
  # Mask: burn pixels with good quality
  mask <- bitwise_extract(qa_img, 0, 0)$eq(1)
  
  # Burned = 1, Unburned = 0
  burned_bin <- burn_img$
    where(burn_img$neq(0), 1)$
    unmask(0)$
    updateMask(mask)$
    rename("Burned")$
    set("system:time_start", start_date)
  
  return(burned_bin)
}

# Build the annual binary image collection
burned_col <- ee$ImageCollection$fromImages(
  year_list$map(ee_utils_pyfunc(function(yr) {
    get_annual_binary(yr)
  }))
)

# Sum the collection to get fire frequency
fire_frequency <- burned_col$sum()$clip(aoi)$round()$toDouble()

# Visualization
vis_params <- list(
  min = 0,
  max = 23,
  palette = c("#ffffff", "#ffffb2", "#fd8d3c", "#e31a1c", "#b10026")  # white to dark red
)
Map$centerObject(aoi, 6)
Map$addLayer(fire_frequency, vis_params, "Fire Frequency")


export_fire_frequency <- ee_image_to_drive(
  image = fire_frequency,
  region = aoi,
  folder = "rgee_backup_fire_frequency",
  description = "fire_frequency",
  scale = 500,
  timePrefix = FALSE,
  maxPixels = 1e13
)
export_fire_frequency$start()

Sys.sleep(300)
monitor_gee_task(pattern = "fire_frequency", path = "rgee_backup_fire_frequency",
                 mail = "jonas.trepel@gmail.com", last_sleep_time = 600)

Sys.sleep(600)
drive_files_fire_frequency <- drive_ls(path = "rgee_backup_fire_frequency", pattern = "fire_frequency") %>%
  dplyr::select(name) %>% 
  unique()


for(filename in unique(drive_files_fire_frequency$name)){
  
  path_name = paste0("data/raw_data/raw_tiles/", filename)
  
  drive_download(file = filename, path = path_name, overwrite = TRUE)
  
}

googledrive::drive_rm(unique(drive_files_fire_frequency$name))
#googledrive::drive_empty_trash()


files <- list.files("data/raw_data/raw_tiles/", full.names = T, pattern = "fire_frequency")

raster_list <- lapply(files, rast)

file_name_merge <- paste0("data/spatial_data/covariates/n_fires_500m_2001_2024.tif")


if(length(raster_list) > 1){
  global_fire_frequency <- do.call(merge, c(raster_list,
                                     list(filename = file_name_merge, overwrite = TRUE)))
} else {
  global_fire_frequency <- raster_list[[1]]
  writeRaster(global_fire_frequency, filename = file_name_merge, overwrite = TRUE)
}

plot(global_fire_frequency)

file.remove(files)

print(paste0("Fire frequency done ", Sys.time()))
