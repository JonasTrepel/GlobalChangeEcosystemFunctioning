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

ee_clean_user_credentials()
ee$Authenticate(auth_mode='notebook')
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

modis_lc <- ee$ImageCollection("MODIS/061/MCD12Q1")$
  select("LC_Type1")


# Function to get annual binary burn map (1 = burned, 0 = unburned)
get_annual_binary <- function(year) {
  year <- ee$Number(year)
  
  # LC_Type1 and QA filtered by year
  start_date <- ee$Date$fromYMD(year, 1, 1)
  end_date <- ee$Date$fromYMD(year, 12, 31)
  
  lc_img <- modis_lc$
    filterDate(start_date, end_date)$
    select("LC_Type1")$
    max()
  

  # Burned = 1, Unburned = 0
  lc_bin <- lc_img$
    remap(
      from = c(12, 13, 14, 15, 17),
      to   = rep(1, 5),
      defaultValue = 0
    )$
    #updateMask(mask)$
    rename("Exclude")$
    set("system:time_start", start_date)
  
  return(lc_bin)
}

# Build the annual binary image collection
lc_col <- ee$ImageCollection$fromImages(
  year_list$map(ee_utils_pyfunc(function(yr) {
    get_annual_binary(yr)
  }))
)


# Sum the collection to get fire frequency
lc_sum <- lc_col$sum()$clip(aoi)$round()$toDouble()

# Visualization
vis_params <- list(
  min = 0,
  max = 23,
  palette = c("#ffffff", "#ffffb2", "#fd8d3c", "#e31a1c", "#b10026")  # white to dark red
)
Map$centerObject(aoi, 6)
Map$addLayer(lc_sum, vis_params, "Years in Which Landcover Is Unsuitable")

lc_binary <- lc_sum$gte(1)$toDouble()

vis_params <- list(
  min = 0,
  max = 1,
  palette = c("#ffffff", "#fd8d3c")  # white to dark red
)
Map$centerObject(aoi, 6)
Map$addLayer(lc_binary, vis_params, "Years in Which Landcover Is Unsuitable")

export_lc <- ee_image_to_drive(
  image = lc_binary,
  region = aoi,
  folder = "rgee_backup_lc",
  description = "lc",
  scale = 500,
  timePrefix = FALSE,
  maxPixels = 1e13
)
export_lc$start()

Sys.sleep(300)
monitor_gee_task(pattern = "lc", path = "rgee_backup_lc",
                 mail = "jonas.trepel@gmail.com", last_sleep_time = 600)

Sys.sleep(600)
drive_files_lc <- drive_ls(path = "rgee_backup_lc", pattern = "lc") %>%
  dplyr::select(name) %>% 
  unique()


for(filename in unique(drive_files_lc$name)){
  
  path_name = paste0("data/raw_data/raw_tiles/", filename)
  
  drive_download(file = filename, path = path_name, overwrite = TRUE)
  
}

googledrive::drive_rm(unique(drive_files_lc$name))
#googledrive::drive_empty_trash()


files <- list.files("data/raw_data/raw_tiles/", full.names = T, pattern = "lc")

raster_list <- lapply(files, rast)

file_name_merge <- paste0("data/spatial_data/covariates/modis_unsuitable_landcover_500m_2001_2024.tif")


if(length(raster_list) > 1){
  global_lc <- do.call(merge, c(raster_list,
                                            list(filename = file_name_merge, overwrite = TRUE)))
} else {
  global_lc <- raster_list[[1]]
  writeRaster(global_lc, filename = file_name_merge, overwrite = TRUE)
}

plot(global_lc)

file.remove(files)

print(paste0("LC done ", Sys.time()))
