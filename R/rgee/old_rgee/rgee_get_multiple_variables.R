#### rgee download / extract variables of interest ####

library(rgee)
library(data.table)
library(tidyverse)
library(googledrive)
library(terra)

ee_Initialize(project = "ee-jonastrepel", drive = TRUE)
drive_auth(email = "jonas.trepel@bio.au.dk")

### define parameters 

start_year <- 2001
end_year <- 2023

start_date <- paste0(start_year, "-01-01")
end_date <- paste0(end_year, "-12-31")

### Fire Frequency ------------------------

## get burn date band 
annual_img_col <- ee$ImageCollection("MODIS/061/MCD64A1")$
  filterDate(start_date, end_date)$
  select("BurnDate")

Map$addLayer(annual_img_col$mean())

# select quality band 
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
   # unmask(0)$ ## think carefully here - leaves non-burned areas NA is not activated (may be ok!)
    updateMask(mask)$
    set("system:time_start", ee$Date$fromYMD(year, 1, 1))
  
  return(fire_sub2)
}

years_ls <- ee$List$sequence(start_year, end_year)  # 

# Create an ImageCollection of annual burned area
annual_col <- ee$ImageCollection$fromImages(
  years_ls$map(ee_utils_pyfunc(function(anno) {
    yearly_fires(anno)
  }))
)

fire_frequency <- annual_col$sum()$divide((end_year-start_year))

# Visualize 

fire_viz_params = list(
  min = 0, max = 2, 
  palette = c('blue', 'yellow', 'orange', 'red'),
  opacity = 0.8 
)

Map$addLayer(fire_frequency, fire_viz_params, 'Fire Frequency');

## Export image: 
fire_frequency


### Burned area  ------------------------

## get burn date band 
annual_img_col <- ee$ImageCollection("MODIS/061/MCD64A1")$
  filterDate(start_date, end_date)$
  select("BurnDate")

# select quality band 
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

years_ls <- ee$List$sequence(start_year, end_year)  # 

# Create an ImageCollection of annual burned area
annual_col <- ee$ImageCollection$fromImages(
  years_ls$map(ee_utils_pyfunc(function(anno) {
    yearly_fires(anno)
  }))
)


burned_area <- annual_col$
    max()$
    setDefaultProjection(
      crs = "SR-ORG:6974",
      crsTransform = c(463.3127165279165, 0, -20015109.354, 
                       0, -463.3127165279167, 7783653.6376640005)
    )$
    reduceResolution(
      reducer = ee$Reducer$mean(), # Fractional cover 
      maxPixels = 255
    )$
    reproject(
      scale = 5000,
      crs = "SR-ORG:6974"
    )$
    rename("BurnedFraction")$
    set("year", year)$
    set("system:time_start", ee$Date$fromYMD(year, 1, 1))

fire_viz_params = list(
  min = 0, max = 1, 
  palette = c('blue', 'yellow', 'orange', 'red'),
  opacity = 0.8 
)

Map$addLayer(burned_area, fire_viz_params, 'Burned Area');

## Export image: 
burned_area


### Days since last fire  ------------------------

# Define reference date and current date (i.e., the date from which the days since the last fire should be calculated)
reference_date <- ee$Date('2001-01-01')
current_date <- ee$Date(end_date)

# Load the MODIS Fire Collection
modis_fire <- ee$ImageCollection("MODIS/061/MCD64A1")$
  filterDate(reference_date, current_date)

# select quality band 
qa_band <- ee$ImageCollection("MODIS/061/MCD64A1")$
  filterDate(reference_date, current_date)$
  select("QA")

bitwise_extract = function(input, fromBit, toBit) {
  maskSize = ee$Number(1)$add(toBit)$subtract(fromBit)
  mask = ee$Number(1)$leftShift(maskSize)$subtract(1)
  return(input$rightShift(fromBit)$bitwiseAnd(mask))
}


# Function to calculate the fire date for each image
add_dependents <- function(image) {
  # Get the image date
  img_date <- ee$Date(image$get("system:time_start"))
  img_year <- img_date$get("year")
  
  # Create a new date object representing January 1st of the given year
  year_start <- ee$Date$fromYMD(img_year, 1, 1)
  
  # Calculate the difference in days from the start of the year to the reference date
  day_of_year <- year_start$difference(reference_date, "day")$add(1)
  
  # Select the 'BurnDate' band and add the day of year to it
  fire_date <- image$select("BurnDate")$add(day_of_year)$int()
  
  # Filter QA band for the same year and create a mask
  qa_sub <- qa_band$filter(ee$Filter$calendarRange(img_year, img_year, "year"))$max()
  mask <- bitwise_extract(qa_sub, 0, 0)$eq(1)
  
  # Mask the fire date with the QA mask
  fire_date_masked <- fire_date$updateMask(mask)
  
  # Return the image with the new 'firedate' band
  image$addBands(fire_date_masked$rename("firedate"))
}

# Get the difference in days since the reference date
days_since_reference_date <- current_date$difference(reference_date, "day")$add(1)

# Apply the function to calculate the fire date for each image
fire_date <- modis_fire$map(add_dependents)$
  select("firedate")$
  max()

# Calculate days since the last fire
days_since_last_fire <- ee$Image(days_since_reference_date)$subtract(fire_date)$int()

# Visualization parameters
viz_params <- list(
  min = 100,
  max = 3000,
  palette = c("red", "yellow", "blue")
)

# Add layers to the map
Map$addLayer(days_since_last_fire, viz_params, "Days Since Last Fire")

## export image: 
days_since_last_fire

### EVI -------------------------

mean_evi <- ee$
  ImageCollection('MODIS/061/MOD13A1')$
  map(function(img) {
    qa <- img$select("SummaryQA")
    img$updateMask(qa$eq(0))})$ ## select only high quality data 
  select('EVI')$
  filterDate(start_date, end_date)$
  mean()

# Visualization parameters
evi_viz_params <- list(
  min = -2000,
  max = 9000,
  palette = c("brown", "yellow", "green", "forestgreen")
)

# Add layers to the map
Map$addLayer(mean_evi, evi_viz_params, "EVI")

## export image 
mean_evi

### NPP -------------------------

mean_npp <- ee$
  ImageCollection('MODIS/061/MOD17A3HGF')$
  select('Npp')$
  filterDate(start_date, end_date)$
  mean()

# Visualization parameters
npp_viz_params <- list(
  min = -5000,max = 30000,
  palette = c("brown", "yellow", "green", "forestgreen")
)

# Add layers to the map
Map$addLayer(mean_npp, npp_viz_params, "NPP")

## export image: 
mean_npp

### Elevation -----------------------

elevation <- ee$
  ImageCollection("COPERNICUS/DEM/GLO30")$
  select('DEM')$
  mean()

ele_viz_params <- list(
  min = -100, max = 8500,
  palette = c("brown", "yellow", "green", "white")
)

# Add layers to the map
Map$addLayer(elevation, ele_viz_params, "Elevation")

## export image: 
mean_npp

### Mean Temperature --------------------

mean_temperature <- ee$
  ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
  select('temperature_2m')$
  filterDate(start_date, end_date)$
  mean()$subtract(273.15)

temp_viz_params <- list(
  min = -10, max = 40,
  palette = c("blue", "yellow", "red")
)

# Add layers to the map
Map$addLayer(mean_temperature, temp_viz_params, "Mean Temp")

## export image: 
mean_temperature


### Max Temperature -------------------------
max_temperature <- ee$
  ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
  select('temperature_2m_max')$
  filterDate(start_date, end_date)$
  max()$subtract(273.15)

temp_viz_params <- list(
  min = -10, max = 40,
  palette = c("blue", "yellow", "red")
)

# Add layers to the map
Map$addLayer(max_temperature, temp_viz_params, "Max Temp")

## export image: 
max_temperature

### Min Temperature -------------------------
min_temperature <- ee$
  ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
  select('temperature_2m_min')$
  filterDate(start_date, end_date)$
  min()$subtract(273.15)

min_temp_viz_params <- list(
  min = -40, max = 20,
  palette = c("blue", "yellow", "red")
)

# Add layers to the map
Map$addLayer(min_temperature, min_temp_viz_params, "Min Temp")

## export image: 
min_temperature

### Precipitation -------------------------
# mean monthly precipitation #
mean_prec <- ee$
  ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
  select('total_precipitation_sum')$
  filterDate(start_date, end_date)$
  mean()$multiply(1000)

prec_viz_params <- list(
  min = 0, max = 500,
  palette = c("brown", "yellow", "lightblue", "blue")
)

# Add layers to the map
Map$addLayer(mean_prec, prec_viz_params, "Prec")

## export image: 
mean_prec

### Land cover ----------------------------

lc <- ee$
  ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global")$
  select("discrete_classification")$
  filterDate("2019-01-01", "2019-12-31")$median()
ee_print(lc)

lc_viz_params <- list(
  min = 0, max = 200,
  palette = c("brown", "yellow", "lightblue", "blue")
)

# Add layers to the map
Map$addLayer(lc, lc_viz_params, "Land cover")

## export image: 
lc
  

### Tolan Canopy Height ------------------

canopy_height <- ee$ImageCollection("projects/meta-forest-monitoring-okw37/assets/CanopyHeight")$
  select('cover_code')$
  mosaic()

ch_viz_params <- list(
  min = 0, max = 40,
  palette = c("brown", "yellow", "green", "forestgreen"), 
  opacity = 0.5
)

Map$addLayer(canopy_height, ch_viz_params, "Canopy Height")

### Export images of interest ####

fire_frequency # scale: 500m 
burned_area # scale: 5000m 
days_since_last_fire # scale: 500m 
mean_evi # 500m
mean_npp # 500m
elevation # 30m
mean_temperature #11132m
max_temperature #11132m
min_temperature #11132m
mean_prec #11132m
lc #100m
canopy_height #1m


### Export code 
#get world extent --> replace by any other region of interest. 
world_ext <- ee$Geometry$Rectangle(
  coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
  #coords = c(20.7, 41.3, 21, 41.5), 
  proj = "EPSG:4326",
  geodesic = FALSE
)

## adjust file path, names etc 

export_task <- ee_image_to_drive(image = , #INSERT IMAGE,
                                 region = world_ext,
                                 folder = "rgee_backup",
                                 description = 'var_export', #adjust as you want 
                                 scale = , ## adjust as wanted  
                                 timePrefix = FALSE, 
                                 maxPixels = 1e13
)
export_task$start()

#monitor task

# Monitor the task status

drive_auth(email = "jonas.trepel@bio.au.dk")

for (i in 1:1000) {
  
  drive_files <- drive_ls(path = "rgee_backup", pattern = "var_export") %>% 
    dplyr::select(name)
  
  # Check if the folder is empty
  if (n_distinct(drive_files) == 0) {
    Sys.sleep(30)  
    print(paste0("Attempt ", i, ": Drive still empty"))
  } else {
    print("Files found:")
    print(drive_files)
    
    if(n_distinct(drive_files) < 2){
      Sys.sleep(150) #to make sure all tiles are there
      drive_files <- drive_ls(path = "rgee_backup", pattern = "var_export") %>% 
        dplyr::select(name)
    }
    #check again
    if(n_distinct(drive_files) < 2){
      Sys.sleep(150) #to make sure all tiles are there
    }
    drive_files <- drive_ls(path = "rgee_backup", pattern = "var_export") %>%  dplyr::select(name)
    print(drive_files)
    
    break  #
  }
}

drive_files <- drive_ls(path = "rgee_backup", pattern = "var_export") %>%
  dplyr::select(name) %>% 
  unique()


for(filename in unique(drive_files$name)){
  
  path_name = paste0("", filename) ## adjust paths
  
  drive_download(file = filename, path = path_name, overwrite = TRUE)
  
}

googledrive::drive_rm(unique(drive_files$name))
#googledrive::drive_empty_trash()


files <- list.files("",
                    full.names = T, 
                    pattern = "var_export")


file_name_merge <- paste0("") #adjust 

if(n_distinct(files) > 1){
  
  raster_list <- lapply(files, rast)
  
  global_ba <- do.call(merge, c(raster_list, list(filename = file_name_merge, overwrite = TRUE)))
  
}else{
  global_ba <- rast(files)
  
  writeRaster(global_ba, filename = file_name_merge, overwrite = TRUE )
  
}


### Extract image information 

#1. load shapefile 
## shapefile <- read_sf(....)
#2. decide on image of interest 
## image <- image of interest 
# ?ee_extract
#3. extraction <- ee_extract(x = image, y = shapefile, scale = ...,  fun = ee$Reducer$mean(), )


