#install.packages("geojsonio")
library(rgee)
library(data.table)
library(tidyverse)

ee_Initialize(project = "ee-jonastrepel", drive = TRUE)

start_year <- 2002
end_year <- 2023

start_date <- paste0(start_year, "-01-01")
end_date <- paste0(end_year, "-12-31")

img_col <- ee$
  ImageCollection('MODIS/061/MOD13A1')$
  map(function(img) {
    qa <- img$select("SummaryQA")
    img$updateMask(qa$eq(0))})$ ## select only high quality data 
  select('EVI')$
  filterDate(start_date, end_date)

vis_params <- list(
  min = 0, max = 9000, palette = c("brown", "yellow", "green")
)
Map$addLayer(img_col$median(), vis_params, "EVI")

years = ee$List$sequence(ee$Number(start_year), ee$Number(end_year))

annualComposite <- function(year, col) {
  # Filter the collection for the given year
  annual <- col$filter(ee$Filter$calendarRange(year, year, "year"))$
    median() # Adjust metric as needed 
  
  # Set year and time start properties
  annual <- annual$set("year", year)$
    set("system:time_start", ee$Date$fromYMD(year, 1, 1))
  
  return(annual)
}

annual_col <- ee$ImageCollection$fromImages(
  ee$List(years)$map(ee_utils_pyfunc(function(year) {
    annualComposite(year, img_col)
  }))
)


new_names <- paste0("EVI_", start_year:end_year)

annual_bands <- annual_col$toBands()$rename(new_names)
ee_print(annual_bands)

vis_params <- list(
  min = 0, max = 9000, palette = c("brown", "yellow", "green")
)
Map$addLayer(annual_bands$select('EVI_2023'), vis_params, "EVI annual meand")


native_scale <- annual_bands$projection()$nominalScale()$getInfo()
print(native_scale)



#define world extent 
world_ext <- ee$Geometry$Rectangle(
 # coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
  coords = c(20.7, 41.3, 21.5, 42.0),  # Smaller bounding box
  proj = "EPSG:4326",
  geodesic = FALSE
)

Map$addLayer(world_ext, vis_params, "EVI annual meand")

pixel_centroids <- annual_bands$sample(
  region = world_ext, 
  scale = 500,
  geometries = TRUE
)


size <- pixel_centroids$size()$getInfo()
print(size)

sampled_data <- annual_bands$sampleRegions(
  collection = pixel_centroids,   # the point collection to sample
  scale = 500, #adjust scale if working on a different resolution 
  geometries = TRUE)   


ee_print(sampled_data$limit(10))
# Flatten the result to make it easier to handle as a table
#extracted_data <- sampled_data$flatten()


task <- ee_table_to_drive(
  collection = sampled_data,
  description = "modis_evi_annual",
  folder = "RGEE",
  timePrefix = FALSE,
  
)

task$start()



dt_raw <- fread("/Users/jonas/Downloads/modis_evi_annual.csv") 

dt <- dt_raw %>% 
  mutate(
    coords_weird = gsub('\\{""geodesic"":false,""type"":""Point"",""coordinates""\\:\\[', '', .geo), 
    coords_weird = gsub('\\]}', '', coords_weird), 
    X = as.numeric(str_split_i(coords_weird, pattern = ",", i = 1)),
    Y = as.numeric(str_split_i(coords_weird, pattern = ",", i = 2))) %>% 
  dplyr::select(-coords_weird, -.geo) %>% unique() 


dt %>% ggplot() +
  geom_tile(aes(x = X, y = Y, fill = EVI_2023, color = EVI_2023)) +
  theme_void()
