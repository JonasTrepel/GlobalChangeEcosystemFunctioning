#install.packages("geojsonio")
library(rgee)
library(data.table)
library(tidyverse)
library(tictoc)
ee_Initialize(project = "ee-jonastrepel", drive = TRUE)
rgee::ee_install_upgrade()
start_year <- 1950
end_year <- 2023

start_date <- paste0(start_year, "-01-01")
end_date <- paste0(end_year, "-12-31")

img_col <- ee$
  ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
  select('temperature_2m')$
  filterDate(start_date, end_date)


years = ee$List$sequence(ee$Number(start_year), ee$Number(end_year))

annualComposite <- function(year, col) {
  # Filter the collection for the given year
  annual <- col$filter(ee$Filter$calendarRange(year, year, "year"))$
    mean()$subtract(273.15) # Adjust metric as needed (e.g., convert from Kelvin to Celsius)
  
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

new_names <- paste0("Temp_", start_year:end_year)

annual_bands <- annual_col$toBands()$rename(new_names)

vis_params <- list(
  min = -10, max = 35, palette = c("blue", "yellow", "red")
)
Map$addLayer(annual_bands$select('Temp_2023'), vis_params, "Temperature 2023")


mean_temp <- annual_col$mean()

Map$addLayer(mean_temp, vis_params, "Temperature 2023")

#define world extent 
world_ext <- ee$Geometry$Rectangle(
  coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
  proj = "EPSG:4326",
  geodesic = FALSE
)

#load template raster (alternatively - build it from scratch)
library(terra)
library(sf)
mean_temp_r <- rast("/Users/jonas/Downloads/MatMean19502023.tif")
plot(mean_temp_r)


cell_coords <- terra::xyFromCell(mean_temp_r, cell = 1:ncell(mean_temp_r)) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "cell_number") %>% 
  mutate(temp = values(mean_temp_r)) %>% 
  filter(!is.na(temp)) %>% 
  dplyr::select(x, y)

coord_v <- st_as_sf(cell_coords, coords = c("x", "y"), crs = 4326)
coords <- st_coordinates(coord_v)


table(rep(1:ceiling(nrow(coords) / 1000), each = 1000)[1:nrow(coords)])

points <- cbind(coord_v, coords) %>% 
  mutate(chunk = rep(1:ceiling(nrow(.) / 1000), each = 1000)[1:nrow(.)]) 

temp_extr <- data.table()


for(ch in unique(points$chunk)){

  points_sub <- points %>% filter(chunk %in% ch) %>% st_as_sf()
  
  tic()
  tmp_extr <- ee_extract(x = annual_bands, y = points_sub, scale = 11132)
  toc()
  
  temp_extr <- rbind(tmp_extr, temp_extr)
  
  print(paste0(ch, " done"))

}




tic()
temp_extr <- ee_extract(annual_bands, coord_small, scale = 11132)
toc()

?ee_extract


plot(coord_v)

#define world extent 
world_ext <- ee$Geometry$Rectangle(
  coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
  proj = "EPSG:4326",
  geodesic = FALSE
)

pixel_centroids <- annual_bands$sample(
  region = world_ext, 
  geometries = TRUE
)

#valid_points <- pixel_centroids$filterBounds(ee$Geometry$Rectangle(-179.5, -89.4, 179.5, 89.4))

#limited_points <- valid_points

# Define the sampling function to extract values from each image
# sample_at_points <- function(image) {
#   image$sampleRegions(
#     collection = pixel_centroids,   # the point collection to sample
#     scale = 11132, #adjust scale if working on a different resolution 
#     geometries = TRUE               # returns geometries along with values
#   )
# }

# Map the sampling function across the image collection
#sampled_data <- annual_bands$map(sample_at_points)

sampled_data <- annual_bands$sampleRegions(
  collection = pixel_centroids,   # the point collection to sample
  scale = 11132, #adjust scale if working on a different resolution 
  geometries = TRUE)   

# Flatten the result to make it easier to handle as a table
#extracted_data <- sampled_data$flatten()


task <- ee_table_to_drive(
  collection = sampled_data,
  description = "era5_temperature_annual",
  folder = "RGEE",
  timePrefix = FALSE,
  
)

task$start()

dt_raw <- fread("/Users/jonas/Downloads/era5_temperature_annual.csv") 
  
dt <- dt_raw %>% 
  mutate(
  coords_weird = gsub('\\{""geodesic"":false,""type"":""Point"",""coordinates""\\:\\[', '', .geo), 
  coords_weird = gsub('\\]}', '', coords_weird), 
  X = as.numeric(str_split_i(coords_weird, pattern = ",", i = 1)),
  Y = as.numeric(str_split_i(coords_weird, pattern = ",", i = 2))) %>% 
  dplyr::select(-coords_weird, -.geo) %>% unique()

dt %>% ggplot() +
  geom_tile(aes(x = X, y = Y, fill = Temp2016)) +
  theme_void()

