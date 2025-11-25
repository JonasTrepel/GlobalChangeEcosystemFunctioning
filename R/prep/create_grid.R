
library(tidyverse)
library(sf)
library(rnaturalearth)
library(mapview)
library(terra)
library(data.table)

pas <- read_sf("data/spatial_data/protected_areas/pa_shapes.gpkg")

pas_t <- pas %>% st_transform(crs = "ESRI:54009") %>% 
  mutate(
    iucn_cat = case_when(
      iucn_cat_ord == 1 ~ "Ia",  
      iucn_cat_ord == 2 ~ "Ib",  
      iucn_cat_ord == 3 ~ "II",  
      iucn_cat_ord == 4 ~ "III",  
      iucn_cat_ord == 5 ~ "IV",  
      iucn_cat_ord == 6 ~ "V",  
      iucn_cat_ord == 7 ~ "VI",  
      iucn_cat_ord == 8 ~ "unknown_or_NA",  
    ))

############# skip this if first grid is already done #########

## Load world and make Grid
world <- rnaturalearth::ne_countries() %>%
  dplyr::select(continent, sovereignt, geometry) %>% 
  st_transform(crs = "ESRI:54009")

## For now 5km
world_grid <- st_make_grid(world, cellsize = c(5000, 5000)) %>% st_as_sf()

world_grid2 <- world_grid %>%
  st_join(world) %>%
  filter(!is.na(sovereignt)) %>% 
  mutate(unique_id = paste0("grid_", 1:nrow(.)))

# extract land_cover ----------------------------
# and remove urban and cultivated pixels
lc_num <-  c(0, 20, 30, 40, 50, 60, 70, 80, 90, 
            100, 111, 112, 113, 114, 115, 116,
            121, 122, 123, 124, 125, 126, 
            200)
land_cover <- c("unknown", 
               "shrubs", 
               "herbacous_veg", 
               "cultivated",
               "urban",
               "bare",
               "snow_ice",
               "water",
               "herbaceous_wetland",
               "moss_and_lichen",
               "closed_forest_evergreen_needle",
               "closed_forest_evergreen_broad",
               "closed_forest_deciduous_needle",
               "closed_forest_deciduous_broad",
               "closed_forest_mixed",
               "closed_forest_other",
               "open_forest_evergreen_needle",
               "open_forest_evergreen_broad",
               "open_forest_deciduous_needle",
               "open_forest_deciduous_broad",
               "open_forest_mixed",
               "open_forest_other",
               "ocean")

lc_leg <- data.frame(land_cover=land_cover, lc_num=lc_num)

lc <- rast("data/spatial_data/covariates/GlobalLandCovercopernicus2019.tif")
plot(lc)

grid_trans_lc <- world_grid2 %>%
  st_transform(crs(lc))

grid_trans_lc <- grid_trans_lc[!is.na(st_is_valid(grid_trans_lc)), ]

lc_extr <- exactextractr::exact_extract(lc,
                                       grid_trans_lc,
                                       append_cols = c("unique_id"),
                                       fun = "mode") #https://rdrr.io/cran/exactextractr/man/exact_extract.html, see under User-defined summary functions



lc_extr_fin <- as.data.table(lc_extr) %>% 
  rename(lc_num = `mode`) %>%
  left_join(lc_leg, by = "lc_num") %>% 
  dplyr::select(unique_id, land_cover) 

world_grid3 <- world_grid2 %>%
  left_join(lc_extr_fin) %>% 
  filter(!is.na(land_cover)) %>% 
  filter(!land_cover %in% c("unknown", "cultivated", "urban",
                            "ocean", "water", "snow_ice")) %>% 
  st_sf()

unique(world_grid3$land_cover)

#get grid coordinates
coords <- world_grid3 %>% st_centroid() %>% st_coordinates()

coords_lat_lon <- world_grid3 %>% st_transform(4326) %>% st_centroid() %>% st_coordinates()


world_grid3$x_moll <- coords[,1]
world_grid3$y_moll <- coords[,2]

world_grid3$lon <- coords_lat_lon[,1]
world_grid3$lat <- coords_lat_lon[,2]

glimpse(world_grid3)

st_write(obj = world_grid3, "data/spatial_data/grids/world_grid.gpkg", append = FALSE)
