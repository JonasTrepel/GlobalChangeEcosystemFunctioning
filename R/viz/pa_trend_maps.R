### make nice maps ####

library(tidyverse)
library(data.table)
library(RColorBrewer)
library(scico)
library(rnaturalearth)
library(sf)

evi_trend <- fread("data/processedData/dataFragments/pa_evi_trends.csv")
burned_area_trend <- fread("data/processedData/dataFragments/pa_burned_area_trends.csv")
greenup_trend <- fread("data/processedData/dataFragments/pa_greenup_trends.csv")



world <- rnaturalearth::ne_countries() %>% filter(!name_en == "Antarctica") %>% st_transform(crs = 'ESRI:54030')

shapes <- st_read("data/spatialData/pas_and_controls.gpkg") %>%
  st_transform(crs = 'ESRI:54030') %>% 
  left_join(evi_trend) %>% 
  left_join(burned_area_trend) %>% 
  left_join(greenup_trend)

coords <- st_coordinates(st_centroid(shapes))
shapes$X <- coords[,1]
shapes$Y <- coords[,2]



### evi maps -------------
quantile(shapes$evi_coef, c(.025, .975), na.rm = T)


p_evi <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey95") +
  geom_sf(data = shapes %>%
            filter(!is.na(evi_coef)) %>% 
            mutate(evi_coef = ifelse(evi_coef > 39.1, 39.1, evi_coef),
                   evi_coef = ifelse(evi_coef < -22.5, -22.5, evi_coef)),
          aes(color = evi_coef, fill = evi_coef)) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  theme_minimal()
p_evi

ggsave(plot = p_evi, "builds/plots/evi_pa_shapes_map.png", dpi = 600)

p_evi_points <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey95") +
  geom_point(data = shapes %>%
            filter(!is.na(evi_coef)) %>% 
              filter(!is.na(evi_coef)) %>% 
              mutate(evi_coef = ifelse(evi_coef > 39.1, 39.1, evi_coef),
                     evi_coef = ifelse(evi_coef < -22.5, -22.5, evi_coef)),
          aes(x = X, y = Y, color = evi_coef, fill = evi_coef, size = area_km2), alpha = 0.7) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  theme_minimal() +
  theme(axis.title = element_blank())
p_evi_points

ggsave(plot = p_evi_points, "builds/plots/evi_pa_points_map.png", dpi = 600)

## burned area maps-------
quantile(shapes$burned_area_coef, c(.025, .975), na.rm = T)


p_burned_area <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey95") +
  geom_sf(data = shapes %>%
            filter(!is.na(burned_area_coef)) %>% 
            mutate(burned_area_coef = ifelse(burned_area_coef > 0.003, 0.003, burned_area_coef),
                   burned_area_coef = ifelse(burned_area_coef < -0.00257, -0.00257, burned_area_coef)),
          aes(color = burned_area_coef, fill = burned_area_coef)) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  theme_minimal()
p_burned_area

ggsave(plot = p_burned_area, "builds/plots/burned_area_pa_shapes_map.png", dpi = 600)

p_burned_area_points <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey95") +
  geom_point(data = shapes %>%
               filter(!is.na(burned_area_coef)) %>% 
               mutate(burned_area_coef = ifelse(burned_area_coef > 0.003, 0.003, burned_area_coef),
                      burned_area_coef = ifelse(burned_area_coef < -0.00257, -0.00257, burned_area_coef)),
             aes(x = X, y = Y, color = burned_area_coef, fill = burned_area_coef, size = area_km2), alpha = 0.7) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  theme_minimal() +
  theme(axis.title = element_blank())
p_burned_area_points

ggsave(plot = p_burned_area_points, "builds/plots/burned_area_pa_points_map.png", dpi = 600)

## greenup maps-------
quantile(shapes$greenup_coef, c(.025, .975), na.rm = T)


p_greenup <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey95") +
  geom_sf(data = shapes %>%
            filter(!is.na(greenup_coef) & !grepl("N", FunctionalBiome)) %>% 
            mutate(greenup_coef = ifelse(greenup_coef > 1.42, 1.42, greenup_coef),
                   greenup_coef = ifelse(greenup_coef < -1.88, -1.88, greenup_coef)),
          aes(color = greenup_coef, fill = greenup_coef)) +
  scale_color_scico(palette = "cork", midpoint = 0) +
  scale_fill_scico(palette = "cork", midpoint = 0) +
  theme_minimal()
p_greenup

ggsave(plot = p_greenup, "builds/plots/greenup_pa_shapes_map.png", dpi = 600)

p_greenup_points <- ggplot() +
  geom_sf(data = world, fill = "grey99", color = "grey95") +
  geom_point(data = shapes %>%
               filter(!is.na(greenup_coef) & !grepl("N", FunctionalBiome)) %>% 
               mutate(greenup_coef = ifelse(greenup_coef > 1.42, 1.42, greenup_coef),
                      greenup_coef = ifelse(greenup_coef < -1.88, -1.88, greenup_coef)),
             aes(x = X, y = Y, color = greenup_coef, fill = greenup_coef, size = area_km2), alpha = 0.7) +
  scale_color_scico(palette = "cork", midpoint = 0) +
  scale_fill_scico(palette = "cork", midpoint = 0) +
  theme_minimal() +
  theme(axis.title = element_blank())
p_greenup_points

ggsave(plot = p_greenup_points, "builds/plots/greenup_pa_points_map.png", dpi = 600)

