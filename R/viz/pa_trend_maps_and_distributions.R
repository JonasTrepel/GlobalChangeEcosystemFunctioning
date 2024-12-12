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


##### correlations ####

climate_trends <- fread("data/processedData/data_with_response_timeseries/pas_and_controls_with_climate_trends.csv") %>% 
  dplyr::select(unique_id, mat_coef, map_coef, max_temp_coef)

dt_corr <- shapes %>% 
  left_join(climate_trends) %>% 
  rename(pa_age = PaAge, 
         pa_area = area_km2, 
         nitrogen_depo = NitrogenDepo, 
         human_modification = HumanModification
         ) %>% 
  as.data.table() %>% 
  dplyr::select(
                evi_coef, burned_area_coef, greenup_coef,
                mat_coef, map_coef, max_temp_coef, 
                nitrogen_depo, human_modification) %>% 
  filter(complete.cases(.))

library(ggcorrplot)
corr <- round(cor(dt_corr), 1)
p_corr <- ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE)
p_corr
ggsave(plot = p_corr, "builds/plots/variable_correlations.png", dpi = 600, height = 10, width = 10)


#### trends in different biomes and protected vs unprotected 
library(ggridges)

#evi ridges
p_evi_prot <- ggplot() +
  geom_density_ridges_gradient(data = shapes %>%
                                 filter(!is.na(evi_coef)) %>% 
                                 mutate(evi_coef = ifelse(evi_coef > 39.1, 39.1, evi_coef),
                                        evi_coef = ifelse(evi_coef < -22.5, -22.5, evi_coef)),
                               aes(x = evi_coef, y = og_layer, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  theme_bw() +
  labs(y = "") + 
  theme(legend.position = "none")

p_evi_prot

p_evi_funbi <- ggplot() +
  geom_density_ridges_gradient(data = shapes %>%
                                 filter(!is.na(evi_coef)) %>% 
                                 mutate(evi_coef = ifelse(evi_coef > 39.1, 39.1, evi_coef),
                                        evi_coef = ifelse(evi_coef < -22.5, -22.5, evi_coef)),
                               aes(x = evi_coef, y = FunctionalBiome, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "bam", midpoint = 0) +
  scale_fill_scico(palette = "bam", midpoint = 0) +
  theme_bw() +
  labs(y = "") + 
  theme(legend.position = "none")

p_evi_funbi

p_evi_dens <- gridExtra::grid.arrange(p_evi_prot, p_evi_funbi, ncol = 1)


# burned area ridges
p_burned_area_prot <- ggplot() +
  geom_density_ridges_gradient(data = shapes %>%
                                 filter(!is.na(burned_area_coef)) %>% 
                                 mutate(burned_area_coef = ifelse(burned_area_coef > 0.003, 0.003, burned_area_coef),
                                        burned_area_coef = ifelse(burned_area_coef < -0.00257, -0.00257, burned_area_coef)),
                               aes(x = burned_area_coef, y = og_layer, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  theme_bw() +
  labs(y = "") + 
  theme(legend.position = "none")

p_burned_area_prot

p_burned_area_funbi <- ggplot() +
  geom_density_ridges_gradient(data = shapes %>%
                                 filter(!is.na(burned_area_coef)) %>% 
                                 mutate(burned_area_coef = ifelse(burned_area_coef > 0.003, 0.003, burned_area_coef),
                                        burned_area_coef = ifelse(burned_area_coef < -0.00257, -0.00257, burned_area_coef)),
                               aes(x = burned_area_coef, y = FunctionalBiome, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "vik", midpoint = 0) +
  scale_fill_scico(palette = "vik", midpoint = 0) +
  theme_bw() +
  labs(y = "") + 
  theme(legend.position = "none")

p_burned_area_funbi

p_burned_area_dens <- gridExtra::grid.arrange(p_burned_area_prot, p_burned_area_funbi, ncol = 1)

## greenup ridges 
p_greenup_prot <- ggplot() +
  geom_density_ridges_gradient(data = shapes %>%
               filter(!is.na(greenup_coef) & !grepl("N", FunctionalBiome)) %>% 
               mutate(greenup_coef = ifelse(greenup_coef > 1.42, 1.42, greenup_coef),
                      greenup_coef = ifelse(greenup_coef < -1.88, -1.88, greenup_coef)),
             aes(x = greenup_coef, y = og_layer, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "cork", midpoint = 0) +
  scale_fill_scico(palette = "cork", midpoint = 0) +
  theme_bw() +
  labs(y = "") + 
  theme(legend.position = "none")
  
p_greenup_prot

p_greenup_funbi <- ggplot() +
  geom_density_ridges_gradient(data = shapes %>%
                                 filter(!is.na(greenup_coef) & !grepl("N", FunctionalBiome)) %>% 
                                 mutate(greenup_coef = ifelse(greenup_coef > 1.42, 1.42, greenup_coef),
                                        greenup_coef = ifelse(greenup_coef < -1.88, -1.88, greenup_coef)),
                               aes(x = greenup_coef, y = FunctionalBiome, fill = ..x..), alpha = 0.7) +
  scale_color_scico(palette = "cork", midpoint = 0) +
  scale_fill_scico(palette = "cork", midpoint = 0) +
  theme_bw() +
  labs(y = "") + 
  theme(legend.position = "none")

p_greenup_funbi

p_greenup_dens <- gridExtra::grid.arrange(p_greenup_prot, p_greenup_funbi, ncol = 1)

p_ridges <- gridExtra::grid.arrange(p_evi_dens, p_burned_area_dens, p_greenup_dens, ncol = 3)
ggsave(plot = p_ridges, "builds/plots/pa_ridges.png", dpi = 600, height = 6, width = 12)
