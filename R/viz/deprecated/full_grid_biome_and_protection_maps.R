library(sf)
library(data.table)
library(tidyverse)
library(scico)

grid <- read_sf("data/spatialData/protectedAreas/gridWithCovs.gpkg") %>% 
  filter(!is.na(functional_biome)) %>% 
  mutate(
    protection_cat_broad = case_when(
      iucn_cat %in% c("Ia", "Ib", "II") ~ "Strict", 
      iucn_cat %in% c("III", "IV", "V", "VI", "unknown_or_NA") ~ "Mixed",
      iucn_cat == "unprotected" ~ "Unprotected"), 
    super_biome = case_when(
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("T", functional_biome) ~ "cold_tall", 
      (grepl("C", functional_biome) | grepl("B", functional_biome)) & grepl("S", functional_biome) ~ "cold_short", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("T", functional_biome) ~ "not_cold_tall", 
      !grepl("C", functional_biome) & !grepl("B", functional_biome) & grepl("S", functional_biome) ~ "not_cold_short")
  ) %>%
  rename(unique_id = gridID)  %>%
  unique() %>% 
  distinct(X, Y, .keep_all = TRUE) %>%
  filter(!abs(lon) > 179.9)

world <- rnaturalearth::ne_countries() %>% filter(!name_en == "Antarctica") %>%
  st_transform(crs = 'ESRI:54030') %>% 
  mutate(geometry = st_make_valid(geometry), 
         world = "world") %>% 
  group_by(world) %>% 
  summarize()

grid_prot <- grid %>% 
  filter(!abs(lon) > 179.9) %>% 
  st_transform(crs = 'ESRI:54030') %>% 
  group_by(protection_cat_broad) %>% 
  summarize()

mapview::mapview(grid_prot)
## protected area figure 

p_prot <- ggplot() +
  geom_sf(data = world, fill = "grey75", color = "grey50") +
  geom_sf(data = grid_prot,
          aes(color = protection_cat_broad, fill = protection_cat_broad), alpha = 1) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  theme_void() +
  labs(color = "Protection", fill = "Protection") +
  theme(axis.title = element_blank())
p_prot

ggsave(plot = p_prot, "builds/plots/full_grid_prot_maps_shapes.png", dpi = 600)

## protected area figure 

grid_biome <- grid %>% 
  filter(!abs(lon) > 179.9) %>% 
  st_transform(crs = 'ESRI:54030') %>% 
  group_by(super_biome) %>% 
  summarize()

mapview::mapview(grid_biome)

p_biome <- ggplot() +
  geom_sf(data = world, fill = "grey75", color = "grey50") +
  geom_sf(data = grid_biome %>% 
            mutate(biome_clean = case_when(
              super_biome == "cold_short" ~ "Cold Limited\nShort Vegetation",
              super_biome == "cold_tall" ~ "Cold Limited\nTall Vegetation",
              super_biome == "not_cold_short" ~ "Not Cold Limited\nShort Vegetation",
              super_biome == "not_cold_tall" ~ "Not Cold Limited\nTall Vegetation")),
          aes(color = biome_clean, fill = biome_clean), alpha = 1) +
  scale_color_scico_d(palette = "batlowK") +
  scale_fill_scico_d(palette = "batlowK") +
  theme_void() +
  labs(color = "Biome", fill = "Biome") +
  theme(axis.title = element_blank())
p_biome

ggsave(plot = p_biome, "builds/plots/full_grid_biome_maps_shapes.png", dpi = 600)

