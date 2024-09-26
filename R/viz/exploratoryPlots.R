library(data.table)
library(tidyverse)
library(ggridges)


pasCovsDTRaw <- fread("data/processedData/cleanData/pasWithCovs.csv") 

biomeN <- pasCovsDTRaw %>% 
  group_by(Biome) %>% 
  summarize(nPerBiome = n())

pasCovsDT <- pasCovsDTRaw %>% 
  left_join(biomeN) %>% 
  mutate(BiomeN = paste0(Biome, " (n = ", nPerBiome, ")"))

quantile(pasCovsDT$EviTrend, na.rm = T, c(.99, .01))

ggplot() +
  geom_density_ridges(data = pasCovsDT[abs(pasCovsDT$EviTrend) < 30,], aes(y = BiomeN, x = EviTrend, fill = Biome)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(legend.position = "none")
  
