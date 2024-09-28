library(data.table)
library(tidyverse)
library(ggridges)
library(scico)
library(ggcorrplot)

pasCovsDTRaw <- fread("data/processedData/cleanData/pasWithCovs.csv") %>% 
  rename(PaAreaKm2 = GIS_AREA, 
         PaYear = STATUS_YR) %>% 
  mutate(PaYear = ifelse(PaYear == 0, NA, PaYear), 
         PaAge = 2023-PaYear) %>% 
  filter(complete.cases(HumanModification, NitrogenDepo,
                        SlopeMaxTemp, SlopeMeanTemp, SlopeMinTemp, 
                        SlopePrec, PaAge, PaAreaKm2))

table(pasCovsDTRaw$PaYear)

biomeN <- pasCovsDTRaw %>% 
  group_by(Biome) %>% 
  summarize(nPerBiome = n()) 

pasCovsDT <- pasCovsDTRaw %>% 
  left_join(biomeN) %>% 
  mutate(BiomeN = paste0(Biome, " (n = ", nPerBiome, ")"))  %>% 
  filter(nPerBiome > 50)
  

### Response variable distribution -------------------------------------
## Check out response variables 

quantile(pasCovsDT$EviTrend, na.rm = T, c(.99, .01))
quantile(pasCovsDT$BurnedAreaTrend, na.rm = T, c(.99, .01))
quantile(pasCovsDT$SOSTrend, na.rm = T, c(.99, .01))
quantile(pasCovsDT$NppTrend, na.rm = T, c(.99, .01))

pasCovsLong <- pasCovsDT %>% 
  filter(abs(EviTrend) < 30 & abs(BurnedAreaTrend) < 0.6 & abs(SOSTrend) < 5.5 & abs(NppTrend) < 150) %>% 
  pivot_longer(cols = c("BurnedAreaTrend", "EviTrend", "SOSTrend", "NppTrend"), 
               names_to = "Trend", values_to = "TrendValues") %>% 
  filter(nPerBiome > 50)

quantile(pasCovsDT$EviTrend, na.rm = T, c(.99, .01))

ggplot() +
  geom_density_ridges(data = pasCovsLong[,], aes(y = BiomeN, x = TrendValues, fill = Biome)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Trend, scales = "free_x", ncol = 2) +
  scale_fill_scico_d(palette = "bamako") +
  theme_classic() +
  theme(legend.position = "none")

names(pasCovsDT)  

### Correlation of Predictors (and the responses)
dtCorr <- pasCovsDT %>% 
  dplyr::select(EviTrend, SOSTrend, BurnedAreaTrend, 
                NppTrend, EviSDTrend, 
                SlopeMeanTemp, SlopeMaxTemp, SlopeMinTemp, SlopePrec, 
                NitrogenDepo, HumanModification, 
                BodyMassLoss, SpeciesLoss, PaAreaKm2, PaYear) %>% 
  filter(complete.cases(.)) %>% 
  filter(!is.infinite(BodyMassLoss))
corr <- round(cor(dtCorr), 1)


ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE)

quantile(pasCovsDT$SlopePrec, na.rm = T, c(.99, .01))
quantile(pasCovsDT$PaAreaKm2, na.rm = T, c(.98, .01))


pasCovsLongPred <- pasCovsDT %>% 
  filter(abs(EviTrend) < 30 & abs(BurnedAreaTrend) < 0.6 &
           abs(SOSTrend) < 5.5 & abs(NppTrend) < 150 & abs(SlopePrec) < 6 & 
           abs(PaAreaKm2)< 10000) %>% 
  pivot_longer(cols = c(SlopeMeanTemp, SlopeMaxTemp, SlopeMinTemp, SlopePrec, 
                        NitrogenDepo, HumanModification, PaAge, PaAreaKm2), 
               names_to = "Trend", values_to = "TrendValues") %>% 
  filter(nPerBiome > 50)

ggplot() +
  geom_density_ridges(data = pasCovsLongPred[,], aes(y = BiomeN, x = TrendValues, fill = Biome)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Trend, scales = "free_x", ncol = 4) +
  scale_fill_scico_d(palette = "bamako") +
  theme_classic() +
  theme(legend.position = "none")

names(pasCovsDT)  
