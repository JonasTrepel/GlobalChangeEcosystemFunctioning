library(data.table)
library(tidyverse)
library(ggridges)
library(scico)
library(ggcorrplot)
library(gridExtra)

pasCovsDTRaw <- fread("data/processedData/cleanData/pasWithCovs.csv") %>% 
  rename(PaAreaKm2 = GIS_AREA, 
         PaYear = STATUS_YR) %>% 
  mutate(PaYear = ifelse(PaYear == 0, NA, PaYear), 
         PaAge = 2023-PaYear) %>% 
  filter(complete.cases(HumanModification, NitrogenDepo,
                        SlopeMaxTemp, SlopeMeanTemp, SlopeMinTemp, 
                        SlopePrec, PaAge, PaAreaKm2))

table(pasCovsDTRaw$PaYear)

ClimaticRegionN <- pasCovsDTRaw %>% 
  group_by(ClimaticRegion) %>% 
  summarize(nPerClimaticRegion = n()) 

pasCovsDT <- pasCovsDTRaw %>% 
  left_join(ClimaticRegionN) %>% 
  mutate(ClimaticRegionN = paste0(ClimaticRegion, " (n = ", nPerClimaticRegion, ")"))  %>% 
  filter(nPerClimaticRegion > 50)
  

### Response variable distribution -------------------------------------
## Check out response variables 

quantile(pasCovsDT$EviTrend, na.rm = T, c(.975, .025))
quantile(pasCovsDT$BurnedAreaTrend, na.rm = T, c(.975, .025))
quantile(pasCovsDT$SOSTrend, na.rm = T, c(.975, .025))
quantile(pasCovsDT$NppTrend, na.rm = T, c(.975, .025))

pasCovsLong <- pasCovsDT %>% 
  filter(EviTrend < 27.7 & EviTrend > -18.7 &
           BurnedAreaTrend < 0.32 & BurnedAreaTrend > -.32 & 
           SOSTrend < 3.12 & SOSTrend > -2.95 & 
           NppTrend < 91.9 & NppTrend > -59.1) %>% 
  pivot_longer(cols = c("BurnedAreaTrend", "EviTrend", "SOSTrend", "NppTrend"), 
               names_to = "TrendName", values_to = "TrendValue") %>% 
  filter(!(TrendName == "BurnedAreaTrend" & BurnedAreaMean == 0)) %>% 
  filter(nPerClimaticRegion > 50)

quantile(pasCovsDT$EviTrend, na.rm = T, c(.99, .01))

pResp <- ggplot() +
  geom_density_ridges(data = pasCovsLong[,], aes(y = ClimaticRegionN, x = TrendValue, fill = ClimaticRegion)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~TrendName, scales = "free_x", ncol = 2) +
  scale_fill_scico_d(palette = "bamako") +
  theme_classic() +
  theme(legend.position = "none")
pResp
names(pasCovsDT)  

### Correlation of Predictors (and the responses)
dtCorr <- pasCovsDT %>% 
  dplyr::select(EviTrend, SOSTrend, BurnedAreaTrend, 
                NppTrend,
                SlopeMeanTemp, SlopeMaxTemp, SlopeMinTemp, SlopePrec, 
                NitrogenDepo, HumanModification, PaAreaKm2) %>% 
  filter(complete.cases(.))
corr <- round(cor(dtCorr), 1)


pCorr <- ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE)

ggsave(plot = pCorr, "builds/plots/supplement/correlationOfVars.png", dpi = 600, height = 10)

quantile(pasCovsDT$SlopePrec, na.rm = T, c(.99, .01))
quantile(pasCovsDT$PaAreaKm2, na.rm = T, c(.98, .01))


pasCovsLongPred <- pasCovsDT %>% 
  filter() %>% 
  pivot_longer(cols = c(SlopeMeanTemp, SlopeMaxTemp, SlopeMinTemp, SlopePrec, 
                        NitrogenDepo, HumanModification, PaAreaKm2), 
               names_to = "TrendName", values_to = "TrendValue") %>% 
  filter(nPerClimaticRegion > 50)

### limit the predictor X axis to the 95 % interval 
predictorLimits <- pasCovsLongPred %>%
  group_by(TrendName) %>%
  summarize(
    lowerQuantile = quantile(TrendValue, 0.025, na.rm = T),
    upperQuantile = quantile(TrendValue, 0.975, na.rm = T),
  ) 

for(clean.var in unique(predictorLimits$TrendName)){
  
  upperLim <- predictorLimits[predictorLimits$TrendName == clean.var, ]$upperQuantile
  lowerLim <- predictorLimits[predictorLimits$TrendName == clean.var, ]$lowerQuantile
  
  pasCovsLongPredSub <- pasCovsLongPred[pasCovsLongPred$TrendName == clean.var,] %>% filter(TrendValue > lowerLim & TrendValue < upperLim)
  pasCovsLongPred <- pasCovsLongPred %>% 
    filter(!TrendName == clean.var) %>% 
    rbind(pasCovsLongPredSub)
  
  print(clean.var)
  
}

pDrivers <- ggplot() +
  geom_density_ridges(data = pasCovsLongPred[,], aes(y = ClimaticRegionN, x = TrendValue, fill = ClimaticRegion)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~TrendName, scales = "free_x", ncol = 4) +
  scale_fill_scico_d(palette = "bamako") +
  theme_classic() +
  theme(legend.position = "none")
pDrivers

pDist <- grid.arrange(pResp, pDrivers)
ggsave(plot = pDist, "builds/plots/supplement/distributionOfVars.png", dpi = 600, height = 10)
names(pasCovsDT)  
