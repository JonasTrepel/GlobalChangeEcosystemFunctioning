Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


library(data.table)
library(tidyverse)
library(ggridges)
library(scico)
library(mgcv)
library(MuMIn)
library(broom)
library(tictoc)
library(gridExtra)



## load data 
pasCovsDTRaw <- fread("data/processedData/cleanData/saReservesWithCovs.csv") %>% 
  rename(PaAreaKm2 = area_ha) %>% 
   filter(complete.cases(HumanModification, NitrogenDepo,
                        SlopeMaxTemp, SlopeMeanTemp, SlopeMinTemp, 
                        SlopePrec, PaAreaKm2, EviTrend))

ClimaticRegionN <- pasCovsDTRaw %>% 
  group_by(ClimaticRegion) %>% 
  summarize(nPerClimaticRegion = n()) 

dtMod <- pasCovsDTRaw %>% 
  left_join(ClimaticRegionN) %>% 
  mutate(ClimaticRegionN = paste0(ClimaticRegion, " (n = ", nPerClimaticRegion, ")"), 
         ClimaticRegion = as.factor(ClimaticRegion))  %>% 
  filter(nPerClimaticRegion > 10) %>% 
  as.data.table() %>% 
  mutate(SlopeMeanTemp_scaled = as.numeric(scale(SlopeMeanTemp)), 
         SlopeMaxTemp_scaled = as.numeric(scale(SlopeMaxTemp)), 
         SlopeMinTemp_scaled = as.numeric(scale(SlopeMinTemp)), 
         SlopePrec_scaled = as.numeric(scale(SlopePrec)), 
         NitrogenDepo_scaled = as.numeric(scale(NitrogenDepo)), 
         HumanModification_scaled = as.numeric(scale(HumanModification)), 
         PaAge_scaled = as.numeric(scale(PaAge)),
         PaAreaKm2_scaled = as.numeric(scale(PaAreaKm2)))

dtModLong <- dtMod %>% 
  pivot_longer(cols = c("SlopeMeanTemp", "SlopeMaxTemp", "SlopeMinTemp", "SlopePrec",
                        "NitrogenDepo", "HumanModification", "PaAreaKm2", 
                        "HerbivoreBiomassKgKm2", "HerbivoreFunDiv", "HerbivoreSpeciesRichness", 
                        "MaxBodyMass", "MeanBodyMassCwm"), 
               names_to = "cleanVar", values_to = "varValue")

table(dtMod$ClimaticRegion)

## load model results 
foreach.results <- readRDS("builds/modelOutputs/exploratoryGamReserveResPredAndSpecOnly.Rds")

bm.spec <- foreach.results$bm.spec %>% unique() %>% 
  group_by(modelGroup) %>% slice_min(aic)

bestFormulas <- unique(bm.spec$formulaID)

pred <- foreach.results$pred %>% 
  unique() %>% 
  filter(formulaID %in% bestFormulas) %>% 
  left_join(bm.spec) %>% 
  mutate(randomEffect = case_when(
    .default = FALSE,
    cleanVar == "SlopeMeanTemp" & grepl("s\\(ClimaticRegion, SlopeMeanTemp_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "SlopeMaxTemp" & grepl("s\\(ClimaticRegion, SlopeMaxTemp_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "SlopeMinTemp" & grepl("s\\(ClimaticRegion, SlopeMinTemp_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "SlopePrec" & grepl("s\\(ClimaticRegion, SlopePrec_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "NitrogenDepo" & grepl("s\\(ClimaticRegion, NitrogenDepo_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "HumanModification" & grepl("s\\(ClimaticRegion, HumanModification_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "PaAge" & grepl("s\\(ClimaticRegion, PaAge_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "PaAreaKm2" & grepl("s\\(ClimaticRegion, PaAreaKm2_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "HerbivoreBiomassKgKm2" & grepl("s\\(ClimaticRegion, HerbivoreBiomassKgKm2_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "HerbivoreFunDiv" & grepl("s\\(ClimaticRegion, HerbivoreFunDiv_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "HerbivoreSpeciesRichness" & grepl("s\\(ClimaticRegion, HerbivoreSpeciesRichness_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "MaxBodyMass" & grepl("s\\(ClimaticRegion, MaxBodyMass_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
    cleanVar == "MeanBodyMassCwm" & grepl("s\\(ClimaticRegion, MeanBodyMassCwm_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE,
  )) %>% 
  mutate(exclude = ifelse(randomEffect == FALSE & ClimaticRegion != "Arid", "fuckit", "keep")) %>% 
  filter(exclude == "keep")

table(pred$ClimaticRegion)
table(pred$randomEffect)
table(pred$exclude)

### limit the predictor X axis to the 95 % interval 
predictorLimits <- dtModLong %>%
  group_by(cleanVar) %>%
  summarize(
    lowerQuantile = quantile(varValue, 0.025, na.rm = T),
    upperQuantile = quantile(varValue, 0.975, na.rm = T),
  ) 

for(clean.var in unique(predictorLimits$cleanVar)){
  
  upperLim <- predictorLimits[predictorLimits$cleanVar == clean.var, ]$upperQuantile
  lowerLim <- predictorLimits[predictorLimits$cleanVar == clean.var, ]$lowerQuantile
  
  dtModLongSub <- dtModLong[dtModLong$cleanVar == clean.var,] %>% filter(varValue > lowerLim & varValue < upperLim)
  dtModLong <- dtModLong %>% 
    filter(!cleanVar == clean.var) %>% 
    rbind(dtModLongSub)
  
  predSub <- pred[pred$cleanVar == clean.var,] %>% filter(varValue > lowerLim & varValue < upperLim)
  pred <- pred %>% 
    filter(!cleanVar == clean.var) %>% 
    rbind(predSub)
  
  print(clean.var)
  
}
# EVI Trend ---------------------------------------

## Full Non Random
unique(bm.spec[bm.spec$response == "EviTrend", ]$modelGroup)
bm.spec[bm.spec$modelGroup == "EviTrendFullNonRandom",]$r_squared
bm.spec[bm.spec$modelGroup == "EviTrendFullNonRandom",]$aic


pEviFullNonR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = EviTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullNonRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullNonRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("EviTrendFullNonRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("EviTrendFullNonRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("EVI Trend (Full Model Without Random Effects; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "EviTrendFullNonRandom",]$r_squared, 2)) * ")")) +
  
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pEviFullNonR

## Full Non Random No Herbivores 
unique(bm.spec[bm.spec$response == "EviTrend", ]$modelGroup)
bm.spec[bm.spec$modelGroup == "EviTrendNonRandomNoHerbivores",]$r_squared
bm.spec[bm.spec$modelGroup == "EviTrendNonRandomNoHerbivores",]$aic


pEviFullNonRNoH <- ggplot() +
 # geom_point(data = dtModLong, aes(x = varValue, y = EviTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendNonRandomNoHerbivores") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendNonRandomNoHerbivores") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("EviTrendNonRandomNoHerbivores") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("EviTrendNonRandomNoHerbivores") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("EVI Trend (Full Model Without Random Effects; " 
                      ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "EviTrendNonRandomNoHerbivores",]$r_squared, 2)) * ")")) +
  
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pEviFullNonRNoH

## Full Random
bm.spec[bm.spec$modelGroup == "EviTrendFullRandom",]$r_squared
bm.spec[bm.spec$modelGroup == "EviTrendFullRandom",]$formula


pEviFullR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = EviTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  # geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == FALSE, ],
  #             aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  # geom_line(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == FALSE, ], 
  #           aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("EVI Trend (Full Model With Random Effects; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "EviTrendFullRandom",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  )
pEviFullR


## Full Random No Herbivores 
bm.spec[bm.spec$modelGroup == "EviTrendRandomNoHerbivores",]$r_squared

pEviFullRNoH <- ggplot() +
 # geom_point(data = dtModLong, aes(x = varValue, y = EviTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendRandomNoHerbivores") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendRandomNoHerbivores") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("EviTrendRandomNoHerbivores") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("EviTrendRandomNoHerbivores") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("EVI Trend (Full Model With Random Effects and No Herbi Vars; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "EviTrendRandomNoHerbivores",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  )

pEviFullRNoH


## Best
bm.spec[bm.spec$modelGroup == "EviTrendBestModel",]$r_squared

pEviBest <- ggplot() +
  # geom_point(data = dtModLong[dtModLong$cleanVar %in% unique(pred[modelGroup %in% c("EviTrendBestModel"), ]$cleanVar), ],
  #            aes(x = varValue, y = EviTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendBestModel") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendBestModel") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("EviTrendBestModel") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("EviTrendBestModel") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("EVI Trend (Best Model; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "EviTrendBestModel",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.key.height = unit(0.3, "cm"),
        legend.box = "vertical", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  ) +
  guides(
    color = guide_legend(ncol = 3), 
    fill = guide_legend(ncol = 3),
    linetype = guide_legend(nrow = 1)
  )

pEviBest




pEvi <- grid.arrange(pEviFullR, pEviBest, ncol = 1, heights = c(2.5, 1))
ggsave(plot = pEvi,  "builds/plots/exploratoryGams/EviTrendGambm.spec.png", dpi = 600, height = 12, width = 9)


# Npp Trend ---------------------------------------

## Full Non Random
bm.spec[bm.spec$modelGroup == "NppTrendFullNonRandom",]$r_squared

pNppFullNonR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = NppTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("NppTrendFullNonRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("NppTrendFullNonRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("NppTrendFullNonRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("NppTrendFullNonRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("Npp Trend (Full Model Without Random Effects; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "NppTrendFullNonRandom",]$r_squared, 2)) * ")")) +
  
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pNppFullNonR

## Full Non Random No Herbis
unique(bm.spec[bm.spec$response == "NppTrend", ]$modelGroup)
bm.spec[bm.spec$modelGroup == "NppTrendNonRandomNoHerbivores",]$r_squared

pNppFullNonRNoH <- ggplot() +
  #geom_point(data = dtModLong, aes(x = varValue, y = NppTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("NppTrendNonRandomNoHerbivores") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("NppTrendNonRandomNoHerbivores") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("NppTrendNonRandomNoHerbivores") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("NppTrendNonRandomNoHerbivores") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("Npp Trend (Full Model Without Random Effects And Herbivores; "
                      ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "NppTrendNonRandomNoHerbivores",]$r_squared, 2)) * ")")) +
  
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pNppFullNonRNoH

## Full Random
bm.spec[bm.spec$modelGroup == "NppTrendFullRandom",]$r_squared

pNppFullR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = NppTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("NppTrendFullRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("NppTrendFullRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("NppTrendFullRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("NppTrendFullRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("Npp Trend (Full Model With Random Effects; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "NppTrendFullRandom",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  )

pNppFullR

## Best
bm.spec[bm.spec$modelGroup == "NppTrendBestModel",]$r_squared

pNppBest <- ggplot() +
  #  geom_point(data = dtModLong[dtModLong$cleanVar %in% unique(pred[modelGroup %in% c("NppTrendBestModel"), ]$cleanVar), ], 
  #            aes(x = varValue, y = NppTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("NppTrendBestModel") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("NppTrendBestModel") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("NppTrendBestModel") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("NppTrendBestModel") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("Npp Trend (Best Model; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "NppTrendBestModel",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.key.height = unit(0.3, "cm"),
        legend.box = "vertical", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  ) +
  guides(
    color = guide_legend(ncol = 3), 
    fill = guide_legend(ncol = 3),
    linetype = guide_legend(nrow = 1)
  )

pNppBest

pNpp <- grid.arrange(pNppFullR, pNppBest, ncol = 1, heights = c(2.5, 2))
ggsave(plot = pNpp,  "builds/plots/exploratoryGams/NppTrendGambm.spec.png", dpi = 600, height = 12, width = 9)

# BurnedArea Trend ---------------------------------------

## Full Non Random
bm.spec[bm.spec$modelGroup == "BurnedAreaTrendFullNonRandom",]$r_squared

pBurnedAreaFullNonR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = BurnedAreaTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendFullNonRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendFullNonRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendFullNonRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendFullNonRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("BurnedArea Trend (Full Model Without Random Effects; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "BurnedAreaTrendFullNonRandom",]$r_squared, 2)) * ")")) +
  
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pBurnedAreaFullNonR

## Full Non Random No Herbivores 
bm.spec[bm.spec$response == "BurnedAreaTrend",]$modelGroup
bm.spec[bm.spec$modelGroup == "BurnedAreaTrendNonRandomNoHerbivores",]$r_squared

pBurnedAreaFullNonRNoH <- ggplot() +
 # geom_point(data = dtModLong, aes(x = varValue, y = BurnedAreaTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendNonRandomNoHerbivores") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendNonRandomNoHerbivores") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendNonRandomNoHerbivores") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendNonRandomNoHerbivores") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("BurnedArea Trend (Full Model Without Random Effects and Herbivores; " 
                      ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "BurnedAreaTrendNonRandomNoHerbivores",]$r_squared, 2)) * ")")) +
  
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pBurnedAreaFullNonRNoH

## Full Random
bm.spec[bm.spec$response == "BurnedAreaTrend",]$modelGroup
bm.spec[bm.spec$modelGroup == "BurnedAreaTrendFullRandom",]$r_squared

pBurnedAreaFullR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = BurnedAreaTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendFullRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendFullRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendFullRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendFullRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("BurnedArea Trend (Full Model With Random Effects; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "BurnedAreaTrendFullRandom",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  )

pBurnedAreaFullR

## Full Random
bm.spec[bm.spec$response == "BurnedAreaTrend",]$modelGroup
bm.spec[bm.spec$modelGroup == "BurnedAreaTrendFullRandom",]$r_squared

pBurnedAreaFullRNoH <- ggplot() +
 # geom_point(data = dtModLong, aes(x = varValue, y = BurnedAreaTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendRandomNoHerbivores") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendRandomNoHerbivores") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendRandomNoHerbivores") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendRandomNoHerbivores") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("BurnedArea Trend (Full Model With Random Effects But No Herbivores; "
                      ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "BurnedAreaTrendRandomNoHerbivores",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  )

pBurnedAreaFullRNoH

## Best
bm.spec[bm.spec$modelGroup == "BurnedAreaTrendBestModel",]$r_squared
bm.spec[bm.spec$modelGroup == "BurnedAreaTrendBestModel",]$formula

pBurnedAreaBest <- ggplot() +
  # geom_point(data = dtModLong[dtModLong$cleanVar %in% unique(pred[modelGroup %in% c("BurnedAreaTrendBestModel"), ]$cleanVar), ], 
  #            aes(x = varValue, y = BurnedAreaTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendBestModel") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("BurnedAreaTrendBestModel") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendBestModel") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("BurnedAreaTrendBestModel") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("BurnedArea Trend (Best Model; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "BurnedAreaTrendBestModel",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.key.height = unit(0.3, "cm"),
        legend.box = "vertical", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  ) +
  guides(
    color = guide_legend(ncol = 3), 
    fill = guide_legend(ncol = 3),
    linetype = guide_legend(nrow = 1)
  )

pBurnedAreaBest

pBurnedArea <- grid.arrange(pBurnedAreaFullR, pBurnedAreaBest, ncol = 1, heights = c(2.5, 2))
ggsave(plot = pBurnedArea,  "builds/plots/exploratoryGams/BurnedAreaTrendGambm.spec.png", dpi = 600, height = 12, width = 9)

# Sos Trend ---------------------------------------

## Full Non Random
bm.spec[bm.spec$modelGroup == "SOSTrendFullNonRandom",]$r_squared
bm.spec[bm.spec$response == "SOSTrend",]$modelGroup


pSOSFullNonR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = SOSTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendFullNonRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendFullNonRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("SOSTrendFullNonRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("SOSTrendFullNonRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("SOS Trend (Full Model Without Random Effects; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "SOSTrendFullNonRandom",]$r_squared, 2)) * ")")) +
  
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pSOSFullNonR


## Full Non Random
bm.spec[bm.spec$modelGroup == "SOSTrendNonRandomNoHerbivores",]$r_squared
bm.spec[bm.spec$response == "SOSTrend",]$modelGroup


pSOSFullNonRNoH <- ggplot() +
  #geom_point(data = dtModLong, aes(x = varValue, y = SOSTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendNonRandomNoHerbivores") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendNonRandomNoHerbivores") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("SOSTrendNonRandomNoHerbivores") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("SOSTrendNonRandomNoHerbivores") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("SOS Trend (Full Model Without Random Effects And Herbivores; "
                      ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "SOSTrendNonRandomNoHerbivores",]$r_squared, 2)) * ")")) +
  
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pSOSFullNonRNoH


## Full Random
bm.spec[bm.spec$modelGroup == "SOSTrendFullRandom",]$r_squared

pSOSFullR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = SOSTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendFullRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendFullRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("SOSTrendFullRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("SOSTrendFullRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("SOS Trend (Full Model With Random Effects; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "SOSTrendFullRandom",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  )

pSOSFullR


## Full Random No Herbivores 
bm.spec[bm.spec$modelGroup == "SOSTrendRandomNoHerbivores",]$r_squared
bm.spec[bm.spec$response == "SOSTrend",]$modelGroup

pSOSFullRNoH <- ggplot() +
 # geom_point(data = dtModLong, aes(x = varValue, y = SOSTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendRandomNoHerbivores") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendRandomNoHerbivores") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("SOSTrendRandomNoHerbivores") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("SOSTrendRandomNoHerbivores") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("SOS Trend (Full Model With Random Effects; "
                      ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "SOSTrendRandomNoHerbivores",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  )

pSOSFullRNoH

## Best
bm.spec[bm.spec$modelGroup == "SOSTrendBestModel",]$r_squared

pSOSBest <- ggplot() +
  # geom_point(data = dtModLong[dtModLong$cleanVar %in% unique(pred[modelGroup %in% c("SOSTrendBestModel"), ]$cleanVar), ], 
  #            aes(x = varValue, y = SOSTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendBestModel") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("SOSTrendBestModel") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("SOSTrendBestModel") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("SOSTrendBestModel") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = ClimaticRegion, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = bquote("SOS Trend (Best Model; " ~ R^2 == .(round(bm.spec[bm.spec$modelGroup == "SOSTrendBestModel",]$r_squared, 2)) * ")")) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.key.height = unit(0.3, "cm"),
        legend.box = "vertical", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 10), 
        panel.grid = element_blank()
  ) +
  guides(
    color = guide_legend(ncol = 3), 
    fill = guide_legend(ncol = 3),
    linetype = guide_legend(nrow = 1)
  )

pSOSBest

pSOS <- grid.arrange(pSOSFullR, pSOSBest, ncol = 1, heights = c(2.5, 1.7))
ggsave(plot = pSOS,  "builds/plots/exploratoryGams/SOSTrendGambm.spec.png", dpi = 600, height = 12, width = 9)


## make manuscript plots 

pa <- grid.arrange(pEviBest, pNppBest, heights = c(1, 1.6))
pb <- grid.arrange(pBurnedAreaBest, pSOSBest, heights = c(1.8, 1))

ggsave(plot = pa, "builds/plots/BestModelsEviAndNppReserves.png", dpi = 600, height = 9, width = 7.5)

ggsave(plot = pb, "builds/plots/BestModelsBurnedAndSosReserves.png", dpi = 600, height = 8.5, width = 7.5)



paSp <- grid.arrange(pEviFullR, pNppFullR)
pbSp <- grid.arrange(pBurnedAreaFullR, pSOSFullR)

ggsave(plot = paSp, "builds/plots/supplement/FullREviAndNppReserves.png", dpi = 600, height = 11, width = 7.5)

ggsave(plot = pbSp, "builds/plots/supplement/FullRBurnedAndSosReserves.png", dpi = 600, height = 11, width = 7.5)



