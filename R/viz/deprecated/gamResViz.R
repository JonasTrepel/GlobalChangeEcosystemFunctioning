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

pasCovsDTRaw <- fread("data/processedData/cleanData/pasWithCovs.csv") %>% 
  rename(PaAreaKm2 = GIS_AREA, 
         PaYear = STATUS_YR) %>% 
  mutate(PaYear = ifelse(PaYear == 0, NA, PaYear), 
         PaAge = 2023-PaYear) %>% 
  filter(complete.cases(HumanModification, NitrogenDepo,
                        SlopeMaxTemp, SlopeMeanTemp, SlopeMinTemp, 
                        SlopePrec, PaAge, PaAreaKm2, EviTrend))

table(pasCovsDTRaw$PaYear)

ClimaticRegionN <- pasCovsDTRaw %>% 
  group_by(ClimaticRegion) %>% 
  summarize(nPerClimaticRegion = n()) 

dtMod <- pasCovsDTRaw %>% 
  left_join(ClimaticRegionN) %>% 
  mutate(ClimaticRegionN = paste0(ClimaticRegion, " (n = ", nPerClimaticRegion, ")"), 
         ClimaticRegion = as.factor(ClimaticRegion))  %>% 
  filter(nPerClimaticRegion > 50) %>% 
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
                        "NitrogenDepo", "HumanModification", "PaAreaKm2"), 
               names_to = "cleanVar", values_to = "varValue")



## load model results 

foreach.results <- readRDS("builds/modelOutputs/exploratoryGamRes.Rds")

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
    cleanVar == "PaAreaKm2" & grepl("s\\(ClimaticRegion, PaAreaKm2_scaled, bs \\= 're', k \\= 4\\)", formula) ~ TRUE
  )) %>% 
  mutate(exclude = ifelse(randomEffect == FALSE & ClimaticRegion != "Cold", "fuckit", "keep")) %>% 
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
bm.spec[bm.spec$modelGroup == "EviTrendFullNonRandom",]$r_squared

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

## Full Random
bm.spec[bm.spec$modelGroup == "EviTrendFullRandom",]$r_squared

pEviFullR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = EviTrend, color = ClimaticRegion), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = ClimaticRegion), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
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

pEvi <- grid.arrange(pEviFullR, pEviBest, ncol = 1, heights = c(1, 1.4))
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

pNpp <- grid.arrange(pNppFullR, pNppBest, ncol = 1, heights = c(1, 1.4))
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

## Full Random
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

pBurnedArea <- grid.arrange(pBurnedAreaFullR, pBurnedAreaBest, ncol = 1, heights = c(1, 1.4))
ggsave(plot = pBurnedArea,  "builds/plots/exploratoryGams/BurnedAreaTrendGambm.spec.png", dpi = 600, height = 12, width = 9)

# Sos Trend ---------------------------------------

## Full Non Random
bm.spec[bm.spec$modelGroup == "SOSTrendFullNonRandom",]$r_squared

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

pSOS <- grid.arrange(pSOSFullR, pSOSBest, ncol = 1, heights = c(1, 1.4))
ggsave(plot = pSOS,  "builds/plots/exploratoryGams/SOSTrendGambm.spec.png", dpi = 600, height = 12, width = 9)


## make manuscript plots 

pa <- grid.arrange(pEviBest, pNppBest)
pb <- grid.arrange(pBurnedAreaBest, pSOSBest)

ggsave(plot = pa, "builds/plots/BestModelsEviAndNpp.png", dpi = 600, height = 11, width = 7.5)

ggsave(plot = pb, "builds/plots/BestModelsBurnedAndSos.png", dpi = 600, height = 11, width = 7.5)



paSp <- grid.arrange(pEviFullR, pNppFullR)
pbSp <- grid.arrange(pBurnedAreaFullR, pSOSFullR)

ggsave(plot = paSp, "builds/plots/supplement/FullREviAndNpp.png", dpi = 600, height = 9, width = 7.5)

ggsave(plot = pbSp, "builds/plots/supplement/FullRBurnedAndSos.png", dpi = 600, height = 9, width = 7.5)



