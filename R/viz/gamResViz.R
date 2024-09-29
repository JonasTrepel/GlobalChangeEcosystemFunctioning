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

biomeN <- pasCovsDTRaw %>% 
  group_by(Biome) %>% 
  summarize(nPerBiome = n()) 

dtMod <- pasCovsDTRaw %>% 
  left_join(biomeN) %>% 
  mutate(BiomeN = paste0(Biome, " (n = ", nPerBiome, ")"), 
         Biome = as.factor(Biome))  %>% 
  filter(nPerBiome > 50) %>% 
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
                        "NitrogenDepo", "HumanModification", "PaAge","PaAreaKm2"), 
               names_to = "cleanVar", values_to = "varValue")

## load model results 

foreach.results <- readRDS("builds/modelOutputs/exploratoryGamRes.Rds")

res <- foreach.results$res %>% unique() %>% 
  group_by(modelGroup) %>% slice_min(aic)

bestFormulas <- unique(res$formulaID)

pred <- foreach.results$pred %>% 
  unique() %>% 
  filter(formulaID %in% bestFormulas) %>% 
  left_join(res) %>% 
  mutate(randomEffect = case_when(
    .default = FALSE,
    cleanVar == "SlopeMeanTemp" & grepl("s\\(Biome, SlopeMeanTemp_scaled, bs \\= 're'\\)", formula) ~ TRUE,
    cleanVar == "SlopeMaxTemp" & grepl("s\\(Biome, SlopeMaxTemp_scaled, bs \\= 're'\\)", formula) ~ TRUE,
    cleanVar == "SlopeMinTemp" & grepl("s\\(Biome, SlopeMinTemp_scaled, bs \\= 're'\\)", formula) ~ TRUE,
    cleanVar == "SlopePrec" & grepl("s\\(Biome, SlopePrec_scaled, bs \\= 're'\\)", formula) ~ TRUE,
    cleanVar == "NitrogenDepo" & grepl("s\\(Biome, NitrogenDepo_scaled, bs \\= 're'\\)", formula) ~ TRUE,
    cleanVar == "HumanModification" & grepl("s\\(Biome, HumanModification_scaled, bs \\= 're'\\)", formula) ~ TRUE,
    cleanVar == "PaAge" & grepl("s\\(Biome, PaAge_scaled, bs \\= 're'\\)", formula) ~ TRUE,
    cleanVar == "PaAreaKm2" & grepl("s\\(Biome, PaAreaKm2_scaled, bs \\= 're'\\)", formula) ~ TRUE
  ))

table(pred$randomEffect)
# EVI Trend ---------------------------------------

## Full Non Random
res[res$modelGroup == "EviTrendFullNonRandom",]$r_squared

pEviFullNonR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = EviTrend, color = Biome), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullNonRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = Biome), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullNonRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("EviTrendFullNonRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("EviTrendFullNonRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = Biome, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = "EVI Trend (Full Model Without Random Effects)", 
       subtitle = bquote(R^2 == .(round(res[res$modelGroup == "EviTrendFullNonRandom",]$r_squared, 2)))) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pEviFullNonR

## Full Random
res[res$modelGroup == "EviTrendFullRandom",]$r_squared

pEviFullR <- ggplot() +
  geom_point(data = dtModLong, aes(x = varValue, y = EviTrend, color = Biome), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = Biome), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("EviTrendFullRandom") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = Biome, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = "EVI Trend (Full Model With Random Effects)", 
       subtitle = bquote(R^2 == .(round(res[res$modelGroup == "EviTrendFullRandom",]$r_squared, 2)))) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pEviFullR

## Best
res[res$modelGroup == "EviTrendBestModel",]$r_squared

pEviBest <- ggplot() +
  geom_point(data = dtModLong[dtModLong$cleanVar %in% unique(pred[modelGroup %in% c("EviTrendBestModel"), ]$cleanVar), ], aes(x = varValue, y = EviTrend, color = Biome), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.1) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendBestModel") & randomEffect == TRUE, ], 
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub, fill = Biome), alpha = 0.5) +
  geom_ribbon(data = pred[modelGroup %in% c("EviTrendBestModel") & randomEffect == FALSE, ],
              aes(x = varValue, ymin =  ci.lb, ymax = ci.ub), alpha = 0.5, fill = "grey50") +
  geom_line(data = pred[modelGroup %in% c("EviTrendBestModel") & randomEffect == FALSE, ], 
            aes(x = varValue, y = fit, linetype = significance), color = "black", linewidth = 1.15) +
  geom_line(data = pred[modelGroup %in% c("EviTrendBestModel") & randomEffect == TRUE, ],
            aes(x = varValue, y = fit, color = Biome, linetype = significance),linewidth = 1.15) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  scale_linetype_manual(values = c("significant" = "solid", "non significant" = "dotted")) + 
  facet_wrap(~ cleanVar, scales = "free_x", ncol = 4) +
  labs(title = "EVI Trend (Best Model)", 
       subtitle = bquote(R^2 == .(round(res[res$modelGroup == "EviTrendBestModel",]$r_squared, 2)))) +
  theme_bw() +
  theme(
        legend.position = c(0.85, 0.2), 
        legend.justification = c(1, 0),
        plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        panel.grid = element_blank()
  )

pEviBest

pEvi <- grid.arrange(pEviFullNonR, pEviFullR, pEviBest, ncol = 1)
ggsave(plot = pEvi,  "builds/plots/exploratoryGams/EviTrendGamRes.png", dpi = 600, height = 16, width = 10)
