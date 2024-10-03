## compare PAs and controls 
library(scico)
library(data.table)
library(tidyverse)
library(gridExtra)


dtControlRaw <- fread("data/processedData/cleanData/paControlsWithCovs.csv") 

ClimaticRegionNControl <- dtControlRaw %>% 
  group_by(ClimaticRegion) %>% 
  summarize(nPerClimaticRegion = n()) 

dtControl <- dtControlRaw %>% 
  dplyr::select(controlFor, unique_id, NppTrend, EviTrend, BurnedAreaTrend, SOSTrend, ClimaticRegion, BurnedAreaMean) %>% 
  mutate(PaOrControl = "Control") %>% 
  pivot_longer(cols = c("NppTrend", "EviTrend", "BurnedAreaTrend", "SOSTrend"), 
               names_to = "TrendName", values_to = "TrendValue") %>% 
  filter(!(TrendName == "BurnedAreaTrend" & BurnedAreaMean == 0)) %>%
  left_join(ClimaticRegionNControl)

dtPaRaw <- fread("data/processedData/cleanData/pasWithCovs.csv") %>% 
  filter(unique_id %in% unique(dtControl$controlFor))

ClimaticRegionNPas <- dtPaRaw %>% 
  group_by(ClimaticRegion) %>% 
  summarize(nPerClimaticRegion = n()) 

dtPa <- dtPaRaw  %>% 
  dplyr::select(unique_id, NppTrend, EviTrend, BurnedAreaTrend, SOSTrend, ClimaticRegion, BurnedAreaMean) %>% 
  mutate(PaOrControl = "Protected Area", controlFor = "Nothing") %>% 
  pivot_longer(cols = c("NppTrend", "EviTrend", "BurnedAreaTrend", "SOSTrend"), 
               names_to = "TrendName", values_to = "TrendValue") %>% 
  filter(!(TrendName == "BurnedAreaTrend" & BurnedAreaMean == 0)) %>% 
  left_join(ClimaticRegionNPas)


dtPlot <- rbind(dtPa, dtControl) %>% 
  filter(nPerClimaticRegion > 10)

table(dtPa$ClimaticRegion)

pa <- ggplot() +
  geom_jitter(data = dtPlot, aes(x = PaOrControl, y = TrendValue, color = ClimaticRegion, fill = ClimaticRegion), alpha = .5) +
 # geom_violin(data = dtPlot, aes(x = PaOrControl, y = TrendValue), alpha = .1) +
  geom_boxplot(data = dtPlot, aes(x = PaOrControl, y = TrendValue), outlier.shape = NA, alpha = .75, size = 1.1)  +
  facet_wrap(~TrendName, scales = "free", ncol = 4) +
  scale_color_scico_d(palette = "bamako") +
  scale_fill_scico_d(palette = "bamako") +
  labs(x = "", y = "Trend Value", subtitle = "a)") +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 0)) +  
  guides(
          color = guide_legend(ncol = 3), 
          fill = guide_legend(ncol = 3),
          linetype = guide_legend(nrow = 1)
        )
pa

pb <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(data = dtPlot[dtPlot$nPerClimaticRegion > 00,], aes(x = PaOrControl, y = TrendValue, color = ClimaticRegion), alpha = .5) +
 # geom_violin(data = dtPlot[dtPlot$nPerClimaticRegion > 50,], aes(x = PaOrControl, y = TrendValue), alpha = .1) +
  geom_boxplot(data = dtPlot[dtPlot$nPerClimaticRegion > 00,], aes(x = PaOrControl, y = TrendValue), outlier.shape = NA, alpha = .75) +
  facet_grid(TrendName~ClimaticRegion, scales = "free", labeller = label_wrap_gen(width = 15)) +
  scale_color_scico_d(palette = "bamako") +
  labs(x = "", y = "Trend Value", subtitle = "b)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
pb 

p <- grid.arrange(pa, pb, heights = c(1, 1.5))
ggsave(plot = p, "builds/plots/comparisonPasAndControls.png", dpi = 600, height = 12, width = 10)
