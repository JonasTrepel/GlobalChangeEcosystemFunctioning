library(data.table)
library(tidyverse)


dt_median_evi <- fread("data/processedData/dataFragments/grid_evi_trends.csv") %>% 
  rename(trend_median_evi = evi_coef)
                       

dt_mean_evi <- fread("data/processedData/dataFragments/grid_envi_trends.csv") %>% 
  rename(trend_mean_evi = envi_coef)

dt <- left_join(dt_median_evi, dt_mean_evi)

cor.test(dt$trend_median_evi, dt$trend_mean_evi)

plot <- dt %>%
  ggplot() +
  geom_abline(linetype = "dashed") +
  geom_point(aes(x = trend_median_evi, y = trend_mean_evi), alpha = 0.75) +
  annotate("text", x = -100, y = 100, label = "cor = 0.94; p < 0.0001") +
  labs(x = "Trend of Median EVI", y = "Trend of Mean EVI") +
  theme_classic()
plot

ggsave(plot = plot, "builds/plots/mean_vs_median_evi.png", dpi = 600, height = 8, width = 8)
