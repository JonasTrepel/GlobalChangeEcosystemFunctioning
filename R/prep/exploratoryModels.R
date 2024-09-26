library(data.table)
library(tidyverse)
library(ggcorrplot)
library(sf)
library(caret)
library(caretEnsemble)
library(pdp)
library(gridExtra)
library(scales)
library(grid)

dt <- fread('data/processedData/cleanData/gridWithCovs25kSample.csv') %>% 
  mutate(BodyMassLoss = ifelse(is.infinite(BodyMassLoss), NA, BodyMassLoss), 
         iucnCatOrd = ifelse(is.na(iucnCatOrd), 9, iucnCatOrd))
  
summary(dt)
table(dt$iucnCatOrd)
dtCont <- dt %>%
  dplyr::select(-c(Biome, iucnCat, LandCover, gridID, sovereignt, continent)) %>% 
  filter(complete.cases(.))

corr <- round(cor(dtCont), 1)
ggcorrplot(corr, type = 'lower', lab = T)
?ggcorrplot

table(dt$Biome)
table(dt$iucnCat)


##### DATA SUBSETS #####
## Biomes w over 1000 cells aka Biomes which we'll analyse separately: 
# Deserts & Xeric Shrublands
# Tropical & Subtropical Grasslands, Savannas & Shrublands
# Tundra
# Boreal Forests/Taiga
# Temperate Broadleaf & Mixed Forests
# Temperate Grasslands, Savannas & Shrublands
# Tropical & Subtropical Moist Broadleaf Forests

## IUCN Cat Ia, Ib and II
## vs Not protected

##### Responses #####
#EviTrend
#NppTrend
#FireFreqTrend
#BurnedAreaTrend
#SOSTrend
#EviSDTrend

##### Predictors #####
#HumanModification
#SpeciesLoss
#BodyMassLoss
#AbsMeanTempDiff
#RelMeanTempDiff
#AbsMaxTempDiff
#AbsMinTempDiff
#RelMinTempDiff
#AbsPrecDiff
#RelPrecDiff
#SlopeMeanTemp
#SlopeMaxTemp
#SlopeMinTemp
#SlopePrec
#NitrogenDepo

##### Requires more thinking: 
#BurnedAreaMean
#FireFrequencyMean
#EviSD
# Current climate /Mean climate

################################################


##############################################################################            
################################## CREATE MODEL GUIDE ########################         
##############################################################################    


subsets <- c("!is.na(gridID)",
             "Biome == 'Deserts & Xeric Shrublands'",
             "Biome == 'Tropical & Subtropical Grasslands, Savannas & Shrublands'", 
             "Biome == 'Tundra'",
             "Biome == 'Boreal Forests/Taiga'",
             "Biome == 'Temperate Broadleaf & Mixed Forests'", 
             "Biome == 'Temperate Grasslands, Savannas & Shrublands'",
             "Biome == 'Tropical & Subtropical Moist Broadleaf Forests'")

tier_labels<- c('Full model',
                'Deserts & Xeric Shrublands',
                'Tropical & Subtropical Grasslands, Savannas & Shrublands',
                'Tundra',
                'Boreal Forests/Taiga',
                'Temperate Broadleaf & Mixed Forests', 
                'Temperate Grasslands, Savannas & Shrublands',
                'Tropical & Subtropical Moist Broadleaf Forests')

tiers <- c('main', 'deserts', 'tropical_open', 'tundra', 'boreal_forest', 'temperate_broadleaf', 
           'temperate_open', 'tropical_moist_broadleaf')

responses <- c('EviTrend', 'NppTrend', 'FireFreqTrend', 'BurnedAreaTrend', 'SOSTrend', 'EviSDTrend')

dt.tier.raw <- data.table(
  subset = subsets, 
  tier = tiers, 
  tier_label = tier_labels, 
  n = 0)

tier.resp <- CJ(response = responses, tier = tiers) %>% 
  mutate(terms = c("'HumanModification', 'AbsMeanTempDiff', 'RelMeanTempDiff', 'AbsMaxTempDiff', 'AbsMinTempDiff', 'RelMinTempDiff', 'AbsPrecDiff', 'RelPrecDiff', 'SlopeMeanTemp', 'SlopeMaxTemp', 'SlopeMinTemp', 'SlopePrec', 'NitrogenDepo'")) %>%
  mutate(
    response_tier = paste0(response, '_', tier)
  )

dt.tier <- dt.tier.raw %>% left_join(tier.resp)

for(i in 1:nrow(dt.tier)){
  subset <- dt.tier[i, ]$subset
  dt.sub <- dt %>% dplyr::filter(eval(parse(text = subset)))
  
  nr <- nrow(dt.sub)
  dt.tier[i, ]$n <- nr
  
}
unique(dt.tier$n)



# color palettes 
palette <- c('HumanModification' = "#fab255",
             'SpeciesLoss' = "#43b284",
             'BodyMassLoss' = "#43b284",
             'AbsMeanTempDiff' = "#dd5129",
             'RelMeanTempDiff' = "#dd5129",
             'AbsMaxTempDiff' = "#dd5129",
             'AbsMinTempDiff' = "#dd5129",
             'RelMinTempDiff' = "#dd5129",
             'AbsPrecDiff' = "#0f7ba2",
             'RelPrecDiff' = "#0f7ba2",
             'SlopeMeanTemp' = "#dd5129",
             'SlopeMaxTemp' = "#dd5129",
             'SlopeMinTemp' = "#dd5129",
             'SlopePrec' = "#0f7ba2",
             'NitrogenDepo' = "#fab255")

#c(met.brewer(name = "Egypt"))
palette.methods <- c("gbm" = "#dd5129", "rf" = "#0f7ba2", "xgbTree" = "#43b284", "ensemble" = "#fab255")

dt.names <- data.table(
  term = c('HumanModification',
           'SpeciesLoss',
           'BodyMassLoss',
           'AbsMeanTempDiff',
           'RelMeanTempDiff',
           'AbsMaxTempDiff',
           'AbsMinTempDiff',
           'RelMinTempDiff',
           'AbsPrecDiff',
           'RelPrecDiff',
           'SlopeMeanTemp',
           'SlopeMaxTemp',
           'SlopeMinTemp',
           'SlopePrec',
           'NitrogenDepo'),
  clean_term = c('HumanModification',
                 'SpeciesLoss',
                 'BodyMassLoss',
                 'AbsMeanTempDiff',
                 'RelMeanTempDiff',
                 'AbsMaxTempDiff',
                 'AbsMinTempDiff',
                 'RelMinTempDiff',
                 'AbsPrecDiff',
                 'RelPrecDiff',
                 'SlopeMeanTemp',
                 'SlopeMaxTemp',
                 'SlopeMinTemp',
                 'SlopePrec',
                 'NitrogenDepo'))

### create Tune Grids 

##### CHANGED IT
### GBM 
tuneGridGbm <- expand.grid(
  shrinkage = c( .001, .005, 0.01),
  interaction.depth = c(2,3,5),
  n.minobsinnode = c(20, 40, 60), 
  n.trees = seq(1000, 10000, 250)
)

### Random forest 
tuneGridRf <- expand.grid(
  mtry = c(1:10))

### XGBoost

tuneGridXgbTree <- expand.grid(
  nrounds = seq(50, 200, 50), 
  eta = c(0.01, 0.1, 0.3),
  gamma = 0,
  max_depth = c(1, 3, 5, 7),
  min_child_weight = c(5, 10, 15), 
  colsample_bytree = 0.9, 
  subsample = 0.8
)

### Fit Control 

fitControl <- trainControl(## 10 fold cross validation
  method = "repeatedcv",
  number = 5, #number of splits 
  repeats = 5, #repeat 10 times
  savePredictions = "final", #keep final model predictions, would otherwise be dumped
  returnResamp = "final")


#for(i in 1:nrow(dt.tier)){

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)
library(parallel)

# Create and register a cluster
clust <- makeCluster(48)
registerDoSNOW(clust)

## progress bar 
iterations <- nrow(dt.tier)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic()
foreach.results <- foreach(i = 1:nrow(dt.tier),
                           .packages = c('data.table', 'tidyverse',
                                         'sf', 'caret',
                                         'caretEnsemble', 'pdp', 
                                         'gridExtra', 'scales', 
                                         'grid', 'Metrics'),
                           .options.snow = opts,
                           .inorder = FALSE,
                           .combine = rbind) %dopar% {
  
  tier <- dt.tier[i, ]$tier
  subset <- dt.tier[i, ]$subset
  tier_label <- dt.tier[i, ]$tier_label
  response <- dt.tier[i, ]$response
  response_tier <- dt.tier[i, ]$response_tier
  
  terms <- dt.tier[i, ]$terms
  
  terms.vector <- strsplit(gsub("'", "", terms), ", ")[[1]]
  
  
  dt.sub <- dt %>% dplyr::filter(eval(parse(text = subset)))
  
  print(paste0("starting ", response,"; tier: ", tier, " (", i, "/", nrow(dt.tier), ")"))
  
  dt.gbm <- dt.sub %>% 
    dplyr::select(all_of(response), all_of(terms.vector)) %>% dplyr::filter(complete.cases(.)) #%>% sample_n(50)
  
  ##############################################################################            
  ################################## RUN ALL MODELS ############################            
  ##############################################################################            
  if(sum(dt.gbm[[response]]) == 0){return(NULL)}
  
  formula <- as.formula(paste(response, "~ ."))
  
  gbmFit <- train(formula,
                  data = dt.gbm, 
                  method = "gbm",
                  trControl = fitControl,
                  verbose = FALSE,
                  tuneGrid = tuneGridGbm)
  
  head(gbmFit$results[order(gbmFit$results$RMSE),])#sort that the best model (lowest root mean squared error) is on top
  
  
  stats <- head(gbmFit$results[order(gbmFit$results$RMSE),], 1)
  
  #### Variable Importance 
  
  relVarImp <- as.data.table(summary(gbmFit)) %>% mutate(var = fct_reorder(var, rel.inf), 
                                                         term = var)
  
  ### prepare plot varImp
  
  p.var <- ggplot(data = relVarImp) +
    geom_col(aes(y = var, x = rel.inf, fill = var),
             position = position_dodge(width = 0.8), alpha = 0.9) +
    scale_fill_manual(values = palette) +
    labs(y = "", x = "Relative Variable Importance") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12))
  
  p.var
  
  #### Partial Dependance Plot
  
  dt.gbm.P <- dt.gbm %>% dplyr::select(-all_of(response))
  
  

    for(p in 1:ncol(dt.gbm.P)){
      predTmp <- partial(gbmFit, p, train = dt.gbm.P)
      term_name <- colnames(predTmp)[1]
      colnames(predTmp) <- c("x", "y")
      predTmp <- predTmp %>% mutate(term = paste0(term_name),
                                    method = "GBM")
      if(p==1){
        marg <- predTmp}else{
          marg <- rbind(marg, predTmp)} }
  
  marg.plot <- marg %>% 
    mutate(
      var = gsub("`", "", term)) %>% 
    left_join(relVarImp)
  

  marg.plot <- marg.plot[!grepl("spatial_predictor", marg.plot$term),]
  
  
  dt.points <- dt.gbm %>% pivot_longer(
    cols = c(-{{response}}), 
    values_to = "x", names_to = "term") %>% 
    mutate(y = median(marg.plot$y, na.rm = T)) %>% 
    left_join(relVarImp) 
  
  dt.points <- dt.points[!grepl("spatial_predictor", dt.points$term),]
  
  
  
  dt.mean.rug <- dt.points %>% dplyr::select(term, x, y) %>% left_join(dt.names)
  
  rects <- dt.mean.rug %>%
    group_by(term) %>%
    mutate(
      lower_quantile_x = quantile(x, 0.05),
      upper_quantile_x = quantile(x, 0.95),
    ) %>%
    ungroup() %>% 
    group_by(term) %>%
    summarize(
      ymin = -Inf,
      ymax = Inf,
      xmin1 = -Inf,
      xmax1 = first(lower_quantile_x),
      xmin2 = first(upper_quantile_x),
      xmax2 = Inf
    ) %>%
    ungroup() %>% left_join(dt.names)
  
  ordered_terms_pred <- relVarImp$var
  marg.plot$clean_term <- factor(marg.plot$term, levels = ordered_terms_pred)
  
  clean.label <- case_when(
    .default = response, 
    response == "A" ~ "A")
  
  p.pd.final <- ggplot()+
    geom_line(data = marg.plot, aes(x=x, y=y, color = clean_term), linewidth = 1.1) +
    #geom_line(data = predsBtPlot, aes(x=x, y=y, group = iteration, color = clean_term), alpha = 0.15, linewidth = 0.5, color = "grey") +
    geom_line(data = marg.plot, aes(x=x, y=y, color = clean_term), linewidth = 1.1) +
    facet_wrap(~factor(clean_term),
               scales="free_x", ncol = 4) +
    scale_x_continuous(breaks = extended_breaks(n = 3)) +
    scale_color_manual(values = palette) +
    geom_rect(data = rects, aes(xmin = xmin1, xmax = xmax1, ymin = ymin, ymax = ymax), 
               fill = "white", alpha = 0.8, inherit.aes = FALSE) +
     geom_rect(data = rects, aes(xmin = xmin2, xmax = xmax2, ymin = ymin, ymax = ymax), 
               fill = "white", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(data = dt.mean.rug, aes(x = x, y = y), sides="b", length = unit(0.03, "npc"), outside = TRUE) +
    coord_cartesian(clip = "off") +
    ylim(min(marg.plot$y, na.rm =T), max(marg.plot$y, na.rm = T)) +
    theme_classic() +
    labs(y = paste0(clean.label), x = "") +
    theme(legend.position = "none",
          panel.grid = element_line(color = "white"), 
          axis.text.y = element_text(size = 12), 
          axis.text.x = element_text(size = 10), 
          axis.title = element_text(size = 12), 
          axis.ticks = element_blank(), 
          panel.spacing = unit(0.6, "lines"), 
          strip.background = element_rect(color = "white", fill= "grey95") 
    )
  p.pd.final
  
  p.comb <- grid.arrange(p.pd.final, p.var, ncol = 2, widths = c(1.5, 1), 
                         top = textGrob(paste0(clean.label, "\n",
                                               tier_label, " (n = ", dt.tier[i, ]$n, "; R-sq = ", round(stats$Rsquared, 2), "; RMSE = ",
                                               round(stats$RMSE, 2),")") ,gp=gpar(fontsize=12)))
  
  
  filename <- paste0("builds/plots/exploratoryGBMs/comb_plot_", response_tier, ".png")
  ggsave(plot = p.comb, filename = filename, dpi = 600, height = 8, width = 14)
  
}

stopCluster(clust)
print("loop done")

toc()
stop(pb)
