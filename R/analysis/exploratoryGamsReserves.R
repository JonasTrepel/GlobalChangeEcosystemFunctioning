library(data.table)
library(tidyverse)
library(ggridges)
library(scico)
library(mgcv)
library(MuMIn)
library(broom)
library(tictoc)

source("R/functions/partialPred.R")



pasCovsDTRaw <- fread("data/processedData/cleanData/saReservesWithCovs.csv") %>% 
  rename(PaAreaKm2 = area_ha) %>% 
  filter(complete.cases(HumanModification, NitrogenDepo,
                        SlopeMaxTemp, SlopeMeanTemp, SlopeMinTemp, 
                        SlopePrec, PaAreaKm2, EviTrend))

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
         PaAreaKm2_scaled = as.numeric(scale(PaAreaKm2)),
         MeanBodyMassCwm_scaled = as.numeric(scale(MeanBodyMassCwm)),
         MaxBodyMass_scaled = as.numeric(scale(MaxBodyMass)),
         HerbivoreSpeciesRichness_scaled = as.numeric(scale(HerbivoreSpeciesRichness)),
         HerbivoreFunDiv_scaled = as.numeric(scale(HerbivoreFunDiv)),
         HerbivoreBiomassKgKm2_scaled = as.numeric(scale(HerbivoreBiomassKgKm2)))

dtCorr <- dtMod %>%
  dplyr::select(SlopeMeanTemp, SlopeMaxTemp, SlopeMinTemp, SlopePrec, 
                NitrogenDepo, HumanModification, PaAge, PaAreaKm2, 
                MeanBodyMassCwm, MaxBodyMass, HerbivoreSpeciesRichness, 
                HerbivoreFunDiv, HerbivoreBiomassKgKm2) %>% filter(complete.cases(.))
corr <- round(cor(dtCorr), 1)
ggcorrplot::ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE)
#Build model guide 

responses<- c("BurnedAreaTrend",
              "EviTrend",
              "NppTrend",
              "SOSTrend")


names(dtMod %>% dplyr::select(contains("scaled")))

vars <- c("s(Biome, SlopeMeanTemp_scaled, bs = 're')",
          "s(Biome, SlopeMaxTemp_scaled, bs = 're')",
          "s(Biome, SlopeMinTemp_scaled, bs = 're')",
          "s(Biome, SlopePrec_scaled, bs = 're')", 
          "s(Biome, NitrogenDepo_scaled, bs = 're')", 
          "s(Biome, HumanModification_scaled, bs = 're')", 
          "s(Biome, PaAreaKm2_scaled, bs = 're')",
          "s(Biome, MaxBodyMass_scaled, bs = 're')",
          "s(Biome, MeanBodyMassCwm_scaled, bs = 're')",
          "s(Biome, HerbivoreSpeciesRichness_scaled, bs = 're')",
          "s(Biome, HerbivoreFunDiv_scaled, bs = 're')",
          "s(Biome, HerbivoreBiomassKgKm2_scaled, bs = 're')",
          "s(PaAreaKm2_scaled)",
          "s(SlopeMeanTemp_scaled)",
          "s(SlopeMaxTemp_scaled)",
          "s(SlopeMinTemp_scaled)",
          "s(SlopePrec_scaled)", 
          "s(NitrogenDepo_scaled)", 
          "s(HumanModification_scaled)", 
          "s(MaxBodyMass_scaled)",
          "s(MeanBodyMassCwm_scaled)",
          "s(HerbivoreSpeciesRichness_scaled)",
          "s(HerbivoreFunDiv_scaled)",
          "s(HerbivoreBiomassKgKm2_scaled)"
)


# Function to generate all combinations of the variables
generateCombinations <- function(vars, max.n) {
  all.combs <- c()  # Initialize an empty vector to store combinations
  
  # Loop through different sizes of combinations 
  for (i in 1:max.n) {
    comb <- combn(vars, i, simplify = FALSE)  # Generate combinations of size i
    all.combs <- c(all.combs, comb)  # Append to all_combinations
  }
  
  return(all.combs)
}

# Generate all combinations
combinations <- generateCombinations(vars, max.n = 12)

c.comb <- sapply(combinations, function(x) paste(x, collapse = " + "))


## build guides
guideRawRaw <- CJ(vars = c.comb, 
                  response = responses) %>% 
  mutate(formulaID = paste0("formula_", 1:nrow(.)), 
         nVar =  sapply(vars, function(x) length(unlist(strsplit(x, " \\+ ")))), 
         formula = paste0(response, " ~ ", vars)) 

#inter.var <- c(1)
# guideInter <-  CJ(vars = inter.var, 
#                    response = responses) %>% 
#   mutate(formulaID = paste0("intercept_formula_", 1:nrow(.)), 
#          nVar = 0, 
#          formula = paste0(response, " ~ ", vars))

guideRaw<- guideRawRaw
#guideRaw <- rbind(guideRawRaw, guideInter)

#### check correlations ---------------
vars
vars.clean <- gsub("s\\(", "", vars)
vars.clean <- gsub("\\)", "", vars.clean)
vars.clean <- gsub(", k = 3", "", vars.clean)
vars.clean <- gsub("Biome, ", "", vars.clean)
vars.clean <- gsub("\\, bs = 're'", "", vars.clean)

vars.clean <- data.table(
  vars.clean = vars.clean) %>%
  filter(!grepl("Biome", vars.clean)) %>%
  pull() %>%
  unique()
vars.clean

pair.combs <- combn(vars.clean, 2, simplify = FALSE)  

dt.corr <- data.table(
  var1 = character(), 
  var2 = character(), 
  corr = numeric(), 
  exclude_if = character()
)

exclusions <- c()

for(i in 1:length(pair.combs)){
  
  dt.tmp <- data.table(
    var1 = pair.combs[[i]][1], 
    var2 = pair.combs[[i]][2], 
    corr = NA,
    exclude_if = NA
  )
  
  var1 <- dtMod %>% dplyr::select(all_of(dt.tmp$var1)) %>% pull()
  var2 <- dtMod %>% dplyr::select(all_of(dt.tmp$var2)) %>% pull()
  
  cor.ob <- cor.test(var1, var2)
  
  corr <- unname(cor.ob$estimate)
  
  dt.tmp$corr <- corr
  
  exclusions.tmp <- c()
  
  if(abs(corr) >= 0.6){
    
    dt.tmp$exclude_if <- paste0("grepl('", dt.tmp$var1, "', formula) & grepl('", dt.tmp$var2,"', formula)")
    
    exclusions.tmp <- guideRaw %>% filter(eval(parse(text = dt.tmp$exclude_if))) %>% dplyr::select(formulaID) %>% pull()
    
  }
  
  dt.corr <- rbind(dt.corr, dt.tmp)
  
  exclusions <- c(exclusions, exclusions.tmp)
  
  exclusions <- unique(exclusions)
  
  print(i)
}

guideRaw2 <- guideRaw %>% filter(!formulaID %in% c(exclusions)) 


#build separate guides for the best model, the full model without random effects and the full model with only random effects (except PA size and age)
subGuideBest <- guideRaw2 %>% 
  mutate(interaction = ifelse(grepl("by", vars), TRUE, FALSE), 
         exclude = case_when(
           .default = "no", 
           grepl("s(NitrogenDepo_scaled)", formula) & grepl("s(Biome, NitrogenDepo_scaled, bs = 're')", formula) ~ "exclude",
           grepl("s(SlopeMeanTemp_scaled)", formula) & grepl("s(Biome, SlopeMeanTemp_scaled, bs = 're')", formula) ~ "exclude",
           grepl("s(SlopeMaxTemp_scaled)", formula) & grepl("s(Biome, SlopeMaxTemp_scaled, bs = 're')", formula) ~ "exclude",
           grepl("s(SlopeMinTemp_scaled)", formula) & grepl("s(Biome, SlopeMinTemp_scaled, bs = 're')", formula) ~ "exclude",
           grepl("s(SlopePrec_scaled)", formula) & grepl("s(Biome, SlopePrec_scaled, bs = 're')", formula) ~ "exclude",
           grepl("s(HumanModification_scaled)", formula) & grepl("s(Biome, HumanModification_scaled, bs = 're')", formula) ~ "exclude", 
           grepl("s(PaAge_scaled)", formula) & grepl("s(Biome, PaAge_scaled, bs = 're')", formula) ~ "exclude", 
           grepl("s(PaAreaKm2_scaled)", formula) & grepl("s(Biome, PaAreaKm2_scaled, bs = 're')", formula) ~ "exclude")
  ) %>% 
  filter(exclude == "no") %>% mutate(modelGroup = paste0(response, "BestModel")) %>% 
  mutate(randomEffect = ifelse(grepl("'re'", vars), TRUE,FALSE)) 


subGuideFullNonRandom <- guideRaw %>% 
  mutate(interaction = ifelse(grepl("by", vars), TRUE, FALSE), 
         exclude = case_when(
           .default = "no", 
           grepl("bs = 're')", formula) ~ "exclude",
           nVar < 12 ~ "exclude")) %>% 
  filter(exclude == "no") %>% mutate(modelGroup = paste0(response, "FullNonRandom")) %>% 
  mutate(randomEffect = ifelse(grepl("'re'", vars), TRUE,FALSE)) 

subGuideFullNonRandomNoHerbivores <- guideRaw %>% 
  mutate(interaction = ifelse(grepl("by", vars), TRUE, FALSE), 
         exclude = case_when(
           .default = "no", 
           grepl("bs = 're')", formula) ~ "exclude",
           grepl("MaxBodyMass_scaled", formula) ~ "exclude",
           grepl("MeanBodyMassCwm_scaled", formula) ~ "exclude",
           grepl("HerbivoreSpeciesRichness_scaled", formula) ~ "exclude",
           grepl("HerbivoreFunDiv_scaled", formula) ~ "exclude",
           grepl("HerbivoreBiomassKgKm2_scaled", formula) ~ "exclude",
           nVar < 7 ~ "exclude")) %>% 
  filter(exclude == "no") %>% mutate(modelGroup = paste0(response, "NonRandomNoHerbivores")) %>% 
  mutate(randomEffect = ifelse(grepl("'re'", vars), TRUE,FALSE)) 


subGuideFullRandom <- guideRaw2 %>%
  mutate(interaction = ifelse(grepl("by", vars), TRUE, FALSE),
         exclude = case_when(
           .default = "no",
           grepl("NitrogenDepo_scaled)", formula) ~ "exclude",
           grepl("SlopeMeanTemp_scaled)", formula) ~ "exclude",
           grepl("SlopeMaxTemp_scaled)", formula) ~ "exclude",
           grepl("SlopeMinTemp_scaled)", formula) ~ "exclude",
           grepl("SlopePrec_scaled)", formula) ~ "exclude",
           grepl("PaAge_scaled)", formula) ~ "exclude",
           grepl("PaAreaKm2_scaled)", formula) ~ "exclude",
           grepl("HumanModification_scaled)", formula) ~ "exclude",
           grepl("MaxBodyMass_scaled)", formula) ~ "exclude",
           grepl("MeanBodyMassCwm_scaled)", formula) ~ "exclude",
           grepl("HerbivoreSpeciesRichness_scaled)", formula) ~ "exclude",
           grepl("HerbivoreFunDiv_scaled)", formula) ~ "exclude",
           grepl("HerbivoreBiomassKgKm2_scaled)", formula) ~ "exclude",
           nVar < 12 ~ "exclude")) %>%
  filter(exclude == "no") %>%
  mutate(modelGroup = paste0(response, "FullRandom")) %>%
  mutate(randomEffect = ifelse(grepl("'re'", vars), TRUE, FALSE))


subGuideRandomNoHerbivores <- guideRaw2 %>%
  mutate(interaction = ifelse(grepl("by", vars), TRUE, FALSE),
         exclude = case_when(
           .default = "no",
           grepl("NitrogenDepo_scaled)", formula) ~ "exclude",
           grepl("SlopeMeanTemp_scaled)", formula) ~ "exclude",
           grepl("SlopeMaxTemp_scaled)", formula) ~ "exclude",
           grepl("SlopeMinTemp_scaled)", formula) ~ "exclude",
           grepl("SlopePrec_scaled)", formula) ~ "exclude",
           grepl("PaAge_scaled)", formula) ~ "exclude",
           grepl("PaAreaKm2_scaled)", formula) ~ "exclude",
           grepl("HumanModification_scaled)", formula) ~ "exclude",
           grepl("MaxBodyMass_scaled", formula) ~ "exclude",
           grepl("MeanBodyMassCwm_scaled", formula) ~ "exclude",
           grepl("HerbivoreSpeciesRichness_scaled", formula) ~ "exclude",
           grepl("HerbivoreFunDiv_scaled", formula) ~ "exclude",
           grepl("HerbivoreBiomassKgKm2_scaled", formula) ~ "exclude",
           nVar < 7 ~ "exclude")) %>%
  filter(exclude == "no") %>%
  mutate(modelGroup = paste0(response, "RandomNoHerbivores")) %>%
  mutate(randomEffect = ifelse(grepl("'re'", vars), TRUE, FALSE))

guide <- rbind(subGuideFullRandom, subGuideFullNonRandom, subGuideBest, subGuideRandomNoHerbivores, subGuideFullNonRandomNoHerbivores) 
unique(guide$modelGroup)


res <- data.table()
pred <- data.table()
pred.int <- data.table()

bm.spec.out <- data.table()
res.out <- data.table()
pred.out <- data.table()
pred.int.out <- data.table()

#modelGroup <- "BurnedAreaTrendFullRandom"

############### create cluster ####################
library(doSNOW)
library(foreach)
library(tictoc)

nCores <- parallel::detectCores()-4
# Create and register a cluster
clust <- makeCluster(nCores)
registerDoSNOW(clust)

rbindLists <- function(x, y) {combined.list <- list(res = rbind(x$res, y$res),
                                                    estimates = rbind(x$estimates, y$estimates),
                                                    pred = rbind(x$pred, y$pred), 
                                                    pred.int = rbind(x$pred.int, y$pred.int)) 
return(combined.list)}

#guide <- guide %>% sample_n(10)
##############################################################################            
################################## LOOOOOOOOOOOOP ############################            
##############################################################################    
tic()

for(modelGroup in unique(guide$modelGroup)){
  
  print(paste0("starting with: ", modelGroup))
  
  model.group <- modelGroup
  
  guide <- guide %>% data.table() %>% as_tibble() %>% as.data.table() 
  
  filter.resp <- unique(guide[modelGroup == model.group, ]$response)
  
  guide.sub <- guide %>% filter(modelGroup %in% c(model.group))
  
  #subset <- unique(guide[modelGroup == model.group, ]$subset)
  #dt.sub <- dt.mod %>% filter(eval(parse(text = subset)))
  
  #guide.sub <- guide.sub %>% sample_n(1)
  
  ## progress bar 
  iterations <- nrow(guide.sub)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  res <- foreach(i = 1:nrow(guide.sub),
                 .packages = c('mgcv', 'ggplot2', 'tidyr', 'data.table', 'tidyverse', 'MuMIn', 'broom'),
                 .options.snow = opts,
                 .inorder = FALSE,
                 .combine = rbind) %dopar% {
                   
                   
                   formula <- as.formula(guide.sub[i,]$formula)
                   
                   m <- tryCatch(
                     {bam(formula, data = dtMod, select = TRUE, method = "fREML")},
                     error = function(e) {cat("Model", i, "failed: ", e$message, "\n") 
                       return(NULL) })
                   
                   if(is.null(m)){next}
                   
                   m.sum <- summary(m)
                   
                   
                   
                   tmp <- data.frame(r_squared = NA, 
                                     aicc = NA, 
                                     aic = NA,
                                     bic = NA,
                                     formula = NA, 
                                     response = NA, 
                                     deviance_expl = NA)
                   
                   
                   tmp <- tmp %>% 
                     mutate(
                       r_squared = m.sum$r.sq, 
                       aicc = AICc(m), 
                       aic = AIC(m),
                       bic = BIC(m),
                       formula = guide.sub[i,]$formula, 
                       response = guide.sub[i,]$response, 
                       modelGroup = guide.sub[i,]$modelGroup,
                       deviance_expl = m.sum$dev.expl,
                       formulaID = guide.sub[i,]$formulaID, 
                       sphere = guide.sub[i,]$sphere,
                       interaction = guide.sub[i,]$interaction,
                       vars = guide.sub[i,]$vars, 
                       random_effect = guide.sub[i,]$random_effect
                     )
                   
                   #res <- rbind(res, tmp)
                   return(tmp)
                   
                   print(paste0(i, "/", nrow(guide.sub), " model group: ", modelGroup))
                 }
  
  
  bm.spec <- res %>% filter(modelGroup == model.group) %>%
    slice_min(aic) %>% 
    unique() 
  
  formula.bm <- as.formula(bm.spec$formula)
  
  
  m <- tryCatch(
    {gam(formula.bm, data = dtMod, select = TRUE, method = "REML")},
    error = function(e) {cat("Model", i, "failed: ", e$message, "\n") 
      return(NULL) })
  
  summary(m)
  bm.spec
  
  var.names <- bm.spec %>% dplyr::select(vars) %>% 
    mutate(sep_vars  = str_split(string = vars, pattern = "\\+")) %>% 
    dplyr::select(sep_vars) %>% 
    pull() %>% 
    map(~ .x[!grepl("by", .x)]) %>% 
    map(~ .x[!grepl("\\|", .x)]) %>% 
    unlist() %>% unique()
  
  vars.clean <- gsub("s\\(", "", var.names)
  vars.clean <- gsub("\\)", "", vars.clean)
  vars.clean <- gsub(", k = 3", "", vars.clean)
  vars.clean <- gsub("Biome, ", "", vars.clean)
  vars.clean <- gsub("\\, bs = 're'", "", vars.clean)
  
  var.names <- data.table(
    vars.clean = vars.clean) %>%
    filter(!grepl("Biome", vars.clean)) %>%
    pull() %>%
    trimws() %>% 
    unique()
  var.names
  
  
  var.names <- var.names %>% as.data.frame() %>% filter(!grepl("'re'", .)) %>% pull()
  
  var.names <- unique(var.names)
  
  if(length(var.names) == 1){if(var.names == ""){var.names <- "intercept"}}
  
  moderators <- bm.spec %>% dplyr::select(vars) %>% 
    mutate(sep_vars  = str_split(string = vars, pattern = "\\+")) %>% 
    dplyr::select(sep_vars) %>% 
    pull() %>%
    map(~ .x[grepl("by", .x)]) %>%
    unlist() %>% 
    as.data.table() %>% 
    mutate(moderators_raw  = str_split(string = ., pattern = "by")) %>%
    dplyr::select(moderators_raw) %>% pull() %>% 
    map(~ .x[!grepl("min_dist_enclosure_scaled", .x)]) %>% unlist()
  
  moderators <-  gsub("\\) ", "", moderators)
  moderators <-  gsub("\\)", "", moderators)
  
  moderators <-  gsub("\\=", "", moderators)
  moderators <-  gsub(" ", "", moderators)
  
  
  ### loop through vars and get pred
  
  if(length(var.names) > 0 & !any(grepl("intercept", var.names))){
    
    for(j in 1:length(var.names)){
      
      var <- var.names[j] 
      
      #filter.resp <- "EviTrend"
      
      marg.tmp <- partialPred(model = m, response = filter.resp,
                              var = var,
                              data = dtMod, newdata = dtMod %>% dplyr::select(-c(all_of(filter.resp)))) 
      
      marg.tmp <- marg.tmp %>% rename(varValue = paste0(gsub("_scaled", "", var))) %>% mutate(term = var,
                                                                                              cleanVar = gsub("log_", "", term),
                                                                                              cleanVar = gsub("_scaled", "", cleanVar), 
                                                                                              responseValue = dtMod %>% dplyr::select(c(all_of(filter.resp)))%>% pull(), 
                                                                                              Biome = dtMod$Biome)
      
      ggplot() + geom_line(aes(x = marg.tmp$varValue, y = marg.tmp$fit))
      
      if(j==1){
        marg <- marg.tmp}else{
          marg <- rbind(marg, marg.tmp)}
    }
    
    tidyNonPara <- tidy(m, parametric = FALSE) %>% 
      mutate(cleanVar = gsub("s\\(Biome,", "", term),   
             cleanVar = gsub("_scaled\\)", "", cleanVar),
             cleanVar = gsub("s\\(", "", cleanVar)) %>% 
      rename(pValue = p.value) %>% 
      dplyr::select(pValue, cleanVar)
    
    
    tidyPara <- tidy(m, parametric = TRUE) %>% 
      mutate(cleanVar = gsub("\\(", "", term),  
             cleanVar = gsub("_scaled\\)", "", cleanVar),
             cleanVar = gsub("\\)", "", cleanVar)) %>% 
      rename(pValue = p.value) %>% 
      dplyr::select(pValue, cleanVar)
    
    tidyM <- rbind(tidyPara, tidyNonPara)
    
    tmp.pred <- marg %>% 
      mutate(modelGroup = model.group, 
             response = filter.resp,
             formulaID = bm.spec$formulaID
      ) %>% left_join(tidyM) %>% 
      mutate(significance = ifelse(pValue < 0.05, "significant", "non significant"))
    
    pred <- rbind(tmp.pred, pred)
    
  }
  
  if(bm.spec$interaction == TRUE){
    
    for(k in 1:length(moderators)){
      
      #var.int <- "min_dist_enclosure_scaled"
      
      moderator = moderators[k]
      
      
      marg.tmp <- partialPred(model = m, response = filter.resp,
                              var = var.int, interaction = TRUE, moderator = moderator,
                              data = dtMod, newdata = dtMod %>% dplyr::select(-c(all_of(filter.resp))))
      
      moderator.clean <-  gsub("_scaled", "", moderator)
      
      marg.tmp.int <- marg.tmp %>% rename(varValue = paste0(gsub("_scaled", "", var.int))) %>% mutate(term = var.int,
                                                                                                      cleanVar = gsub("log_", "", term),
                                                                                                      cleanVar = gsub("_scaled", "", cleanVar),
                                                                                                      moderator = moderator.clean,
                                                                                                      Biome = dtMod$Biome,
                                                                                                      responseValue = dtMod %>% dplyr::select(c(all_of(filter.resp)))%>% pull(),
                                                                                                      moderatorValue = dtMod %>% dplyr::select(c(all_of(moderator.clean)))%>% pull())
      if(k==1){
        marg.int <- marg.tmp.int}else{
          marg.int <- rbind(marg.int, marg.tmp.int)}
    }
    
    
    tmp.pred.int <- marg.int %>% 
      mutate(modelGroup = model.group, 
             response = filter.resp,
             formulaID = bm.spec$formulaID)
    
    pred.int <- rbind(tmp.pred.int, pred.int)
    
  }
  
  bm.spec.out <- rbind(bm.spec, bm.spec.out)
  res.out <- rbind(res, res.out)
  pred.out <- rbind(pred, pred.out)
  pred.int.out <- rbind(pred.int, pred.int.out)
  
  print(paste0(modelGroup, " done"))
  
  
}


print("loop done")

toc()
stopCluster(clust)

modelResults <- list(bm.spec = bm.spec.out, res = res.out, pred = pred.out, pred.int = pred.int.out)
saveRDS(modelResults, "builds/modelOutputs/exploratoryGamReserveRes.Rds")

foreach.results
stop(pb)







