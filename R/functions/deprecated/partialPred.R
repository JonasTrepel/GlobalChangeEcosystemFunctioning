## function for partial dependence plots 

## will only work if you have a dataframe with unscaled and untransformed variables in the same row order than the predicted outcome
## also naming scheme must be: scaled variables contain "_scaled" and log-transformed variables contain "log_"


#### IMPORTANT: All variables used in the model must have the ending "_scaled" - even it they are not. Otherwise the code won#t take the mean/Mode


partialPred <- function(newdata, var, data, response, model,
                        interaction = FALSE, moderator){
  
  ## load packages 
  
  library(tidyverse)
  library(mgcv)
  library(data.table)
  library(MetBrewer)
  
  Mode <- function(x, na.rm = FALSE) {
    if(na.rm){
      x = x[!is.na(x)]
    }
    
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }
  
  ### interaction n ------------------------
  if(!interaction == TRUE){  
    
    for(name in names(newdata)){
      
      #don't take the mean if it's our variable of interest
      if(name == var){next}

      
      #only take the mean if the name doesn't contain scaled
      if(!grepl("scaled", name)){next}
      
      x_col <- newdata %>% 
        dplyr::select({{name}}) %>% 
        as.data.table()
      
      vect <- c(x_col[,1])
      
      vect <- c(vect[[1]])
      
      if(is.list(vect)){
        
        vect <- vect[[1]]
      }
      
      if(sum(is.na(vect)) == nrow(data)){next}
      
      if(is.character(vect[[1]])){
        
        tmp <- newdata %>% 
          mutate(mean_col = Mode(vect, na.rm = TRUE)) %>%
          as.data.table() %>% 
          dplyr::select(mean_col)
        
      }else{
        
        tmp <- newdata %>% 
          mutate(mean_col = mean(vect, na.rm = TRUE)) %>%
          as.data.table() %>% 
          dplyr::select(mean_col)
      }
      
      var_name <- paste0(name)
      
      setnames(tmp, "mean_col", (var_name))
      
      newdata <- newdata %>% 
        dplyr::select(-c({{var_name}})) %>% 
        cbind(tmp)
      
    }
    
    ## now new data should contain the mean of all columns, expect of the variable of interest 
    
    ## predict -----------
    pred <- as.data.frame(predict(model, newdata, se.fit = TRUE))
    
    newdata.pred <- newdata %>% 
      cbind(pred) %>% 
      mutate(
        ci.ub = fit + 1.96*se.fit, 
        ci.lb = fit - 1.96*se.fit)
    
    
    ## get clean variable name 
    clean.var.scaled <- gsub("log_", "", var)
    clean.var <- gsub("_scaled", "", clean.var.scaled)
    dt.pred <- newdata.pred %>% dplyr::select(all_of(clean.var), fit, ci.ub, ci.lb)
    
    
    ## interaction = T -------------
  }else{
    
    newdata.high <-  newdata
    
    for(name in names(newdata.high)){
      
      #don't take the mean if it's our variable of interest
      if(name == var){next}
      
      #only take the mean if the name contains scaled
      if(!grepl("scaled", name)){next}

      if(name == moderator){
        
        x_col <- newdata.high %>% 
          dplyr::select({{name}}) %>% 
          as.data.table()
        
        vect <- c(x_col[,1])
        
        if(is.character(vect[[1]])){next}
        
        vect <- c(vect[[1]])
        
        if(is.list(vect)){
          
          vect <- vect[[1]]
        }
        
        if(sum(is.na(vect)) == nrow(data)){next}
        
        tmp.high <- newdata.high %>% 
          mutate(mean_col = mean(vect, na.rm = TRUE) + sd(vect, na.rm = TRUE)) %>%
          as.data.table() %>% 
          dplyr::select(mean_col)
        
      }else{
        
        x_col <- newdata.high %>% 
          dplyr::select({{name}}) %>% 
          as.data.table()
        
        vect <- c(x_col[,1])
        
        vect <- c(vect[[1]])
        
        if(is.list(vect)){
          
          vect <- vect[[1]]
        }
        
        if(sum(is.na(vect)) == nrow(data)){next}
        
        if(is.character(vect[[1]])){
          
          tmp.high <- newdata.high %>% 
            mutate(mean_col = Mode(vect, na.rm = TRUE)) %>%
            as.data.table() %>% 
            dplyr::select(mean_col)
          
        }else{
          
          tmp.high <- newdata.high %>% 
            mutate(mean_col = mean(vect, na.rm = TRUE)) %>%
            as.data.table() %>% 
            dplyr::select(mean_col)
        }
        
      }
      
      
      var_name <- paste0(name)
      
      setnames(tmp.high, "mean_col", (var_name))
      
      newdata.high <- newdata.high %>% 
        dplyr::select(-c({{var_name}})) %>% 
        cbind(tmp.high)
      
    } ### end for loop high names 
    
    ## now new data should contain the mean of all columns, expect of the variable of interest 
    
    ## predict for high -----------
    pred.high <- as.data.frame(predict(model, newdata.high, se.fit = TRUE))
    
    
    
    newdata.high.pred <- newdata.high %>% 
      cbind(pred.high) %>% 
      mutate(
        ci.ub = fit + 1.96*se.fit, 
        ci.lb = fit - 1.96*se.fit,
        ci.ub.high = fit + 1.96*se.fit, 
        ci.lb.high = fit - 1.96*se.fit,
        fit.high = fit,
        linetype = "high")
    
    
    ### low level of moderatpr 
    newdata.low <-  newdata
    
    
    for(name in names(newdata.low)){
      
      #don't take the mean if it's our variable of interest
      if(name == var){next}
      
      #only take the mean if the name contains scaled
      if(!grepl("scaled", name)){next}
      
      
      if(name == moderator){
        
        x_col <- newdata.low %>% 
          dplyr::select({{name}}) %>% 
          as.data.table()
        
        vect <- c(x_col[,1])
        
        if(is.character(vect[[1]])){next}
        
        vect <- c(vect[[1]])
        
        if(is.list(vect)){
          
          vect <- vect[[1]]
        }
        
        if(sum(is.na(vect)) == nrow(data)){next}
        
        tmp.low <- newdata.low %>% 
          mutate(mean_col = mean(vect, na.rm = TRUE) - sd(vect, na.rm = TRUE)) %>%
          as.data.table() %>% 
          dplyr::select(mean_col)
        
      }else{
        
        x_col <- newdata.low %>% 
          dplyr::select({{name}}) %>% 
          as.data.table()
        
        vect <- c(x_col[,1])
        
        vect <- c(vect[[1]])
        
        if(is.list(vect)){
          
          vect <- vect[[1]]
        }
        
        if(sum(is.na(vect)) == nrow(data)){next}
        
        if(is.character(vect[[1]])){
          
          tmp.low <- newdata.low %>% 
            mutate(mean_col = Mode(vect, na.rm = TRUE)) %>%
            as.data.table() %>% 
            dplyr::select(mean_col)
          
        }else{
          
          tmp.low <- newdata.low %>% 
            mutate(mean_col = mean(vect, na.rm = TRUE)) %>%
            as.data.table() %>% 
            dplyr::select(mean_col)
        }
        
      }
      
      
      var_name <- paste0(name)
      
      setnames(tmp.low, "mean_col", (var_name))
      
      newdata.low <- newdata.low %>% 
        dplyr::select(-c({{var_name}})) %>% 
        cbind(tmp.low)
      
    } ### end for loop low names 
    
    ## now new data should contain the mean of all columns, expect of the variable of interest and the moderator 
    
    ####predict for low ------------------
    pred.low <- as.data.frame(predict(model, newdata.low, se.fit = TRUE))
    
    
    
    newdata.low.pred <- newdata.low %>% 
      cbind(pred.low) %>% 
      mutate(
        ci.ub = fit + 1.96*se.fit, 
        ci.lb = fit - 1.96*se.fit,
        ci.ub.low = fit + 1.96*se.fit, 
        ci.lb.low = fit - 1.96*se.fit,
        fit.low = fit,
        linetype = "low")
    
    ## get clean variable name 
    clean.var.scaled <- gsub("log_", "", var)
    clean.var <- gsub("_scaled", "", clean.var.scaled)
    
    dt.pred.low <- newdata.low.pred %>% dplyr::select(all_of(clean.var), fit.low, ci.ub.low, ci.lb.low)
    dt.pred.high <- newdata.high.pred %>% dplyr::select(fit.high, ci.ub.high, ci.lb.high)
    dt.pred <- cbind(dt.pred.low, dt.pred.high)
    print("return interaction") 
    
    
  } 
  
  return(dt.pred)
  
}