extract_gls_estimates <- function(mod, part = TRUE, resp = NA){
  
  library(tidyverse)
  
  if(part == FALSE){
    dt_est <- data.table(
      term = names(mod$coefficients), 
      p_value = mod$pval_t,
      estimate = mod$coefficients, 
      std_error = mod$SE, 
      response = str_split_i(mod$formula, " ", 1)[1]
    ) %>% mutate(ci_lb = estimate - 1.96*std_error,
                 ci_ub = estimate + 1.96*std_error, 
                 sig = ifelse(p_value < 0.05, "significant", "non-significant"))
  }else if(part == TRUE){
    
    dt_est <- as.data.frame(mod$overall$t.test) %>% 
      rownames_to_column(var = "term") %>% 
      rename(estimate = Est, 
             std_error = SE, 
             p_value = pval.t) %>% mutate(ci_lb = estimate - 1.96*std_error,
                                          ci_ub = estimate + 1.96*std_error, 
                                          sig = ifelse(p_value < 0.05, "significant", "non-significant"), 
                                          response = resp)
    
  }
  
  return(dt_est)
}