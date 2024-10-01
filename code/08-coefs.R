library(tidyverse)
library(sdmTMB)


spec <- readRDS("survey_data/spec_list.rds")
coef_list = list()


for(i in 1:length(spec)){
x = readRDS(paste0("fitted_models/lecren/",spec[i],"_pred_lecren.rds")) 

coef_df = tidy(x, effects = "ran_pars", conf.int = TRUE)[c(3,4),1:2] %>% 
  mutate(common_name = spec[i])

coef_list[[i]] = coef_df

}


coef_df = bind_rows(coef_list) 


condition_sigma = coef_df %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  mutate(fig_names = sub("_", " ", common_name)) 
#saveRDS(condition_sigma, paste0("output/index/condition_sigma.rds"))
