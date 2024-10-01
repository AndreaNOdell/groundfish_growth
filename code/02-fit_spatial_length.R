library(sdmTMB)
library(rfishbase)
library(tidyverse)
library(nwfscSurvey)

# Load in WCGBTS data (from 01_create_wc_data.R)
dat <- readRDS("survey_data/data_for_length_models.rds")

# Load in Fishbase dataset with species and Lmax information
sp.index = fb_tbl("species")[,c("SpecCode", "Genus", "Species", "Length", "LTypeMaxM", "LengthFemale", "LTypeMaxF")] %>% 
  rename("Lmax" = "Length", "LmaxF" = "LengthFemale") 
sp.index$Genus = tolower(sp.index$Genus)

# Remove unspecific or missing species information
dat = dat %>% 
  filter(!scientific_name %in% c("sebastes sp. (aleutianus / melanostictus)", 
                                 "sebastes sp. (miniatus / crocotulus)" )) %>%  
  tidyr::separate_wider_delim(scientific_name, " ", names = c("Genus", "Species")) %>% 
  filter(!is.na(Genus))

# Combine WCGBTS and Fishbase datasets
dat = left_join(dat, sp.index)
dat$common_name = tolower(dat$common_name) 
dat$common_name <- sub(" ", "_", dat$common_name) # replace spaces with underscores in common name

spec <- readRDS("survey_data/spec_list.rds")

# Generate table of Lmax values for species
Lmax_tbl = dat %>% 
  filter(common_name %in% spec) %>% 
  distinct(common_name, Lmax, LTypeMaxM, LmaxF, LTypeMaxF)
# saveRDS(Lmax_tbl, paste0("output/lmax_tbl.rds"))

for(i in 2:length(spec)) {
  
  # Filter dataset to specific species and female only. Remove NAs
  sub <- dplyr::filter(dat, common_name == spec[i], 
                       !is.na(age), !is.na(length_cm),
                       sex == "F")
  
  # Calculate linearized response y_lin (use female Lmax if available) -- data extracted from fishbase 
  sub = sub %>%
    mutate(y_lin = if_else(!is.na(sub$LmaxF), 
                           log(sub$LmaxF - sub$length_cm), # female-specific Lmax
                           log(sub$Lmax - sub$length_cm))) # general Lmax (if no female estimate)
  
  # Create the mesh
  sub <- add_utm_columns(sub, ll_names = c("longitude_dd", "latitude_dd"))
  mesh <- make_mesh(sub, xy_cols = c("X","Y"), cutoff = 30)
  
  # center and scale data
  year_pars <- c(mean(sub$year), sd(sub$year))
  sub$year_scaled <- (sub$year - year_pars[1])/year_pars[2]
  age_pars <- c(mean(sub$age), sd(sub$age))
  sub$age_scaled <- (sub$age - age_pars[1])/age_pars[2]
  sub$fyear <- as.factor(sub$year)
  sub$doy <- lubridate::yday(sub$yday) 
  doy_pars <- c(mean(sub$doy), sd(sub$doy))
  sub$doy_scaled <- (sub$doy - doy_pars[1])/doy_pars[2]
  sub$doy_scaled2 <- sub$doy_scaled^2
   
                           
# Model
map_list <- list(ln_tau_Z = factor(rep(1, times = length(unique(sub$fyear))))) 
fit_s_agexyear_st_2 <- sdmTMB(y_lin ~  age + fyear + doy_scaled + doy_scaled2, 
                             spatial_varying = ~ age:fyear,
                             mesh = mesh,
                             spatiotemporal = "off", # off because random fields in age interaction
                             data = sub,
                             time="year",
                             control = sdmTMBcontrol(newton_loops = 1,
                                                     map = map_list))
# check fit
sanity(fit_s_agexyear_st_2) 
  
  # create spatial grid
  res <- 4   # 4x4km grid
  spat.pred = add_utm_columns(availablecells[, c("Cent.Lat", "Cent.Long")], ll_names = c("Cent.Long", "Cent.Lat"))  %>% 
    mutate(X = floor(X/res)*res, Y = floor(Y/res)*res) %>% 
    distinct(X, Y)
  
  # create prediction grid
  pred <- expand.grid(age = unique(sub$age),
                      year = unique(sub$year),
                      sex = unique(sub$sex))
  pred$fyear = as.factor(pred$year)
  pred$doy = 182 #july 1st
  pred$doy_scaled = sub[sub$doy == 182, ]$doy_scaled[1]
  pred$doy_scaled2 = sub[sub$doy == 182, ]$doy_scaled2[1]
  
  # merge spatial and predition grid
  pred = cross_join(pred, spat.pred)
  
  # make predictions
  pred_len <- predict(fit_s_agexyear_st_2, pred)

  # save data
  saveRDS(fit_s_agexyear_st_2, paste0("fitted_models/growth/",sub$common_name[1],"_pred_length.rds"))
  saveRDS(pred_len, paste0("output/growth/",sub$common_name[1],"_pred_length.rds"))
}



