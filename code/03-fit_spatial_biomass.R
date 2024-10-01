library(sdmTMB)
library(tidyverse)
library(sp)

# Load in WCGBTS catch data
dat <- readRDS("survey_data/wcbts_catch_2024-01-09.rds")
spec <- readRDS("survey_data/spec_list.rds")

# Clean and organize data
names(dat) <- tolower(names(dat))
dat$common_name = tolower(dat$common_name) 
dat$common_name <- sub(" ", "_", dat$common_name) # replace spaces with underscores in common name

# center and scale data
dat = dat[dat$common_name %in% spec,] %>% 
  mutate(fyear = as.factor(year))
dat$doy <- lubridate::yday(dat$date)
doy_pars <- c(mean(dat$doy), sd(dat$doy))
dat$doy_scaled <- (dat$doy - doy_pars[1])/doy_pars[2]
dat$doy_scaled2 <- dat$doy_scaled^2

# run model for each species
for(i in 1:length(spec)) {

  # filter data for specific species
  sub <- dplyr::filter(dat, 
                       common_name == spec[i],
                       !is.na(cpue_kg_km2)) # remove CPUE NA's if any
  
  # Create the mesh
  sub = add_utm_columns(sub, ll_names = c("longitude_dd", "latitude_dd"))
  mesh <- make_mesh(sub, xy_cols = c("X","Y"), cutoff = 30)
  
  # fit biomass model
  fit_CPUE_st <- sdmTMB(cpue_kg_km2 ~ fyear + doy_scaled + doy_scaled2, 
                        mesh = mesh,
                        spatiotemporal = "iid", 
                        spatial = "on",
                        data = sub,
                        time = "year",
                        family = tweedie(link = "log"),
                        #extra_time = missing_years,
                        control = sdmTMBcontrol(newton_loops = 1))
  
  sanity(fit_CPUE_st)

  # create spatial grid
  res <- 4
  spat.pred = add_utm_columns(availablecells[, c("Cent.Lat", "Cent.Long")], ll_names = c("Cent.Long", "Cent.Lat")) %>% 
    mutate(X = floor(X/res)*res, Y = floor(Y/res)*res) %>% 
    distinct(X, Y)
  
  # create prediction grid
  pred = as.data.frame(cbind(year = unique(sub$year)))
  pred$fyear = as.factor(pred$year)
  pred$doy = 182 # July 1st
  pred$doy_scaled = sub[sub$doy == 182, ]$doy_scaled[1]
  pred$doy_scaled2 = sub[sub$doy == 182, ]$doy_scaled2[1]
  
  # merge spatial and prediction grid
  pred = cross_join(pred, spat.pred)
  
  # make predictions
  pred_cpue_spatial = predict(fit_CPUE_st, pred, return_tmb_object = FALSE)
  pred_cpue <- predict(fit_CPUE_st, pred, return_tmb_object = TRUE)
  index <- get_index(pred_cpue, bias_correct = TRUE)
  
  saveRDS(fit_CPUE_st, paste0("fitted_models/biomass/",sub$common_name[1],"_pred_cpue.rds"))
  saveRDS(index, paste0("output/biomass/",sub$common_name[1],"_pred_index.rds"))
  saveRDS(pred_cpue_spatial, paste0("output/biomass/spatial_preds/",sub$common_name[1],"_pred_index.rds"))
}





