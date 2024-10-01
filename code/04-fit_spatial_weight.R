library(sdmTMB)
library(nwfscSurvey)

# Load in WCGBTS data (from 01_create_wc_data.R)
dat <- readRDS("survey_data/data_for_length_models.rds")

# Remove unspecific or missing species information
dat = dat %>% 
  filter(!scientific_name %in% c("sebastes sp. (aleutianus / melanostictus)", 
                                 "sebastes sp. (miniatus / crocotulus)" )) %>%  
  tidyr::separate_wider_delim(scientific_name, " ", names = c("Genus", "Species")) %>% 
  filter(!is.na(Genus))

# Combine WCGBTS and Fishbase datasets
dat$common_name = tolower(dat$common_name) 
dat$common_name <- sub(" ", "_", dat$common_name) # replace spaces with underscores in common name

spec <- readRDS("survey_data/spec_list.rds")

# run model for each species
for(i in 1:length(spec)) {
  
  sub <- dplyr::filter(dat, common_name == spec[i], 
                       !is.na(weight), !is.na(length_cm),
                       sex == "F") %>% 
    mutate(ln_weight_g = log(weight), ln_length_cm = log(length_cm)) # log transform data
  
  # Create the mesh
  sub <- add_utm_columns(sub, ll_names = c("longitude_dd", "latitude_dd"))
  mesh <- make_mesh(sub, xy_cols = c("X","Y"), cutoff = 50)
  
  # Center and scale data
  year_pars <- c(mean(sub$year), sd(sub$year))
  sub$year_scaled <- (sub$year - year_pars[1])/year_pars[2]
  sub$fyear <- as.factor(sub$year)
  sub$fyear_scaled <- as.factor(sub$year_scaled)
  sub$length_scaled <- scale(sub$length_cm)
  sub$doy <- lubridate::yday(sub$yday) 
  doy_pars <- c(mean(sub$doy), sd(sub$doy))
  sub$doy_scaled <- (sub$doy - doy_pars[1])/doy_pars[2]
  sub$doy_scaled2 <- sub$doy_scaled^2

# estimated weight for Le Cren's calculation
null_model <- sdmTMB(
    formula = ln_weight_g ~ ln_length_cm,
    data = sub,
    time = NULL,
    mesh = mesh,
    family = student(link = "identity", df = 5),
    spatiotemporal = "off",
    spatial = "off",
    silent = TRUE,
    reml = FALSE
  ) 

# extract a & b parameters  
a <- as.numeric(tidy(null_model)[1, "estimate"])
b <- as.numeric(tidy(null_model)[2, "estimate"])

# predicted weight for each length
sub$weight_g_pred <- exp(a)*sub$length_cm^b
  
# visual checks
plot(sub$weight ~ sub$length_cm)
plot(sub$weight_g_pred ~ sub$length_cm)
plot(sub$weight_g_pred ~ sub$weight); abline(a = 0, b = 1, col = "red")

sub$le_cren <- sub$weight / sub$weight_g_pred
sub$log_le_cren <- log(sub$le_cren)
hist(sub$le_cren)
 
# fit spatiotemporal model 
mfull <- sdmTMB(
  formula = log_le_cren ~ fyear_scaled + doy_scaled + doy_scaled2, 
  data = sub,
  time = "year",
  mesh = mesh, 
  family = student(link = "identity", df = 5),
  spatiotemporal = "iid",
  spatial = "on",
  silent = TRUE,
  reml = FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

# check model fit
sanity(mfull)

# Make spatial grid
res <- 4 # in km
spat.pred = add_utm_columns(availablecells[, c("Cent.Lat", "Cent.Long")], ll_names = c("Cent.Long", "Cent.Lat"))  %>% 
  mutate(X = floor(X/res)*res, Y = floor(Y/res)*res) %>% 
  distinct(X, Y)

# make prediction grid
pred <- expand.grid(year = unique(sub$year))
pred$fyear_scaled = as.factor((pred$year - year_pars[1])/year_pars[2])
pred$doy = 182 #july 1st
pred$doy_scaled = sub[sub$doy == 182, ]$doy_scaled[1]
pred$doy_scaled2 = sub[sub$doy == 182, ]$doy_scaled2[1]

# merge spatial and prediction grid
pred = cross_join(pred, spat.pred)

# make predictions
pred_lecren <- predict(mfull, pred) %>% 
  mutate(common_name = sub$common_name[1])

saveRDS(mfull, paste0("fitted_models/lecren/",sub$common_name[1],"_pred_lecren.rds"))
saveRDS(pred_lecren, paste0("output/lecren/",sub$common_name[1],"_pred_lecren.rds"))
}
