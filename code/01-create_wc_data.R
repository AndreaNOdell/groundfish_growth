library(nwfscSurvey)
library(sdmTMB)
library(dplyr)
library(stringr)




# all catch data
  # catch <- nwfscSurvey::PullCatch.fn(SurveyName = "NWFSC.Combo", sample_types = c(NA, "NA", #"Life Stage", "Size")[1:2])
  # catch_hake <- nwfscSurvey::PullCatch.fn(SurveyName = "NWFSC.Combo", Name = "Pacific hake", #sample_types = c(NA, "NA", "Life Stage", "Size")) # hake catch data - more complete info
  # updated_catch = rbind(catch[!catch$Common_name == "Pacific hake",], catch_hake) # merge all catch data without hake, and then the complete hake catch data.
  # saveRDS(updated_catch, "survey_data/wcbts_catch_2024-01-09.rds")

# Haul 
  # haul <- nwfscSurvey::PullHaul.fn(SurveyName = "NWFSC.Combo")
  # saveRDS(haul, "survey_data/wcbts_haul_2023-11-28.rds")

# Bio
  # bio <- nwfscSurvey::PullBio.fn(SurveyName = "NWFSC.Combo")
  # saveRDS(bio, "survey_data/wcbts_bio_2023-11-28.rds")

# spec_names = group_by(DataPull, common_name) %>%
  # dplyr::summarize(scientific_name = scientific_name[1])
# spec_names$scientific_name = tolower(spec_names$scientific_name)
# saveRDS(spec_names,"nwfsc_lookup.rds")

spec_names <- readRDS("survey_data/nwfsc_species_lookup.rds")

bio <- readRDS("survey_data/wcbts_bio_2023-11-28.rds")
names(bio) <- tolower(names(bio))

haul <- readRDS("survey_data/wcbts_haul_2023-11-28.rds")
haul$year <- as.numeric(substr(haul$date_yyyymmdd, 1, 4))
haul$month <- as.numeric(substr(haul$date_yyyymmdd, 5, 6))
haul$day <- as.numeric(substr(haul$date_yyyymmdd, 7, 8))
haul$yday <- lubridate::parse_date_time(haul$date_yyyymmdd, orders="ymd")

haul <- dplyr::rename(haul,
                      o2 = o2_at_gear_ml_per_l_der,
                      degc = temperature_at_gear_c_der,
                      depthm = depth_hi_prec_m
) %>%
  dplyr::select(
    o2, degc, depthm, 
    performance, trawl_id, yday
  )
haul$trawl_id <- as.character(haul$trawl_id)
  
dat <- dplyr::left_join(bio, haul)

# get rid of fish without ages
dat <- dplyr::filter(dat, !is.na(age))
# get rid of fish with unknown sex
dat <- dplyr::filter(dat, sex !="U")

# filter to only include species with > 500 samples historically
thresh <- 500
dat <- dplyr::group_by(dat, scientific_name) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(n > thresh) %>%
  dplyr::select(-n)

# join in spec names
dat$scientific_name <- tolower(dat$scientific_name)
dat <- dplyr::left_join(dat, spec_names)

saveRDS(dat, "survey_data/data_for_length_models.rds")

