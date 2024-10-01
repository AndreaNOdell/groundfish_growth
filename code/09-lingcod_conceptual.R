library(sdmTMB)
library(rfishbase)
library(tidyverse)
library(nwfscSurvey)
library(glmmTMB)
library(viridis)
library(tidyr)
library(ggpubr)

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

for(i in 4:4) {
  
  # Filter dataset to specific species and female only. Remove NAs
  sub <- dplyr::filter(dat, common_name == spec[i], 
                       !is.na(age), !is.na(length_cm),
                       sex == "F")
  
  # Calculate linearized response y_lin (use female Lmax if available) -- data extracted from fishbase 
  #sub$LmaxF = sub$Lmax = 93.4
  #sub$Lmax = max(sub$length_cm + 1)
  sub = sub %>%
    mutate(y_lin = if_else(!is.na(sub$LmaxF), 
                           log(sub$LmaxF - sub$length_cm), 
                           log(sub$Lmax - sub$length_cm)))
  
  # Create the mesh
  #sub <- add_utm_columns(sub, ll_names = c("longitude_dd", "latitude_dd"))
  #mesh <- make_mesh(sub, xy_cols = c("X","Y"), cutoff = 30)
  
  # set up data
  sub$fyear <- as.factor(sub$year)

  fit <- glmmTMB(y_lin ~ -1 + (age|fyear), data = sub)
  coefs <- ranef(fit)[[1]][[1]]               
  
  dat <- data.frame(
    year = c(2003:2019),
    intercept = coefs[,1],
    slope = coefs[,2]
  )
  
  # Create a range of x values
  x_values <- seq(0,20,by=0.1)
  
  # Function to calculate predictions for a given year
  predict_for_year <- function(intercept, slope, x_values) {
    sapply(x_values, function(x) intercept + slope * x)
  }
  
  # Calculate predictions for each year
  predictions <- sapply(1:nrow(dat), function(i) {
    predict_for_year(dat$intercept[i], dat$slope[i], x_values)
  })
  predictions_df <- as.data.frame(predictions)
  predictions_long <- predictions_df %>%
    mutate(x_value = x_values) %>% # Add x values as a column
    gather(key = "year", value = "prediction", -x_value) # Gather year columns into long format
  
  yeardf <- data.frame(year = paste0("V",seq(1,17)), years = dat$year)
  # Convert the 'year' column to numeric, if it's stored as a character
  predictions_long <- dplyr::left_join(predictions_long, yeardf)
}

# Use recent pub on lingcod
#https://www.researchgate.net/profile/Gary-Longo/publication/351631903_Geographic_variability_in_lingcod_Ophiodon_elongatus_life-history_and_demography_along_the_US_West_Coast_oceanographic_drivers_and_management_implications/links/618d887b07be5f31b76e8ff7/Geographic-variability-in-lingcod-Ophiodon-elongatus-life-history-and-demography-along-the-US-West-Coast-oceanographic-drivers-and-management-implications.pdf

Linf <- 93.4 # coastwide estimates
k <- 0.266
t0 <- 0

df = data.frame(ages = seq(0,20,by=0.1))
df$pred_len <- Linf * (1 - exp(-k*(df$ages - t0)))
df$growth_inc <- c(NA, diff(log(df$pred_len)))
df$logLinf <- log(Linf - df$pred_len)

p1 <- ggplot(df, aes(ages, pred_len)) + 
  geom_line(col = viridis(1), linewidth=1.2) + 
  xlab("Age (years)") + ylab("Length (cm)") + theme_bw()
p2 <- ggplot(df, aes(ages, logLinf)) + 
  geom_line(col = viridis(1), linewidth=1.2) + xlab("Age (years)") + 
  ylab(expression(L[infinity] ~ "- length (cm)")) + theme_bw()

predictions_long$Year <- predictions_long$years
p3 <- ggplot(predictions_long, aes(x_value, prediction, col = Year, group=Year)) + 
  geom_line() + scale_color_viridis(option="magma") + theme_bw() + 
  xlab("Age (years)") + ylab(expression(L[infinity] ~ "- length (cm)"))

coefs$year <- as.numeric(rownames(coefs))
p4 <- ggplot(coefs, aes(year, -age)) + 
  geom_line(col = viridis(1), linewidth=1.2, alpha=0.5) + 
  geom_point(col = viridis(1),size=3) + 
  scale_color_viridis() + 
  theme_bw() + xlab("Year") + ylab("Estimated k")

final_p = ggarrange(p1,p2,p3,p4,
                    nrow=2, ncol = 2,
                    labels = c("a)", "b)", "c)", "d)"))

ggsave(final_p, file=paste0("plots/final_figures/fig1.jpg"), width =170, height = 106, dpi = 400, units = "mm")




