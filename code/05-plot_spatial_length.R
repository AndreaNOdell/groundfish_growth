library(ggplot2)
library(tidyverse)
library(sf)
library(broom)

spec <- readRDS("survey_data/spec_list.rds")
ref_year = 2007

# generate map
map_data <- rnaturalearth::ne_states(
  returnclass = "sf", country = "united states of america")

# Crop the polygon for plotting and efficiency:
usa <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = -145, ymin = 25, xmax = -117, ymax = 48.5))))

coast <- sf::st_transform(usa, crs = 3157)

# create empty list to save data
dat_list = list()


for(i in 1:length(spec)) {
  # load model oututs
  d <- readRDS(paste0("output/growth/",spec[i],"_pred_length.rds"))
  age_coef <- readRDS(paste0("fitted_models/growth/",spec[i],"_pred_length.rds")) %>% 
    tidy() %>% 
    filter(term == "age") %>% 
    select(estimate) %>% 
    as.numeric()
  
  # organize dataframe
  d <- dplyr::filter(d, age == min(d$age)) %>% 
    pivot_longer(cols = starts_with("zeta_s_age:fyear"),
                 names_to = "zeta_s_fyear",
                 names_prefix = "zeta_s_age:fyear",
                 values_to = "zeta_s_vals")  %>% 
    mutate(zeta_vals_age_adjusted = zeta_s_vals + age_coef) %>% 
    mutate(k_f = -(zeta_vals_age_adjusted)) %>% # calculate k 
    filter(year == ref_year) %>% 
    # global_mean = the average k_f across all years and all spatial coordinates
    mutate(global_mean = mean(k_f)) %>% 
    group_by(X, Y) %>% 
    # mean_k_f = for each spatial coordinate, take the average k_f across all years
    # anom_k_f = for each spatial coordinate, it is the difference between k_f for that year and the average k_f across all years (mean_k_f) at that specific spatial coordinate
    mutate(mean_k_f = mean(k_f), sd_k_f = sd(k_f)) %>% 
    mutate(anom_k_f = k_f - mean_k_f) %>% 
    mutate(common_name = spec[i])
 
  # save data frame
  dat_list[[i]] = d
  
  # average k_f spatial map
  p1 <- ggplot(coast) + geom_sf() + 
    geom_tile(data = d, aes(x = floor(X*1000), y = floor(Y*1000), fill = mean_k_f)) + 
    scale_fill_viridis_c(name=expression(k[f])) + 
    theme_bw() + 
    xlab("") + ylab("") +
    labs(subtitle = paste0(spec[i])) +
    theme(strip.background = element_rect(fill="white"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # annual k_f anomaly spatial map
  p2 <- ggplot(coast) + geom_sf() + 
    geom_tile(data = d, aes(x = X*1000, y = Y*1000, fill = anom_k_f)) +
    scale_fill_gradient2(name=expression(k[f] ~ "anom")) + 
    theme_bw() + 
    facet_wrap(~ zeta_s_fyear) +
    xlab("") + ylab("") + 
    labs(subtitle = paste0(spec[i])) +
    theme(strip.background = element_rect(fill="white"),
          axis.text.x = element_text(angle = 45, hjust = 1))
    

# make plots
 #ggsave(p1, file=paste0("plots/map/average/k_st_avg_length_", spec[i], ".png"))
 #ggsave(p2, file=paste0("plots/map/annual/k_st_annual_length_", spec[i], ".png"), width =7, height = 7)
}

# plot variability in k_f --------
pal = c("#69b27b","#c86dbb","#7a87d4","#76bc50","#d99431","#ac57c5", "#51bdc6","#d85e73","#b2a84f","#ca764f")

dat_df = bind_rows(dat_list) %>% # bind k_f predictions for each species into one df
  group_by(common_name) %>% 
  mutate(norm_k_f = ((mean_k_f - min(mean_k_f))/(max(mean_k_f)- min(mean_k_f))),
         norm_sd = ((sd_k_f - min(sd_k_f))/(max(sd_k_f)- min(sd_k_f)))) 
dat_df$fig_names <- str_to_sentence(sub("_", " ", dat_df$common_name)) # replace spaces with underscores in common name

#saveRDS(dat_df, paste0("output/index/dat_df.rds"))

dat_df = dat_df %>%
  mutate(fig_names = factor(fig_names,
                            levels = c("Dover sole","Pacific sanddab",  "Petrale sole", "Aurora rockfish","Darkblotched rockfish", "Splitnose rockfish", "Lingcod", "Pacific hake",  "Sablefish")))

p1 = ggplot(coast) + geom_sf() + 
  geom_tile(data = dat_df, aes(x = X*1000, y = Y*1000, fill = norm_k_f)) +
  scale_fill_viridis_c(name=expression("Normalized mean" ~ k[s])) + 
  theme_bw() + 
  xlab("") + ylab("") +
  facet_wrap(~fig_names, ncol = 3) +
  theme(strip.background = element_rect(fill="lightgray"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8.5))


p2 = dat_df %>% 
  group_by(common_name, zeta_s_fyear) %>% # calculate mean for each year (temporal variability)
  mutate(year_k_f = mean(k_f)) %>% 
  ungroup() %>% 
  group_by(fig_names) %>% 
  summarise(spat_cv = sd(mean_k_f)/mean(mean_k_f), year_cv = sd(year_k_f)/mean(year_k_f)) %>% 
  ggplot(aes(x = spat_cv, y = year_cv)) +
  geom_point(aes(fill = fig_names), shape = 21, col = "black", size = 5) +
  scale_fill_manual(values= pal) +
  theme_classic() +
  lims(x = c(0, 0.83), y = c(0, 0.25)) +
  labs(x = "spatial CV in k", y = "temporal CV in k", fill = "Species")

p5 = ggplot(coast) + geom_sf() + 
  geom_tile(data = dat_df, aes(x = X*1000, y = Y*1000, fill = mean_k_f)) +
  scale_fill_viridis_c(name=expression(k[f])) + 
  theme_bw() + 
  xlab("") + ylab("") +
  facet_wrap(~common_name) +
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 5))
  

ggsave(p1, file=paste0("plots/final_figures/fig3.jpg"), width = 170, height = 255, dpi = 500, units = "mm")
ggsave(p2, file=paste0("plots/spat_temp_CV.png"), width =8, height = 5)
#ggsave(p5, file=paste0("plots/spatial_k_all_spp_wide_unnormalized.png"), width = 8, height = 12)


p3 = ggplot(coast) + geom_sf() + 
  geom_tile(data = dat_df, aes(x = X*1000, y = Y*1000, fill = norm_sd)) +
  scale_fill_viridis_c(name=expression("Normalized sd" ~ k)) + 
  theme_bw() + 
  xlab("") + ylab("") +
  facet_wrap(~fig_names, ncol = 3) +
  theme(strip.background = element_rect(fill="lightgray"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10))

