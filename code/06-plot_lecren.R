library(ggplot2)
library(tidyverse)
library(sf)
library(patchwork)
library(colorspace)

spec <- readRDS("survey_data/spec_list.rds")


# generate map
map_data <- rnaturalearth::ne_states(
  returnclass = "sf", country = "united states of america")
# Crop the polygon for plotting and efficiency:
usa <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = -145, ymin = 25, xmax = -117, ymax = 48.5))))
coast <- sf::st_transform(usa, crs = 3157)

# create empty list
lecren_list = list()


for(i in 1:length(spec)) {
# load data
biomass_df <- readRDS(paste0("output/biomass/spatial_preds/", spec[i], "_pred_index.rds")) %>% 
   mutate(biomass = exp(est), zeta_s_fyear = fyear)
d <- readRDS(paste0("output/lecren/",spec[i],"_pred_lecren.rds"))

# merge dataframes
d = left_join(d, biomass_df[,c("year", "X", "Y", "biomass")])

# organize data frame
d_est = d %>% 
  select(-contains("zeta_s_fyear_scaled")) %>% 
  mutate(global_mean = mean(exp(est)), condition = exp(est)) %>% 
  mutate(common_name = replace(common_name, common_name == "lingcod" & Y < 4440.148, "lingcod_south"),
         common_name = replace(common_name, common_name == "lingcod" & Y >= 4440.148, "lingcod_north")) %>% 
  group_by(year, common_name) %>% 
  mutate(tot_biomass = sum(biomass), b_weight = biomass/tot_biomass) %>% 
  ungroup() %>% 
  group_by(X, Y) %>% 
  # mean_lecren = for each spatial coordinate, take the average k_f across all years
  # anom_lecren = for each spatial coordinate, it is the difference between k_f for that year and the average k_f across all years (mean_k_f) at that specific spatial coordinate
  mutate(mean_lecren = mean(condition)) %>% 
  mutate(anom_lecren = condition - mean_lecren)

# save dataframe
lecren_list[[i]] = d_est

}


pal = c("#69b27b","#c86dbb","#7a87d4","#76bc50","#d99431","#ac57c5", "#51bdc6","#d85e73","#b2a84f","#ca764f")

# bind k_f predictions for each species into one df
lecren_df = bind_rows(lecren_list) 
#saveRDS(lecren_df, paste0("output/index/lecren_df.rds"))

# calculate condition indices for each species
condition_index = lecren_df %>% 
  group_by(common_name, year) %>% 
  summarise(condition_index = sum(condition*b_weight), unweighted_cond_ind = mean(condition))  %>% 
  pivot_longer(c(condition_index, unweighted_cond_ind), names_to = "biomass_weighted", values_to = "condition_index") %>% 
  mutate(biomass_weighted=replace(biomass_weighted, biomass_weighted== "condition_index", "biomass-weighted"),
         biomass_weighted=replace(biomass_weighted, biomass_weighted== "unweighted_cond_ind", "unweighted")) %>% 
  group_by(common_name, biomass_weighted) %>% 
  mutate(global_avg = mean(condition_index))
condition_index$fig_names <- str_to_sentence(sub("_", " ", condition_index$common_name)) # replace spaces with underscores in common name
#saveRDS(condition_index, paste0("output/index/condition_index.rds"))


# merge condition and growth rate dataframes
cond_growth_index_df = left_join(condition_index, index_df[,c("common_name", "year", "biomass_weighted", "growth_index")])


# plot 1 - map of condition 

# merge lingcod subpopulations into one
lecren_map_df = lecren_df %>% 
  mutate(common_name = replace(common_name, common_name == "lingcod_south", "lingcod"),
         common_name = replace(common_name, common_name == "lingcod_north", "lingcod"))

# assign fish information
fish_info = as.data.frame(cbind(common_name = c("aurora_rockfish","darkblotched_rockfish","splitnose_rockfish","dover_sole", "pacific_sanddab","petrale_sole", "lingcod","pacific_hake", "sablefish"), fish_taxa = c(rep("Rockfish", 3), rep("Flatfish", 3), rep("Roundfish", 3)), depth_dist = c("slope", "shelf-slope", "slope", "shelf-slope", "shelf-slope", "shelf", "shelf", "shelf-slope", "shelf-slope")))

# merge data and fish information
lecren_map_df = left_join(lecren_map_df, fish_info)
lecren_map_df$fig_names <- str_to_sentence(sub("_", " ", lecren_map_df$common_name)) 
lecren_map_df = lecren_map_df %>%
  mutate(fig_names = factor(fig_names,
                            levels = c("Dover sole","Pacific sanddab",  "Petrale sole", "Aurora rockfish","Darkblotched rockfish", "Splitnose rockfish", "Lingcod", "Pacific hake",  "Sablefish")))


# spatial map of body condition
p1 = ggplot(coast) + geom_sf() + 
  geom_tile(data = lecren_map_df, aes(x = X*1000, y = Y*1000, fill = mean_lecren)) +
  theme_bw() + 
  xlab("") + ylab("") +
  facet_wrap(~fig_names) +
  scale_fill_continuous_divergingx(palette = "RdYlBu", mid = 1, rev = TRUE) +
  theme(strip.background =element_rect(fill="lightgray"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8.5),
        axis.title=element_text(size=13),
        axis.text=element_text(size=10),
        legend.text=element_text(size=8),
        legend.title = element_text(size = 10)) +
  labs(fill = expression("Mean" ~ "K"^n*""[s]))

# plot 2 - relationship between condition index and growth index
p4 = cond_growth_index_df[cond_growth_index_df$biomass_weighted == "biomass-weighted",] %>%
  ggplot(aes(x = condition_index, y = growth_index)) +
  geom_point(aes(fill = year), shape = 21, col = "black", size = 5) +
  scale_fill_gradient(low = "gold", high = "red3") +
  theme_classic() +
  facet_wrap(~fig_names, scales = "free_y") +
  labs(x = "Biomass-Weighted Condition Index", y = "Biomass-Weighted Growth Index", fill = "Year") +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 10))

ggsave(p1, file=paste0("plots/final_figures/fig4.jpeg"), width = 170, height = 255, dpi = 500, units = "mm")
ggsave(p2, file=paste0("plots/lecren/condition_growth_corr.jpeg"), width =9, height = 5)




