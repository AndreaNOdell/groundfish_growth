library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)

# plot 1 - Growth and Condition Indices --------------------

# load data sets
growth_index = readRDS("output/index/growth_index.rds")
condition_index = readRDS("output/index/condition_index.rds")
lecren_df = readRDS("output/index/lecren_df.rds")
growth_df = readRDS("output/index/dat_df.rds")
condition_sigma = readRDS("output/index/condition_sigma.rds") # run code/08-coefs.R

# Merge two datasets
growth_cond_df = left_join(growth_index[,c("common_name", "year", "biomass_weighted", "growth_index")], condition_index[,c("common_name", "year", "biomass_weighted", "condition_index")])
growth_cond_df$fig_names <- str_to_sentence(sub("_", " ", growth_cond_df$common_name)) # create figure names
growth_cond_df = growth_cond_df %>% 
  mutate(fig_names = str_replace(fig_names, "Lingcod north", "Lingcod (north)"),
         fig_names = str_replace(fig_names, "Lingcod south", "Lingcod (south)"))

# label species by fish type
fish_info = as.data.frame(cbind(fig_names = c("Aurora rockfish","Darkblotched rockfish","Splitnose rockfish","Dover sole", "Pacific sanddab","Petrale sole", "Lingcod (north)", "Lingcod (south)", "Pacific hake", "Sablefish"), fish_taxa = c(rep("Rockfish", 3), rep("Flatfish", 3), rep("Roundfish", 4)), depth_dist = c("slope", "shelf-slope", "slope", "shelf-slope", "shelf-slope", "shelf", "shelf", "shelf", "shelf-slope", "shelf-slope"), lifespan = c("long-lived", "long-lived", "long-lived", "long-lived", "short-lived", "short-lived", "short-lived", "short-lived", "short-lived", "long-lived"), maturity = c("late-maturation", "late-maturation", "late-maturation", "early-maturation", "early-maturation", "early-maturation", "early-maturation", "early-maturation","early-maturation","early-maturation"))) %>% 
  mutate(life_history = paste(lifespan, maturity))

# add fish type to dataset and convert to long format for plotting
growth_cond_df = left_join(growth_cond_df, fish_info) %>% 
  filter(biomass_weighted == "biomass-weighted")  %>% 
  pivot_longer(cols = growth_index:condition_index, names_to = "index", values_to = "value") %>% 
  mutate(index_fct = factor(index, levels = c("growth_index", "condition_index")))
growth_cond_df$index_fct =  sub("_", " ", growth_cond_df$index_fct)

# create column with growth rate from most recent year to calculate relative change (to show which years a species growth rate increases versus decreases.)
relchange_df = growth_cond_df %>% 
  filter(index == "growth_index") %>% 
  select(common_name, biomass_weighted, value, index, year) %>% 
  arrange(common_name, year) %>% 
  group_by(common_name) %>% 
  mutate(prev_yr = lag(value), start_val = first(value)) 

# merge relative change information to data set
growth_cond_df = left_join(growth_cond_df, relchange_df)


# order species for figure and calculate relative change
growth_cond_df = growth_cond_df %>% 
  mutate(fig_names = factor(fig_names,
                            levels = c("Lingcod (north)", "Lingcod (south)" ,"Petrale sole", "Dover sole", "Pacific sanddab", "Darkblotched rockfish",  "Pacific hake", "Sablefish", "Aurora rockfish", "Splitnose rockfish"))) %>% 
  group_by(index_fct) %>% 
  mutate(growth_relchange = (value - prev_yr)/prev_yr, growth_change = 100*((value - start_val)/start_val)) %>% 
  replace(is.na(.), 0) 

# color palette 
pal = c("#69b27b","#c86dbb","#7a87d4","#76bc50","#d99431","#ac57c5", "#51bdc6","#d85e73","#b2a84f","#ca764f")

# Create growth plot
p1 = ggplot(growth_cond_df[growth_cond_df$index_fct == "growth index",], aes( x = year, y = growth_change, col = fig_names)) +
  geom_hline(aes(yintercept = 0), linewidth = 0.5, col = "darkgray") +
  geom_line(aes(col = fig_names), size = 1, alpha =0.8) + 
  geom_point(aes(shape = life_history, fill = fig_names), col = "gray45", size = 1.75, alpha = 0.9) + 
  ggh4x::facet_grid2(depth_dist ~ index_fct, independent = "y", scales = "free", space = "free_x") +
  theme_bw() +
  scale_color_manual(values= pal) +
  scale_fill_manual(values= pal, guide = "none") +
  scale_shape_manual(values=c(21, 24, 22)) +
  theme(strip.background =element_rect(fill="lightgray"),
        strip.text = element_text(size = 9),
        strip.text.y = element_blank(),
        axis.title=element_text(size=9),
        axis.text=element_text(size=8)) +
  labs(col = "Species", x = "Year", y = expression(Percent~change~"in"~k~relative~to~t[1])) + 
  guides(shape = "none", col = "none", linetype = "none")

# create condition plot
p2 = ggplot(growth_cond_df[growth_cond_df$index_fct == "condition index",], aes( x = year, y = value)) +
  geom_hline(aes(yintercept = 1), linewidth = 0.5, col = "darkgray") +
  geom_line(aes(col = fig_names), size = 1, alpha = 0.8) +
  geom_point(aes(shape = life_history, fill = fig_names), col = "gray45", size = 1.75, alpha = 0.9) + 
  ggh4x::facet_grid2(depth_dist ~ index_fct, independent = "y", scales = "free", space = "free_x") +
  theme_bw() +
  scale_color_manual(values= pal) +
  scale_fill_manual(values= pal, guide = "none") +
  scale_shape_manual(values=c(21, 24, 22)) +
  #lims(y = c(0.93,1.05)) +
  theme(strip.background =element_rect(fill="lightgray"),
        strip.text = element_text(size = 9),
        axis.title=element_text(size=9),
        axis.text=element_text(size=8),
        legend.text=element_text(size=8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(4.5, 'mm'),
        guides(fill = "none")) +
  labs(y = expression(Relative~condition~factor~K[n]), col = "Species", x = "Year", shape = "Life history")


plot = ggarrange(p1,p2, common.legend = TRUE, legend = "right", font.label = list(size = 11))

ggsave(plot, file=paste0("plots/final_figures/fig5.jpeg"), width =170, height = 110, dpi = 400, units = "mm")



# plot 2 growth vs condition relationship local ----------------------------
growth_df$zeta_s_fyear = as.integer(growth_df$zeta_s_fyear)
lecren_df = lecren_df %>% 
  rename(zeta_s_fyear = year) %>% 
  mutate(common_name = replace(common_name, common_name == "lingcod_south", "lingcod"),
         common_name = replace(common_name, common_name == "lingcod_north", "lingcod"))

full_dat = left_join(growth_df[,c( "zeta_s_fyear", "X", "Y", "k_f", "mean_k_f", "common_name", "fig_names")], 
          lecren_df[,c("zeta_s_fyear",  "X", "Y", "common_name", "biomass", "condition", "tot_biomass", "b_weight", "mean_lecren")])
full_dat$fig_names =  str_to_sentence(sub("_", " ", full_dat$common_name))

p4 = full_dat %>% 
  distinct(X, Y, fig_names, mean_lecren, mean_k_f) %>% 
ggplot(aes(x = mean_lecren, y = mean_k_f)) +
  geom_point() +
  facet_wrap(~fig_names) +
  theme_bw() +
labs(alpha = "biomass weighting", x = expression(Average~body~condition~"in"~location~italic(s)), y = expression(Average~growth~rate~"in"~location~italic(s)))
  
ggsave(p4, file=paste0("plots/final_figures/figS2.jpeg"), width =170, height = 100, dpi = 400, units = "mm")



# Plot 3 Spatial vs Temporal c.v. Growth and Condition ----------------
condition_cv = lecren_df %>% 
  group_by(common_name, year) %>% # calculate mean for each year (temporal variability)
  mutate(year_lecren = mean(condition)) %>% 
  ungroup()  %>%
  mutate(fig_names = sub("_", " ", lecren_df$common_name)) %>% 
  group_by(fig_names) %>% 
  summarise(spat_cv = sd(mean_lecren)/mean(mean_lecren), temp_cv = sd(year_lecren)/mean(year_lecren)) %>%
  mutate(index = "condition index") 

# merge growth and condition data
growth_cv = full_dat %>% 
  group_by(fig_names) %>% 
  summarise(st_cv = sd(k_f)/mean(k_f), 
            spat_cv = sd(mean_k_f)/mean(mean_k_f))

condition_sigma$fig_names = str_to_sentence(condition_sigma$fig_names)

growth_condition_var = left_join(growth_cv, condition_sigma[,2:4]) %>% 
  pivot_longer(!fig_names, names_to = "metric", values_to = "value") %>% 
  mutate(index = case_when(metric == "st_cv" ~ "Growth index",
                           metric == "spat_cv" ~ "Growth index",
                           metric == "sigma_O" ~ "Condition index",
                           metric == "sigma_E" ~ "Condition index"),
         metric_fig = case_when(metric == "st_cv" ~ "spatiotemporal",
                                metric == "spat_cv" ~ "spatial",
                                metric == "sigma_O" ~ "spatial",
                                metric == "sigma_E" ~ "spatiotemporal")) %>% 
  select(fig_names, value, index, metric_fig) %>% 
  pivot_wider(names_from = metric_fig, values_from = value) 
#growth_condition_var$fig_names = str_to_sentence(growth_condition_var$fig_names)

fish_info[fish_info$fig_names == "Lingcod (south)",]$fig_names = "Lingcod"

growth_condition_var = left_join(growth_condition_var, fish_info) %>% 
  mutate(index = factor(index, levels = c("Growth index", "Condition index"))) %>% 
  mutate(fig_names = factor(fig_names,
                            levels = c("Lingcod", "", "Petrale sole", "Dover sole", "Pacific sanddab", "Darkblotched rockfish",  "Pacific hake", "Sablefish", "Aurora rockfish", "Splitnose rockfish")))

pal_onelingcod = c("#69b27b","#7a87d4","#76bc50","#d99431","#ac57c5", "#51bdc6","#d85e73","#b2a84f","#ca764f")

p1 = growth_condition_var %>% 
  filter(index == "Growth index") %>% 
ggplot(aes(x = spatial, y = spatiotemporal)) +
  geom_abline(slope = 1, intercept = 0, col = "darkgray") +
  geom_point(aes(fill = fig_names, shape = life_hist), size = 4, col = "gray31") +
  scale_fill_manual(values= pal_onelingcod) +
  scale_shape_manual(values=c(21, 24, 22)) +
  theme_bw() +
  facet_wrap(~ index, labeller = labeller(index = 
                                            c("Growth index" = "Growth rate"))) +
  theme(strip.background =element_rect(fill="lightgray"),
        strip.text = element_text(size = 10),
        axis.title=element_text(size=9),
        axis.text=element_text(size=8),
        legend.text=element_text(size=8),
        legend.title = element_text(size=10)) +
  guides(fill= "none", shape = "none") +
  labs(x = "Range of spatial variability (c.v.)", y = "Range of spatiotemporal variability (c.v.)") +
  xlim(0, NA) + ylim(0, NA)

p2 = growth_condition_var %>% 
  filter(index == "Condition index") %>% 
  ggplot(aes(x = spatial, y = spatiotemporal)) +
  geom_abline(slope = 1, intercept = 0, col = "darkgray") +
  geom_point(aes(fill = fig_names, shape = life_hist), size = 4) +
  scale_fill_manual(values= pal_onelingcod) +
  scale_shape_manual(values=c(21, 24, 22))+
  theme_bw() +
  facet_wrap(~ index, labeller = labeller(index = 
                                            c("Condition index" = "Body condition"))) +
  theme(strip.background =element_rect(fill="lightgray"),
        strip.text = element_text(size = 10),
        axis.title=element_text(size=9),
        axis.text=element_text(size=8),
        legend.text=element_text(size=8),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=10)) +
  guides(fill=guide_legend(override.aes=list(shape=21, alpha = 1), order = 1),
         shape = guide_legend(override.aes=list(alpha = 1))) +
  labs(x = expression(paste("Range of spatial variability (", Sigma[omega], ")")), y = expression(paste("Range of spatiotemporal variability (", Sigma[epsilon], ")")), fill = "Species", shape = "Life history") +
  xlim(0, NA) + ylim(0, NA)

plot = ggarrange(p1,p2, labels = c("a)", "b)"), common.legend = TRUE, legend = "right", font.label = list(size = 11))

ggsave(plot, file=paste0("plots/final_figures/fig2.jpeg"), width =170, height = 80, dpi = 400, units = "mm")


