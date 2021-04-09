# ERP and topo plots for hs/per data

# load packages
library(tidyverse)
library(eegUtils)
library(haven)
library(here)
library(patchwork)
#'
#' read in data
dat <- read_csv(here("data", "paper_three", "hs_per_temp_fac_score_erp.csv")) %>%
  filter(comp %in% c("RC2", "RC3", "RC5", "RC7", "RC8", "RC9")) # retain only those components identified in study two
                                                         # that are sensitive to image valence
## peaks from hs_per_pca script:
# [1] "The maximum timepoint for RC2 is 370.800031605892"
# [2] "The maximum timepoint for RC3 is 736.346627189118"
# [3] "The maximum timepoint for RC5 is 134.269881522629"
# [4] "The maximum timepoint for RC7 is 194.868515015035"
# [5] "The maximum timepoint for RC12 is 102.993167462032"
# [6] "The maximum timepoint for RC8 is 257.421943136229"
# [7] "The maximum timepoint for RC17 is 1176.17541866626"
# [8] "The maximum timepoint for RC9 is 159.682211696864"

############################################
## electrodes for each component ordered chronologically
rc5_elec <- c("A29", "B26")
rc7_elec <- c("A29", "B26", "A31", "B30")
neg_rc8_elec <- c("B21", "B28")
pos_rc8_elec <- c("A29", "B26", "A31", "B30")
rc2_elec <- c("A29", "B26")
rc3_elec <- c("A29", "B26", "A26", "B23",
              "B28",
              "A30", "B27", "A25", "B22")
rc9_elec <- c("A28", "B25")
###########################################
## list of component sites for map function
component_list <- list("RC5",
                       "RC7",
                       "RC8",
                       "RC8",
                       "RC2",
                       "RC3",
                       "RC9")

# list of component selections
elec_selections <- list(rc5_elec,
                        rc7_elec,
                        neg_rc8_elec,
                        pos_rc8_elec,
                        rc2_elec,
                        rc3_elec,
                        rc9_elec)

# titles for each of the plots
titles <- list("134 ms Component Waveforms",
               "195 ms Component Waveforms",
               "Positive 257 ms Component Waveforms",
               "Negative 257 ms Component Waveforms",
               "371 ms Component Waveforms",
               "736 ms Component Waveforms",
               "160 ms Component Waveforms")

# change block variable to factor, reorder, and rename
dat$block <- factor(dat$block, levels = c("Neg_Watch",
                                          "Pos_Watch",
                                          "Neu_Watch"))

levels(dat$block) <- c("Negative",
                       "Positive",
                       "Neutral")

# iterate over each component and electrode selection to create ERP plots
watch_plots <-
  pmap(list(component_list,
            elec_selections,
            titles), ~ {
              # watch plots
              dat %>%
                filter(elec %in% ..2,
                       comp == ..1,
                       block %in% c("Positive", "Neutral", "Negative")) %>%
                group_by(block, ms) %>%
                mutate(ms = as.numeric(ms)) %>%
                summarise(mv = mean(mv)) %>%
                ggplot(., aes(ms, mv, color = block)) +
                geom_line(size = 1.1, alpha = 0.9) +
                geom_vline(xintercept = 0, linetype = "dashed") +
                geom_hline(yintercept = 0, linetype = "dashed") +
                labs(x = "Time (ms)",
                     y = expression(paste("Amplitude (",mu,"V)")),
                     #title = paste("Average", ..3),
                     color = "Block") +
                theme_classic() +
                theme(axis.title = element_text(size = 16),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 16),
                      legend.text = element_text(size = 12),
                      legend.key.size = unit(2, "line"),
                    #  plot.title = element_text(hjust = 0.5),
                      title = element_text(size = 16)) +
                scale_color_manual(breaks = c("Negative", "Positive", "Neutral"),
                                   values=c("red", "blue", "gray"))
            })

######### Topo plots ###########
# read in EEG coordinate data
elec_loc <- read_csv(here("data", "paper_two", "Equidistant Layout.csv"))
elec_loc <- elec_loc %>%
  rename("channel" = `channel name`) %>%
  filter(channel != "CMS", channel != "DRL") %>%
  select(-number)

elec_loc$radian_phi <- pi/180 * elec_loc$phi

## data frame with electrode coordinates
elec_loc <- elec_loc %>%
  mutate(x = theta * cos(radian_phi),
         y = theta * sin(radian_phi))

temp_dat_fac_scores <- read_csv(here("data", "paper_three", "hs_per_fac_score_dat.csv")) %>%
  rename("electrode" = "elec") %>%
  full_join(elec_loc, by = c("electrode" = "channel")) %>%
  filter(electrode != "EXG1",
         electrode != "EXG2") # not used in plotting, so omit

# define functions for topo plots
topo_plot_fun <- function(comp, elec) {
  values <-
    temp_dat_fac_scores %>%
    group_by(block, electrode) %>%
    summarise(comp = mean(comp, na.rm = TRUE)) %>%
    select(comp) %>%
    pull()

  maximum <- max(values)
  minimum <- min(values)

  neu <-
    temp_dat_fac_scores %>%
    filter(block == "Neu_Watch") %>%
    rename("amplitude" = comp) %>%
    topoplot(interp_limit = "head",
             highlights = elec,
             limits = c(minimum, maximum),
             scaling = 0.8) +
    ggtitle("Neutral") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  pos <-
    temp_dat_fac_scores %>%
    filter(block == "Pos_Watch") %>%
    rename("amplitude" = comp) %>%
    topoplot(interp_limit = "head",
             highlights = elec,
             limits = c(minimum, maximum),
             scaling = 0.8) +
    ggtitle("Positive") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  neg <-
    temp_dat_fac_scores %>%
    filter(block == "Neg_Watch") %>%
    rename("amplitude" = comp) %>%
    topoplot(interp_limit = "head",
             highlights = elec,
             limits = c(minimum, maximum),
             scaling = 0.8) +
    ggtitle("Negative") +
    theme(plot.title = element_text(hjust = 0.5))

  list(pos, neu, neg)
}

# iterate and create topo plots
topo_list <- map2(component_list, elec_selections, ~ topo_plot_fun(.x, .y))

# save objects
saveRDS(watch_plots, file = "data/paper_three/erp_plots.rds")
saveRDS(topo_list, file = "data/paper_three/topo_plots.rds")

# compose the final images
layout <- '
AAABBBCCC
DDDDDDDDD
'

rc5_component_plots <-
  topo_list[[1]][[2]] +
  topo_list[[1]][[1]] +
  topo_list[[1]][[3]] +
  watch_plots[[1]] +
plot_layout(design = layout,
            heights = c(1, 1.4),
            widths = 2,
            guides = "auto") +
  plot_annotation(title = "134 ms Component",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "rc5_plots.png"),
       plot = rc5_component_plots,
       height = 5,
       width = 8)

rc7_component_plots <-
  topo_list[[2]][[2]] +
  topo_list[[2]][[1]] +
  topo_list[[2]][[3]] +
  watch_plots[[2]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "195 ms Component",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "rc7_plots.png"),
       plot = rc7_component_plots,
       height = 5,
       width = 8)

neg_rc8_component_plots <-
  topo_list[[3]][[2]] +
  topo_list[[3]][[1]] +
  topo_list[[3]][[3]] +
  watch_plots[[3]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "Negative 257 ms Component (EPN)",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "neg_rc8_plots.png"),
       plot = neg_rc8_component_plots,
       height = 5,
       width = 8)

pos_rc8_component_plots <-
  topo_list[[4]][[2]] +
  topo_list[[4]][[1]] +
  topo_list[[4]][[3]] +
  watch_plots[[4]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "Positive 257 ms Component (EPN)",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "pos_rc8_plots.png"),
       plot = pos_rc8_component_plots,
       height = 5,
       width = 8)


rc2_component_plots <-
  topo_list[[5]][[2]] +
  topo_list[[5]][[1]] +
  topo_list[[5]][[3]] +
  watch_plots[[5]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "371 ms Component (LPP)",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "rc2_plots.png"),
       plot = rc2_component_plots,
       height = 5,
       width = 8)

rc3_component_plots <-
  topo_list[[6]][[2]] +
  topo_list[[6]][[1]] +
  topo_list[[6]][[3]] +
  watch_plots[[6]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "736 ms Component (LPP)",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "rc3_plots.png"),
       plot = rc3_component_plots,
       height = 5,
       width = 8)

rc9_component_plots <-
  topo_list[[7]][[2]] +
  topo_list[[7]][[1]] +
  topo_list[[7]][[3]] +
  watch_plots[[7]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "160 ms Component",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "rc9_plots.png"),
       plot = rc9_component_plots,
       height = 5,
       width = 8)
