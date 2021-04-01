# ERPs for non-PCA data and components

library(tidyverse)
library(eegUtils)
library(haven)
library(here)
library(patchwork)

# read in data
dat <- read_csv(here("data", "paper_three", "total_erp_dat.csv")) %>%
  filter(block %in% c("Pos_Watch", "Neg_Watch", "Neu_Watch"),
         ms < 2000)

############################################
## electrodes for each component ordered chronologically
neg_rc8_elec <- c("B21", "B28")
pos_rc8_elec <- c("A29", "B26", "A31", "B30")
rc2_elec <- c("A29", "B26")
rc3_elec <- c("A29", "B26", "A26", "B23",
              "B28",
              "A30", "B27", "A25", "B22")
###########################################
# list of component selections
elec_selections <- list(neg_rc8_elec,
                        pos_rc8_elec,
                        rc2_elec,
                        rc3_elec
                        )

# change block variable to factor, reorder, and rename
dat$block <- factor(dat$block, levels = c("Neg_Watch",
                                          "Pos_Watch",
                                          "Neu_Watch"))

levels(dat$block) <- c("Negative",
                       "Positive",
                       "Neutral")

# iterate over each component and electrode selection to create ERP plots
watch_plots <-
  map(elec_selections, ~ {
              # watch plots
              dat %>%
                select(pid:ms, .x) %>%
                pivot_longer(cols = .x,
                             names_to = "elec",
                             values_to = "mv") %>%
                filter(block %in% c("Positive", "Neutral", "Negative")) %>%
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

dat_loc <- dat %>%
  pivot_longer(cols = c(A1:EXG2),
               names_to = "electrode",
               values_to = "mv") %>%
  full_join(elec_loc, by = c("electrode" = "channel")) %>%
  filter(electrode != "EXG1",
         electrode != "EXG2") # not used in plotting, so omit

# define functions for topo plots
topo_plot_fun <- function(elec) {
  values <-
    dat_loc %>%
    group_by(block, electrode) %>%
    summarise(mv = mean(mv, na.rm = TRUE)) %>%
    select(mv) %>%
    pull()

  maximum <- max(values)
  minimum <- min(values)

  neu <-
    dat_loc %>%
    filter(block == "Neutral") %>%
    rename("amplitude" = mv) %>%
    topoplot(interp_limit = "head",
             highlights = elec,
             limits = c(minimum, maximum),
             scaling = 0.8) +
    ggtitle("Neutral") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  pos <-
    dat_loc %>%
    filter(block == "Positive") %>%
    rename("amplitude" = mv) %>%
    topoplot(interp_limit = "head",
             highlights = elec,
             limits = c(minimum, maximum),
             scaling = 0.8) +
    ggtitle("Positive") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  neg <-
    dat_loc %>%
    filter(block == "Negative") %>%
    rename("amplitude" = mv) %>%
    topoplot(interp_limit = "head",
             highlights = elec,
             limits = c(minimum, maximum),
             scaling = 0.8) +
    ggtitle("Negative") +
    theme(plot.title = element_text(hjust = 0.5))

  list(pos, neu, neg)
}

# iterate and create topo plots
topo_list <- map(elec_selections, ~ topo_plot_fun(.x))

# compose the final images
layout <- '
AAABBBCCC
DDDDDDDDD
'

neg_rc8_component_plots <-
  topo_list[[1]][[2]] +
  topo_list[[1]][[1]] +
  topo_list[[1]][[3]] +
  watch_plots[[1]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "CPz and Pz Activity",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "pre_PCA", "neg_rc8_plots.png"),
       plot = neg_rc8_component_plots,
       height = 5,
       width = 8)

pos_rc8_component_plots <-
  topo_list[[2]][[2]] +
  topo_list[[2]][[1]] +
  topo_list[[2]][[3]] +
  watch_plots[[2]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "PO7, PO8, O1, and O2 Activity",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "pre_PCA", "pos_rc8_plots.png"),
       plot = pos_rc8_component_plots,
       height = 5,
       width = 8)

rc2_component_plots <-
  topo_list[[3]][[2]] +
  topo_list[[3]][[1]] +
  topo_list[[3]][[3]] +
  watch_plots[[3]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "PO7 and PO8 Activity",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "pre_PCA", "rc2_plots.png"),
       plot = rc2_component_plots,
       height = 5,
       width = 8)

rc3_component_plots <-
  topo_list[[4]][[2]] +
  topo_list[[4]][[1]] +
  topo_list[[4]][[3]] +
  watch_plots[[4]] +
  plot_layout(design = layout,
              heights = c(1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "PO7, PO8, PO3, PO4, P1, P2,\nP5/P3, P6/P4 and Pz Activity",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(here("images", "paper_3", "ERP Topo Images", "pre_PCA", "rc3_plots.png"),
       plot = rc3_component_plots,
       height = 5,
       width = 8)