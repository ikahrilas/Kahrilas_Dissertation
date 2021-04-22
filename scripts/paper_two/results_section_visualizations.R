# results section visualizations

## read in packages
library(tidyverse)
library(gridExtra)
library(patchwork)
library(eegUtils)
library(here)
library(grid)

# read in data
fac_score_dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-18.csv")

glimpse(fac_score_dat)

# change data to long form and reorder variables
fac_score_dat_long <- fac_score_dat %>%
  pivot_longer(cols = RC2:RC12,
               names_to = "component",
               values_to = "fac_score") %>%
  relocate(pid, block, component, fac_score, everything())

# change component variable to a factor type
fac_score_dat_long$component <- factor(fac_score_dat_long$component,
                                       levels = c("RC5", "RC11", "RC12", "RC3", "RC2"))

# rename the factor levels
levels(fac_score_dat_long$component) <- c("RC5", "RC11", "neg_RC12", "RC3", "RC2")

levels(fac_score_dat_long$component) <- c("124 ms Pos Peak",
                                          "162 ms Pos Peak",
                                          "259 ms Neg Peak",
                                          "381 ms Pos Peak",
                                          "740 ms Pos Peak")

# reorder block factor for proper plotting of legend
fac_score_dat_long$block <- as.factor(fac_score_dat_long$block)
levels(fac_score_dat_long$block) <- c("Negative Decrease",
                                      "Negative Increase",
                                      "Negative Watch",
                                      "Neutral Watch",
                                      "Positive Decrease",
                                      "Positive Increase",
                                      "Positive Watch")

# envisioning 3 grouped box plots, one with watch conditions, one with positive conditions
# and one with negative conditions, with simple topo plots with electrode regions highlighted.
# use this stackoverflow post: https://stackoverflow.com/questions/29263046/how-to-draw-the-boxplot-with-significant-level
# for significance bars

# make separate data frames for each facet, then create facet variable for actual facetting
watch_cases <- fac_score_dat_long %>%
  filter(str_detect(block, "Watch")) %>%
  mutate(facet = 1, # facet variable
         block = factor(block, levels = c("Negative Watch", "Neutral Watch", "Positive Watch")))
pos_cases <- fac_score_dat_long %>%
  filter(str_detect(block, "Positive")) %>%
  mutate(facet = 2, # facet variable
         block = factor(block, levels = c("Positive Decrease", "Positive Watch", "Positive Increase")))
neg_cases <- fac_score_dat_long %>%
  filter(str_detect(block, "Negative")) %>%
  mutate(facet = 3, # facet variable
         block = factor(block, levels = c("Negative Decrease", "Negative Watch", "Negative Increase")))

# first plot the neutral conditions
p_1 <-
  ggplot(watch_cases, aes(x = component, y = fac_score)) +
  geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
  scale_fill_manual(values = c(`Negative Watch` = "purple",
                               `Neutral Watch` = "grey",
                               `Positive Watch` ="blue")) +
  geom_boxplot(aes(group = interaction(block, component)), fatten = 0.75, outlier.size = 1,
               width = 0.2, fill = "white", position = position_dodge(width = .75)) +
  facet_wrap(~ facet, ncol = 1) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(x = NULL,
       y = expression(paste("Amplitude (",mu,"V)")),
       fill = "Block") +
  coord_cartesian(
    xlim = NULL,
    ylim = c(-4.5, 6),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  annotate(geom = "segment", # 170 ms peak annotations
           x = c(1.75, 2.025),
           xend = c(1.975, 2.25),
           y = 5.6,
           yend = 5.6,
           color = "black") +
  annotate(geom = "segment",
           x = c(1.75, 1.975, 2.025, 2.25),
           xend = c(1.75, 1.975, 2.025, 2.25),
           y = 5.6,
           yend = c(5.2, 4.1, 4.1, 5.1),
           color = "black") +
  annotate(geom = "segment",
           x = 1.75,
           xend = 2.25,
           y = -1.85,
           yend = -1.85,
           color = "black") +
  annotate(geom = "segment",
           x = c(1.75, 2.25),
           xend = c(1.75, 2.25),
           y = c(-1.85, -1.85),
           yend = c(-0.25, -1.05)) +
  annotate(geom = "text",
           x = c(1.8625, 2, 2.1375),
           y = c(5.7, -2.9, 5.7),
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 250 ms negative peak annotations
           x = c(2.75, 2.75),
           xend = c(3, 3.25),
           y = c(2.7, -3.9),
           yend = c(2.7, -3.9),
           color = "black") +
  annotate(geom = "segment",
           x = c(2.75, 3, 2.75, 3.25),
           xend = c(2.75, 3, 2.75, 3.25),
           y = c(2.7, 2.7, -3.9, -3.9),
           yend = c(2.4, 1.5, -3.65, -3.1),
           color = "black") +
  annotate(geom = "text",
           x = c(2.875, 3),
           y = c(2.7, -4.9),
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 375 ms positive peak annotations
           x = c(3.75, 4.025, 3.75),
           xend = c(4, 4.275, 4.275),
           y = c(3.45, 3.45, -1.5),
           yend = c(3.45, 3.45, -1.5),
           color = "black") +
  annotate(geom = "segment",
           x = c(3.75, 4, 4.025, 4.275, 3.75, 4.275),
           xend = c(3.75, 4, 4.025, 4.275, 3.75, 4.275),
           y = c(3.45, 3.45, 3.45, 3.45, -1.5, -1.5),
           yend = c(3.25, 1.8, 1.8, 2, -0.4, -0.5),
           color = "black") +
  annotate(geom = "text",
           x = c(3.875, 4.15, 4.0125),
           y = c(3.55, 3.55, -2.5),
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 800 ms peak annotations
           x = c(4.75, 5.025, 4.75),
           xend = c(4.975, 5.25, 5.25),
           y = c(5, 5, -1.7),
           yend = c(5, 5, -1.7),
           color = "black") +
  annotate(geom = "segment",
           x = c(4.75, 4.975, 5.025, 5.25, 4.75, 5.25),
           xend = c(4.75, 4.975, 5.025, 5.25, 4.75, 5.25),
           y = c(5, 5, 5, 5, -1.7, -1.7),
           yend = c(4.7, 2.85, 2.85, 3.6, -0.2, -0.7)) +
  annotate(geom = "text",
           x = c(4.8625, 5.1375, 5),
           y = c(5.1, 5.1, -2.7),
           label = "*",
           size = 5)

p_2 <-
  ggplot(pos_cases, aes(x = component, y = fac_score)) +
  geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
  scale_fill_manual(values = c(`Positive Decrease` = "cadetblue1",
                               `Positive Watch` = "blue",
                               `Positive Increase` = "springgreen")) +
  geom_boxplot(aes(group = interaction(block, component)), fatten = 0.75, outlier.size = 1,
               width = 0.2, fill = "white", position = position_dodge(width = .75)) +
  facet_wrap(~ facet, ncol = 1) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
  ) +
  labs(x = NULL,
       y = expression(paste("Amplitude (",mu,"V)")),
       fill = "Block") +
  ylim(-4.5, 6) +
  annotate(geom = "segment", # 250 ms negative peak annotations
           x = c(2.75, 3),
           xend = c(3.25, 3.25),
           y = c(-3.5, 1.75),
           yend = c(-3.5, 1.75),
           color = "black") +
  annotate(geom = "segment",
           x = c(2.75, 3.25, 3, 3.25),
           xend = c(2.75, 3.25, 3, 3.25),
           y = c(-3.5, -3.5, 1.75, 1.75),
           yend = c(-2.5, -3.25, 1.5, 1.5),
           color = "black") +
  annotate(geom = "text",
           x = c(3, 3.125),
           y = c(-4.5, 1.85),
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 375 ms peak annotations
           x = c(4.75, 5),
           xend = c(5.25, 5.25),
           y = c(-1.2, 2.65),
           yend = c(-1.2, 2.65),
           color = "black") +
  annotate(geom = "segment",
           x = c(4.75, 5.25, 5, 5.25),
           xend = c(4.75, 5.25, 5, 5.25),
           y = c(-1.2, -1.2, 2.65, 2.65),
           yend = c(-.95, -.7, 2.05, 2.3)) +
  annotate(geom = "text",
           x = c(5, 5.1275),
           y = c(-2.2, 2.75),
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 800 ms peak annotations
           x = c(5.75),
           xend = c(6.25),
           y = c(-1),
           yend = c(-1),
           color = "black") +
  annotate(geom = "segment",
           x = c(5.75, 6.25),
           xend = c(5.75, 6.25),
           y = c(-1, -1),
           yend = c(-0.3, -0.7)) +
  annotate(geom = "text",
           x = c(6),
           y = c(-2),
           label = "*",
           size = 5)

p_3 <-
  ggplot(neg_cases, aes(x = component, y = fac_score)) +
  geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
  scale_fill_manual(values = c(`Negative Decrease` = "plum",
                               `Negative Watch` = "purple",
                               `Negative Increase` = "red")) +
  geom_boxplot(aes(group = interaction(block, component)), fatten = 0.75, outlier.size = 1,
               width = 0.2, fill = "white", position = position_dodge(width = .75)) +
  facet_wrap(~ facet, ncol = 1) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x = element_text(size = 16)
  ) +
  labs(x = NULL,
       y = expression(paste("Amplitude (",mu,"V)")),
       fill = "Block") +
  ylim(-4.5, 6) +
  coord_cartesian(clip = "off") +
  annotate(geom = "segment", # 125 ms peak annotations
           x = c(1),
           xend = c(1.25),
           y = 5.85,
           yend = 5.85,
           color = "black") +
  annotate(geom = "segment",
           x = c(1, 1.25),
           xend = c(1, 1.25),
           y = 5.85,
           yend = c(4.4, 5.5),
           color = "black") +
  annotate(geom = "text",
           x = c(1.125),
           y = 5.95,
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 250 ms negative peak annotations
           x = c(3),
           xend = c(3.25),
           y = 2.8,
           yend = 2.8,
           color = "black") +
  annotate(geom = "segment",
           x = c(3, 3.25),
           xend = c(3, 3.25),
           y = 2.8,
           yend = c(2.45, 1.0),
           color = "black") +
  annotate(geom = "text",
           x = c(3.125),
           y = 2.9,
           label = "*",
           size = 5)

##########################################
# plots grouped by regulation condition #
##########################################

# make separate data frames for each facet, but this time by increase/decrease
inc_cases <- fac_score_dat_long %>%
  filter(str_detect(block, "Increase")) %>%
  mutate(facet = 2, # facet variable
         block = factor(block, levels = c("Negative Increase", "Positive Increase")))

dec_cases <- fac_score_dat_long %>%
  filter(str_detect(block, "Decrease")) %>%
  mutate(facet = 3, # facet variable
         block = factor(block, levels = c("Negative Decrease", "Positive Decrease")))

 p_4 <-
  ggplot(inc_cases, aes(x = component, y = fac_score)) +
  geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
  scale_fill_manual(values = c(`Negative Increase` = "red",
                               `Positive Increase` = "springgreen")) +
  geom_boxplot(aes(group = interaction(block, component)), fatten = 0.75, outlier.size = 1,
               width = 0.2, fill = "white", position = position_dodge(width = .75)) +
  facet_wrap(~ facet, ncol = 1) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(x = NULL,
       y = expression(paste("Amplitude (",mu,"V)")),
       fill = "Block") +
  coord_cartesian(
    xlim = NULL,
    ylim = c(-4.5, 6),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  annotate(geom = "segment", # 250 ms negative peak annotations
           x = 2.815,
           xend = 3.185,
           y = 1.6,
           yend = 1.6,
           color = "black") +
  annotate(geom = "segment",
           x = c(2.815, 3.185),
           xend = c(2.815, 3.185),
           y = c(1.6, 1.6),
           yend = c(0.7, 1.35),
           color = "black") +
  annotate(geom = "text",
           x = 3,
           y = 1.7,
           label = "*",
           size = 5) +
  annotate(geom = "segment", # 250 ms positive peak annotations
           x = 3.815,
           xend = 4.185,
           y = 4.1,
           yend = 4.1,
           color = "black") +
  annotate(geom = "segment",
             x = c(3.815, 4.185),
             xend = c(3.815, 4.185),
             y = c(4.1, 4.1),
             yend = c(3.4, 3.8),
             color = "black") +
  annotate(geom = "text",
             x = 4,
             y = 4.2,
             label = "*",
             size = 5) +
    annotate(geom = "segment", # 800 ms positive peak annotations
             x = 5.815,
             xend = 6.185,
             y = 4.3,
             yend = 4.3,
             color = "black") +
    annotate(geom = "segment",
             x = c(5.815, 6.185),
             xend = c(5.815, 6.185),
             y = c(4.3, 4.3),
             yend = c(4.05, 3.8),
             color = "black") +
    annotate(geom = "text",
             x = 6,
             y = 4.4,
             label = "*",
             size = 5)

p_5 <-
  ggplot(dec_cases, aes(x = component, y = fac_score)) +
    geom_violin(aes(fill = block),
              position = position_dodge(width = .75),
              trim = TRUE) +
    scale_fill_manual(values = c(`Negative Decrease` = "plum",
                                 `Positive Decrease` = "cadetblue1")) +
    geom_boxplot(aes(group = interaction(block, component)), fatten = 0.75, outlier.size = 1,
                 width = 0.2, fill = "white", position = position_dodge(width = .75)) +
    facet_wrap(~ facet, ncol = 1) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
    ) +
    labs(x = "Component",
         y = expression(paste("Amplitude (",mu,"V)")),
         fill = "Block") +
    ylim(-4.5, 6) +
    annotate(geom = "segment", # 250 ms negative peak annotations
           x = 2.815,
           xend = 3.185,
           y = 1.8,
           yend = 1.8,
           color = "black") +
    annotate(geom = "segment",
             x = c(2.815, 3.185),
             xend = c(2.815, 3.185),
             y = c(1.8, 1.8),
             yend = c(1.3, 1.55),
             color = "black") +
    annotate(geom = "text",
             x = 3,
             y = 1.9,
             label = "*",
             size = 5) +
    annotate(geom = "segment", # 375 ms negative peak annotations
             x = 4.815,
             xend = 5.185,
             y = 3.2,
             yend = 3.2,
             color = "black") +
    annotate(geom = "segment",
             x = c(4.815, 5.185),
             xend = c(4.815, 5.185),
             y = c(3.2, 3.2),
             yend = c(2.9, 2.9),
             color = "black") +
    annotate(geom = "text",
             x = 5,
             y = 3.3,
             label = "*",
             size = 5) +
    annotate(geom = "segment", # 800 ms negative peak annotations
             x = 5.815,
             xend = 6.185,
             y = 4.5,
             yend = 4.5,
             color = "black") +
    annotate(geom = "segment",
             x = c(5.815, 6.185),
             xend = c(5.815, 6.185),
             y = c(4.5, 4.5),
             yend = c(4.2, 4),
             color = "black") +
    annotate(geom = "text",
             x = 6,
             y = 4.6,
             label = "*",
             size = 5)

############################################
##### Topographical plots ##################
############################################

### Starting with raw voltage topoplots

# load electrode layout for topo plots
# load electrode layout
elec_loc <- read_csv(here("data", "paper_two", "Equidistant Layout.csv"))
elec_loc <- elec_loc %>%
  rename("channel" = `channel name`) %>%
  filter(channel != "CMS", channel != "DRL")

elec_loc$radian_phi <- pi/180 * elec_loc$phi

elec_loc <-
  elec_loc %>%
  mutate(x = theta * cos(radian_phi),
         y = theta * sin(radian_phi))

# merge coordinates with ERP data
eeg_dat <- read_csv(here("data", "paper_two", "created_data", "erp_avr.csv"))
eeg_dat_long <- pivot_longer(eeg_dat, cols = A1:EXG2, names_to = "electrode", values_to = "amplitude")
eeg_dat_loc <- left_join(eeg_dat_long, elec_loc, by = c("electrode" = "channel"))
topo_elec <- c(paste0("A", 1:32), paste0("B", 1:32))

# define lists of time windows for iteration
lower_time_windows <- list(62.5,
                           120,
                           235,
                           70,
                           500)

upper_time_windows <- list(187.5,
                           240,
                           350,
                           1000,
                           1200)

# iterate over time windows and create list of topoplots
topo_plot_lst <- map2(lower_time_windows,
                      upper_time_windows, ~ {
  eeg_dat_loc %>%
  filter(!is.na(amplitude),
         electrode %in% topo_elec,
         between(ms, .x, .y)) %>%
  group_by(electrode, x, y) %>%
  dplyr::summarize(amplitude = mean(amplitude, na.rm = TRUE)) %>%
  topoplot(interp_limit = "head",
           scaling = 0.25) +
  # coord_fixed(clip = 'off') +
  ggtitle(paste(as.character(.x), "-", as.character(.y), "ms")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  size = 10))
})

###  Now the component  topoplots
# read in factor score data with ms variable
temp_dat_fac_scores <- read_csv(here("data", "paper_two", "temp_fac_score_dat.csv")) %>%
  rename("electrode" = "elec") %>%
  full_join(elec_loc, by = c("electrode" = "channel"))

# define vectorS for iteration
## component names
comp_list_topo <- c("RC5", "RC11", "RC12", "RC12", "RC2", "RC3")

## electrodes to highlight
electrode_highlight_lst <- list(c("A29", "B26"),
                                c("A29", "B26"),
                                c("B21", "B28"),
                                c("A29", "B26"),
                                c("A29", "B26"),
                                c("A29", "B26", "A26", "B23","B28", "A30", "B27", "A25", "B22")
                                )

# iterate and create list of topo plots
comp_topo_lst <- map2(comp_list_topo,
                      electrode_highlight_lst,
                      ~ {
temp_dat_fac_scores %>%
  rename("amplitude" = .x) %>%
    topoplot(interp_limit = "head",
             highlights = .y,
             scaling = 0.25) +
    # coord_fixed(clip = 'off') +
    theme(legend.position = "none")
    # guides(fill = guide_colorbar(title = expression(paste("Average ",
    #                                                       mu, "V")),
    #                              title.position = "right",
    #                              barwidth = rel(1),
    #                              barheight = rel(6),
    #                              title.theme = element_text(angle = 270)
    #                              )
    #        )
})

arrow_1 <-
  ggplot() +
  theme_classic() +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  lims(x = c(0, 1), y = c(0, 1)) +
  annotate(geom = "segment",
           x = 1,
           xend = 0.3,
           y = 1,
           yend = 0,
           size = 0.5,
           arrow = arrow(type = "open",
                         length = unit (.1, "inches")))

arrow_2 <-
  ggplot() +
  theme_classic() +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  lims(x = c(0, 1), y = c(0, 1)) +
  annotate(geom = "segment",
           x = 0,
           xend = 0.7,
           y = 1,
           yend = 0,
           size = 0.5,
           arrow = arrow(type = "open",
                         length = unit (.1, "inches")))

arrow_3 <-
ggplot() +
  theme_classic() +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  lims(x = c(0, 1), y = c(0, 1)) +
  annotate(geom = "segment",
           x = 0.5,
           xend = 0.5,
           y = 1,
           yend = 0,
           size = 0.5,
           arrow = arrow(type = "open",
                         length = unit (.1, "inches")))

arrow_4 <-
  ggplot() +
  theme_classic() +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  lims(x = c(0, 1), y = c(0, 1)) +
  annotate(geom = "segment",
           x = 0.5,
           xend = 0.5,
           y = 1,
           yend = 0,
           size = 0.5,
           arrow = arrow(type = "open",
                         length = unit (.1, "inches")))

arrow_5 <-
  ggplot() +
  theme_classic() +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  lims(x = c(0, 1), y = c(0, 1)) +
  annotate(geom = "segment",
           x = 0.5,
           xend = 0.5,
           y = 1,
           yend = 0,
           size = 0.5,
           arrow = arrow(type = "open",
                         length = unit (.1, "inches")))

arrow_6 <-
  ggplot() +
  theme_classic() +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  lims(x = c(0, 1), y = c(0, 1)) +
  annotate(geom = "segment",
           x = 0.5,
           xend = 0.5,
           y = 1,
           yend = 0,
           size = 0.5,
           arrow = arrow(type = "open",
                         length = unit (.1, "inches")))

layout <- "
AAABBB##CC##DDDEEE
FFFGGG##HI##JJJKKK
LLLMMMNNNOOOPPPQQQ
RRRRRRRRRRRRRRRRRR
SSSSSSSSSSSSSSSSSS
TTTTTTTTTTTTTTTTTT
UUUUUUUUUUUUUUUUUU
VVVVVVVVVVVVVVVVVV"

(topo_plot_lst[[1]] +
  topo_plot_lst[[2]] +
  topo_plot_lst[[3]] +
  topo_plot_lst[[4]] +
  topo_plot_lst[[5]] +
  arrow_3 +
  arrow_5 +
  arrow_1 +
  arrow_2 +
  arrow_4 +
  arrow_6 +
  comp_topo_lst[[1]] +
  comp_topo_lst[[2]] +
  comp_topo_lst[[3]] +
  comp_topo_lst[[4]] +
  comp_topo_lst[[5]] +
  comp_topo_lst[[6]] +
  p_1 +
  p_2 +
  p_3 +
  p_4 +
  p_5 +
  plot_layout(design = layout,
              heights = c(1, 0.5, 1, 1.5, 1.5, 1.5, 1.5, 1.5)))

ggsave(here("images", "paper_2", "results_images", "contrast_plot.png"),
       plot = last_plot(),
       height = 11,
       width = 9)



# unusued code for blank topos
# elec_loc <- elec_loc %>%
#   mutate(rc2_color = if_else(electrode %in% rc2_elec, TRUE, FALSE),
#          rc3_color = if_else(electrode %in% rc3_elec, TRUE, FALSE),
#          rc5_color = if_else(electrode %in% rc5_elec, TRUE, FALSE),
#          rc11_color = if_else(electrode %in% rc11_elec, TRUE, FALSE),
#          rc12_color = if_else(electrode %in% rc12_elec, TRUE, FALSE),
#          pos_rc12_color = if_else(electrode %in% pos_rc12_elec, TRUE, FALSE)
#   )
#
# var_to_iter <- c("rc5_color", "rc11_color", "rc12_color", "pos_rc12_color", "rc2_color", "rc3_color")
#
#     ggplot(elec_loc, aes(x = x, y = y)) +
#       geom_mask(size = 6) +
#       geom_head(interp_limit = "head") +
#       geom_channels(aes(color = elec_loc[[.x]])) +
#       scale_color_manual(values = c("black", "red")) +
#       coord_equal() +
#       theme_void() +
#       theme(legend.position = "none")
