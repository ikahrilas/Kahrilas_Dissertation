# LPP positive topos

library(devtools)
library(eegUtils) # remotes::install_github("craddm/eegUtils")
library(tidyverse)
library(here)
library(patchwork)

# load electrode layout
elec_loc <- read_csv(here("data", "paper_two", "Equidistant Layout.csv"))
elec_loc <- elec_loc %>%
  rename("channel" = `channel name`) %>%
  filter(channel != "CMS", channel != "DRL")

elec_loc$radian_phi <- pi/180 * elec_loc$phi

elec_loc <- elec_loc %>%
  mutate(x = theta * cos(radian_phi),
         y = theta * sin(radian_phi))

Pz_elec <- c("B28")

# merge coordinates with ERP data
eeg_dat <- read_csv(here("data", "paper_two", "created_data", "erp_mast.csv"))
eeg_dat$B28[eeg_dat$pid == 206201823 & eeg_dat$block == "Pos_Inc"] <- NA
eeg_dat$A5[eeg_dat$pid == 206201843 & eeg_dat$block == "Neg_Watch"] <- NA
eeg_dat_long <- pivot_longer(eeg_dat, cols = A1:EXG2, names_to = "electrode", values_to = "amplitude")
eeg_dat_loc <- left_join(eeg_dat_long, elec_loc, by = c("electrode" = "channel"))
topo_elec <- c(paste0("A", 1:32), paste0("B", 1:32))

late_tmp <- eeg_dat_loc %>%
  filter(!is.na(amplitude), block %in% c("Pos_Dec", "Pos_Watch", "Pos_Inc"),
         electrode %in% topo_elec, between(ms, 800, 2000)) %>%
  group_by(block, electrode, x, y) %>%
  dplyr::summarize(amplitude = mean(amplitude, na.rm = TRUE)) %>%
  mutate(highlight = if_else(electrode %in% Pz_elec, "white", "black"),
         block = if_else(block == "Pos_Dec", "Decrease",
                         if_else(block == "Pos_Watch", "Watch", "Increase")),
         block = factor(block, levels = c("Decrease", "Watch", "Increase")),
         comp = "Late")

reg_topos <- eeg_dat_loc %>%
  filter(!is.na(amplitude), block %in% c("Pos_Dec", "Pos_Watch", "Pos_Inc",
                                         "Neg_Dec", "Neg_Watch", "Neg_Inc"),
         electrode %in% topo_elec, between(ms, 400, 800)) %>%
  group_by(block, electrode, x, y) %>%
  dplyr::summarize(amplitude = mean(amplitude, na.rm = TRUE)) %>%
  mutate(highlight = if_else(electrode %in% Pz_elec, "white", "black")) %>%
  separate(block, c("valence_cond", "regulation_cond"), "_") %>%
  mutate(valence_condition = if_else(valence_cond == "Neg", "Negative",
                                     if_else(valence_cond == "Pos", "Positive", "Neutral")),
         regulation_condition = if_else(regulation_cond == "Dec", "Decrease",
                                        if_else(regulation_cond == "Inc", "Increase", "Watch"))) %>%
  ggplot(., aes(x = x, y = y, fill = amplitude)) +
  stat_scalpmap() +
  geom_mask(scale_fac = 1.6, size = .01) +
  geom_head() +
  geom_channels(aes(color = highlight), size = 0.125) +
  scale_fill_viridis_c(limits = c(-8, 8), oob = scales::squish) +
  scale_color_manual(breaks = c("black", "white"),
                     values = c("black", "white"),
                     guide = FALSE) +
  labs(fill = expression(paste("Amplitude (", mu,"V)"))) +
  coord_equal() +
  theme_void() +
  facet_grid(vars(valence_condition), vars(regulation_condition), switch = "y") +
  theme(strip.text.x = element_text(size = 12, vjust = 1),
        strip.text.y = element_text(size = 12, angle = 90, vjust = 1))

reg_lpp <- eeg_dat %>%
  filter(block %in% c("Pos_Dec", "Pos_Watch", "Pos_Inc",
                      "Neg_Dec", "Neg_Watch", "Neg_Inc")) %>%
  select(all_of(Pz_elec),  block:prop_trials) %>%
  filter(ms <= 1500) %>%
  pivot_longer(., cols = Pz_elec, names_to = "electrode", values_to = "mv") %>%
  group_by(block, ms) %>%
  dplyr::summarize(mv = mean(mv, na.rm = TRUE)) %>%
  separate(block, c("valence_cond", "regulation_cond"), "_") %>%
  mutate(valence_condition = if_else(valence_cond == "Neg", "Negative",
                                     if_else(valence_cond == "Pos", "Positive", "Neutral")),
         regulation_condition = if_else(regulation_cond == "Dec", "Decrease",
                                        if_else(regulation_cond == "Inc", "Increase", "Watch"))) %>%
  ggplot(., aes(ms, mv, color = regulation_condition)) +
  geom_line(size = 1.1) +
  facet_wrap(~ valence_condition, ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(400, 800), linetype = "solid", size = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = 400, xmax = 800, ymin = -Inf, ymax = Inf, alpha = .2) +
  labs(x = "Time (ms)",
       y = expression(paste("Amplitude (",mu,"V)"))) +
  #      title = paste("Average LPP Waveforms")) +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(2, "line"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, hjust = 0)) +
  # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
  scale_color_manual(values = c("red", "green", "blue"),
                     name = "Regulation",
                     breaks = c("Increase", "Watch", "Decrease"))

reg_topos / reg_lpp + plot_layout(widths = c(2))

ggsave("images/paper_2/lpp_reg_topos.png")
