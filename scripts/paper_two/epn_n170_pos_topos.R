#epn n170 pos topos

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

n170_epn_elec <- c("A29", "B26")

# merge coordinates with ERP data
eeg_dat <- read_csv(here("data", "paper_two", "created_data", "erp_avr.csv"))
eeg_dat_long <- pivot_longer(eeg_dat, cols = A1:EXG2, names_to = "electrode", values_to = "amplitude")
eeg_dat_loc <- left_join(eeg_dat_long, elec_loc, by = c("electrode" = "channel"))
topo_elec <- c(paste0("A", 1:32), paste0("B", 1:32))

epn <- eeg_dat_loc %>%
  filter(!is.na(amplitude), block %in% c("Pos_Dec", "Pos_Watch", "Pos_Inc"),
         electrode %in% topo_elec, between(ms, 224, 326)) %>%
  group_by(block, electrode, x, y) %>%
  dplyr::summarize(amplitude = mean(amplitude, na.rm = TRUE)) %>%
  mutate(highlight = if_else(electrode %in% n170_epn_elec, "white", "black"),
         block = if_else(block == "Pos_Dec", "Decrease",
                         if_else(block == "Pos_Watch", "Watch", "Increase")),
         block = factor(block, levels = c("Decrease", "Watch", "Increase")),
         comp = "EPN")

n170_epn_topos <- eeg_dat_loc %>%
  filter(!is.na(amplitude), block %in% c("Pos_Dec", "Pos_Watch", "Pos_Inc"),
         electrode %in% topo_elec, between(ms, 113, 218)) %>%
  group_by(block, electrode, x, y) %>%
  dplyr::summarize(amplitude = mean(amplitude, na.rm = TRUE)) %>%
  mutate(highlight = if_else(electrode %in% n170_epn_elec, "white", "black"),
         block = if_else(block == "Pos_Dec", "Decrease",
                         if_else(block == "Pos_Watch", "Watch", "Increase")),
         block = factor(block, levels = c("Decrease", "Watch", "Increase")),
         comp = "N170") %>%
  bind_rows(., epn) %>%
  mutate(., comp = factor(comp)) %>%
  ggplot(., aes(x = x, y = y, fill = amplitude)) +
  stat_scalpmap() +
  geom_mask(scale_fac = 1.6, size = .01) +
  geom_head() +
  geom_channels(aes(color = highlight), size = 0.125) +
  scale_fill_viridis_c(limits = c(-12, 12), oob = scales::squish) +
  scale_color_manual(breaks = c("black", "white"),
                     values = c("black", "white"),
                     guide = FALSE) +
  labs(fill = expression(paste("Amplitude (", mu,"V)"))) +
  coord_equal() +
  theme_void() +
  facet_grid(vars(comp), vars(block), switch = "y") +
  theme(strip.text.x = element_text(size = 10, vjust = 1),
        strip.text.y = element_text(size = 10, angle = 90))

pos_n170_epn <- eeg_dat %>%
  filter(block %in% c("Pos_Dec", "Pos_Watch", "Pos_Inc")) %>%
  select(all_of(n170_epn_elec),  block:prop_trials) %>%
  filter(ms <= 1000) %>%
  pivot_longer(., cols = n170_epn_elec, names_to = "electrode", values_to = "mv") %>%
  group_by(block, ms) %>%
  dplyr::summarize(mv = mean(mv, na.rm = TRUE)) %>%
  ggplot(., aes(ms, mv, color = block)) +
  geom_line(size = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(113, 218), linetype = "solid", size = 1.05) +
  geom_vline(xintercept = c(224, 326), linetype = "solid", size = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = 113, xmax = 218, ymin = -Inf, ymax = Inf, alpha = .2) +
  annotate("rect", xmin = 224, xmax = 326, ymin = -Inf, ymax = Inf, alpha = .2) +
  labs(x = "Time (ms)",
       y = expression(paste("Amplitude ( ",mu,"V)"))) +
  #      title = paste("Average LPP Waveforms")) +
  theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(2, "line"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 10)) +
  # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
  scale_color_manual(values = c("blue", "green", "red"),
                     name = "Block",
                     breaks = c("Pos_Dec", "Pos_Watch", "Pos_Inc"),
                     labels = c("Positive Decrease", "Positive Watch", "Positive Increase"))

n170_epn_topos / pos_n170_epn + plot_layout(widths = 2)

ggsave("images/paper_2/epn_n170_pos_topos.png")
