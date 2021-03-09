# ERPs for non-PCA data and components

library(tidyverse)
library(eegUtils)
library(haven)
library(here)
library(patchwork)

# read in data
pca_dat <- read_csv(here("data", "paper_two", "temp_fac_score_erp.csv"))
nor_dat <- read_csv(here("data", "paper_two", "created_data", "erp_avr.csv"))

po7_po8_plot <-
  nor_dat %>%
  select(A29, B26,  block:prop_trials) %>%
  filter(ms < 2000) %>%
  pivot_longer(., cols = c(A29, B26), names_to = "electrode", values_to = "mv") %>%
  group_by(ms) %>%
  summarize(mv = mean(mv, na.rm = TRUE)) %>%
  ggplot(., aes(ms, mv)) +
  geom_line(size = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = c(65, 187), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 65, xmax = 187,), fill = "blue", alpha = 0.1) +
  # geom_vline(xintercept = c(125, 240), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 125, xmax = 240), fill = "red", alpha = 0.1) +
  # geom_vline(xintercept = c(242, 356), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 242, xmax = 356), fill = "green", alpha = 0.1) +
  # geom_vline(xintercept = c(65, 1000), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 65, xmax = 1000), fill = "purple", alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Time (ms)",
       y = expression(paste("Amplitude (",mu,"V)")),
       title = paste("Average PO7/PO8 Waveform")) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 14))

pz_cpz_plot <-
  nor_dat %>%
  select(B21, B28,  block:prop_trials) %>%
  filter(ms < 2000) %>%
  pivot_longer(., cols = c(B21, B28), names_to = "electrode", values_to = "mv") %>%
  group_by(ms) %>%
  summarize(mv = mean(mv, na.rm = TRUE)) %>%
  ggplot(., aes(ms, mv)) +
  geom_line(size = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = c(65, 187), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 65, xmax = 187,), fill = "blue", alpha = 0.1) +
  # geom_vline(xintercept = c(125, 240), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 125, xmax = 240), fill = "red", alpha = 0.1) +
  # geom_vline(xintercept = c(242, 356), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 242, xmax = 356), fill = "green", alpha = 0.1) +
  # geom_vline(xintercept = c(65, 1000), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 65, xmax = 1000), fill = "purple", alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Time (ms)",
       y = expression(paste("Amplitude (",mu,"V)")),
       title = paste("Average CPz and Pz Waveform")) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 14))

lpp_plot <-
  nor_dat %>%
  select(A29, B26, A26, B23,B28, A30, B27, A25, B22,  block:prop_trials) %>%
  filter(ms < 2000) %>%
  pivot_longer(., cols = c(A29, B26, A26, B23,B28, A30, B27, A25, B22), names_to = "electrode", values_to = "mv") %>%
  group_by(ms) %>%
  summarize(mv = mean(mv, na.rm = TRUE)) %>%
  ggplot(., aes(ms, mv)) +
  geom_line(size = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = c(65, 187), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 65, xmax = 187,), fill = "blue", alpha = 0.1) +
  # geom_vline(xintercept = c(125, 240), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 125, xmax = 240), fill = "red", alpha = 0.1) +
  # geom_vline(xintercept = c(242, 356), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 242, xmax = 356), fill = "green", alpha = 0.1) +
  # geom_vline(xintercept = c(65, 1000), linetype = "solid", size = 1.05, alpha = 0.5) +
  # geom_ribbon(aes(xmin = 65, xmax = 1000), fill = "purple", alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Time (ms)",
       y = expression(paste("Amplitude (",mu,"V)")),
       title = paste("Average PO3/PO4, PO7/PO8, P1/P2,\nP3/P5, P6/P4, and Pz Waveform")) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 14))

# list of electrode sites
elec_selections <- list(c("A29", "B26"),
                        c("A29", "B26"),
                        c("B21", "B28"),
                        c("A29", "B26"),
                        c("A29", "B26"),
                        c("A29", "B26", "A26", "B23","B28", "A30", "B27", "A25", "B22"))

# list of component sites for map function ordered chronologically
component_list <- list("RC5",
                       "RC11",
                       "RC12",
                       "RC12",
                       "RC2",
                       "RC3")

# name
comp_names <- list("124 ms",
                   "162 ms",
                   "Negative 259 ms",
                   "Positive 259 ms",
                   "381 ms",
                   "740 ms")

comp_erps <-
  pmap(list(component_list, elec_selections, comp_names), ~ {
  pca_dat %>%
  filter(comp == ..1,
         elec %in% ..2) %>%
  group_by(ms) %>%
  summarize(mv = mean(mv, na.rm = TRUE)) %>%
  ggplot(., aes(ms, mv)) +
  geom_line(size = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Time (ms)",
       y = expression(paste("Amplitude (",mu,"V)")),
       title = paste("Average", ..3, "\nComponent Waveform")) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 14))
})

# arrow
arrow <-
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
           size = 1,
           arrow = arrow(type = "closed",
                         length = unit (.1, "inches"))) +
  annotate(geom = "text",
           x = 0.47,
           y = 0.51,
           label = "PCA",
           angle = 90,
           size = 7)

arrow_2 <-
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
           size = 1,
           arrow = arrow(type = "closed",
                         length = unit (.1, "inches"))) +
  annotate(geom = "text",
           x = 0.46,
           y = 0.52,
           label = "PCA",
           angle = 90,
           size = 7)
# compose plots for manuscript
layout <- '
#AA#
BBBB
CCDD
EEFF
'

po7_po8_plot +
  arrow +
  comp_erps[[1]] +
  comp_erps[[2]] +
  comp_erps[[4]] +
  comp_erps[[5]] +
  plot_layout(design = layout,
              heights = c(1.5, 0.8, 1, 1))

ggsave(here("images", "paper_2", "erp_raw_comp", "po7_po8_erps.png"),
       plot = last_plot(),
       height = 8,
       width = 8)

(pz_cpz_plot | lpp_plot) /
  (arrow_2 | arrow_2) /
  (comp_erps[[3]] | comp_erps[[6]]) +
  plot_layout(heights = c(1.5, 1, 1.5))

ggsave(here("images", "paper_2", "erp_raw_comp", "pz_cpz_lpp_erps.png"),
       plot = last_plot(),
       height = 6,
       width = 9)

