#' ---
#' title: "ERP Plots"
#' author: "Ian J. Kahrilas"
#' date: "2020/5/8"
#' output: "html_document"
#' ---
#+ Load packages, include = FALSE
library(tidyverse)
library(eegUtils)
library(haven)
library(here)
library(patchwork)
#'
#' read in data
#+ read in eeg data, incude = FALSE
eeg_df_mast <- read_csv(here("data", "paper_two", "created_data", "erp_mast.csv"))
eeg_df_avr <- read_csv(here("data", "paper_two", "created_data", "erp_avr.csv"))
spn_df <- read_csv(here("data", "paper_two", "created_data", "spn.csv")) %>%
  mutate(block = str_remove(block, "Pre_"))

#'
#' define clusters of electrodes and time windows for each component
#+ electrode clusters and time windows
# clusters
lpp_elec <- "B28" # look at 400 - 800 ms and 800 - 2000 ms windows
epn_elec <- c("A29", "B26") # 250 - 375 ms
N170_elec <- c("A29", "B26") # 160 - 200 ms
spn_elec <- "A5" # 450 - 1250 ms
#'
#' Omit extreme cases as indicated by histograms and QQ plots in Univarite Exploration script
#+ Omit extreme cases
eeg_df_mast$B28[eeg_df_mast$pid == 206201823 & eeg_df_mast$block == "Pos_Inc"] <- NA

#' Create plots for each component with all conditions
#+ plot creation
erp_plot_fun <- function(dat, cluster, comp_name, time_window_low, time_window_high) {
  dat %>%
    select(all_of(cluster),  block:prop_trials) %>%
    filter(ms < 2000) %>%
    pivot_longer(., cols = cluster, names_to = "electrode", values_to = "mv") %>%
    group_by(block, ms) %>%
    summarize(mv = mean(mv, na.rm = TRUE)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = c(time_window_low, time_window_high), linetype = "solid", size = 1.05) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Time (ms)",
         y = expression(paste("Amplitude ( ",mu,"V)")),
         title = paste("Average", comp_name, "Waveforms")) +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16))
}
    # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
    # scale_color_viridis_d(name = "Group",
    #                       breaks = c("Neg_Dec", "Neg_Watch", "Neg_Inc",
    #                                  "Neu_Watch",
    #                                  "Pos_Dec", "Pos_Watch", "Pos_Inc"),
    #                       labels = c("Negative Decrease", "Negative Watch", "Negative Increase",
    #                                  "Neutral Watch",
    #                                  "Positive Decrease", "Positive Watch", "Positive Increase"))

#' plots for just passive conditions
#+ passive plot function
erp_plot_fun_passive <- function(dat, cluster, comp_name, time_window_low, time_window_high) {
  dat %>%
    filter(block %in% c("Neg_Watch", "Neu_Watch", "Pos_Watch")) %>%
    select(all_of(cluster),  block:prop_trials) %>%
    filter(ms < 2000) %>%
    pivot_longer(., cols = cluster, names_to = "electrode", values_to = "mv") %>%
    group_by(block, ms) %>%
    summarize(mv = mean(mv, na.rm = TRUE)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = c(time_window_low, time_window_high), linetype = "solid", size = 1.05) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Time (ms)",
         y = expression(paste("Amplitude ( ",mu,"V)")),
         title = paste("Average", comp_name, "Waveforms")) +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16))
}
    # # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
    # scale_color_viridis_d(name = "Group",
    #                       breaks = c("Neg_Watch", "Neu_Watch", "Pos_Watch"),
    #                       labels = c("Negative Watch", "Neutral Watch", "Positive Watch"))

#' plots for positive conditions
#+ positive plot function
erp_plot_fun_positive <- function(dat, cluster, comp_name, time_window_low, time_window_high) {
  dat %>%
    filter(block %in% c("Pos_Dec", "Pos_Watch", "Pos_Inc")) %>%
    select(all_of(cluster),  block:prop_trials) %>%
    filter(ms < 2000) %>%
    pivot_longer(., cols = cluster, names_to = "electrode", values_to = "mv") %>%
    group_by(block, ms) %>%
    summarize(mv = mean(mv, na.rm = TRUE)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = c(time_window_low, time_window_high), linetype = "solid", size = 1.05) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Time (ms)",
         y = expression(paste("Amplitude ( ",mu,"V)")),
         title = paste("Average", comp_name, "Waveforms")) +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16))
}
    # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
    # scale_color_viridis_d(name = "Group",
    #                       breaks = c("Pos_Dec", "Pos_Watch", "Pos_Inc"),
    #                       labels = c("Positive Decrease", "Positive Watch", "Positive Increase"))

#' plot function for negative condition
#+ negative plot function
erp_plot_fun_negative <- function(dat, cluster, comp_name, time_window_low, time_window_high) {
  dat %>%
    filter(block %in% c("Neg_Dec", "Neg_Watch", "Neg_Inc")) %>%
    select(all_of(cluster),  block:prop_trials) %>%
    filter(ms < 2000) %>%
    pivot_longer(., cols = cluster, names_to = "electrode", values_to = "mv") %>%
    group_by(block, ms) %>%
    summarize(mv = mean(mv, na.rm = TRUE)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = c(time_window_low, time_window_high), linetype = "solid", size = 1.05) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Time (ms)",
         y = expression(paste("Amplitude ( ",mu,"V)")),
         title = paste("Average", comp_name, "Waveforms")) +
    theme_classic() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16))
}
    # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
    # scale_color_viridis_d(name = "Group",
    #                       breaks = c("Neg_Dec", "Neg_Watch", "Neg_Inc"),
    #                       labels = c("Negative Decrease", "Negative Watch", "Negative Increase"))

#'
#' Use pmap to iterate plotting function over list of parameters.
#+ iterate and plot
# all
plots_all <- pmap(list(dat = list(eeg_df_mast,
                                  eeg_df_mast,
                                  eeg_df_avr,
                                  eeg_df_avr,
                                  spn_df),
                   cluster = list(lpp_elec,
                                  lpp_elec,
                                  epn_elec,
                                  N170_elec,
                                  spn_elec),
                   comp_name = c("LPP",
                                 "Late LPP",
                                 "EPN",
                                 "N170",
                                 "SPN"),
                   time_window_low = c(400,
                                       800,
                                       225,
                                       160,
                                       315),
                   time_window_high = c(800,
                                        2000,
                                        375,
                                        200,
                                        1250)),
                   .f = erp_plot_fun)

# passive blocks
plots_passive <- pmap(list(dat = list(eeg_df_mast,
                                      eeg_df_mast,
                                      eeg_df_avr,
                                      eeg_df_avr,
                                      spn_df),
                           cluster = list(lpp_elec,
                                          lpp_elec,
                                          epn_elec,
                                          N170_elec,
                                          spn_elec),
                           comp_name = c("LPP",
                                         "Late LPP",
                                         "EPN",
                                         "N170",
                                         "SPN"),
                           time_window_low = c(450,
                                               800,
                                               225,
                                               160,
                                               450),
                           time_window_high = c(800,
                                                2000,
                                                375,
                                                200,
                                                1250)),
                      .f = erp_plot_fun_passive)
# positive blocks
plots_positive <- pmap(list(dat = list(eeg_df_mast,
                                       eeg_df_mast,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       spn_df),
                            cluster = list(lpp_elec,
                                           lpp_elec,
                                           epn_elec,
                                           N170_elec,
                                           spn_elec),
                            comp_name = c("LPP",
                                          "Late LPP",
                                          "EPN",
                                          "N170",
                                          "SPN"),
                            time_window_low = c(400,
                                                800,
                                                225,
                                                160,
                                                450),
                            time_window_high = c(800,
                                                 2000,
                                                 375,
                                                 200,
                                                 1250)),
                       .f = erp_plot_fun_positive)
# negative blocks
plots_negative <- pmap(list(dat = list(eeg_df_mast,
                                       eeg_df_mast,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       spn_df),
                            cluster = list(lpp_elec,
                                           lpp_elec,
                                           epn_elec,
                                           N170_elec,
                                           spn_elec),
                            comp_name = c("LPP",
                                          "Late LPP",
                                          "EPN",
                                          "N170",
                                          "SPN"),
                            time_window_low = c(400,
                                                800,
                                                225,
                                                160,
                                                450),
                            time_window_high = c(800,
                                                 2000,
                                                 375,
                                                 200,
                                                 1250)),
                       .f = erp_plot_fun_negative)
#'
#' save images to workspace
#+ save the images
# all
map2(plots_all, c("LPP", "Late_LPP", "EPN", "N170", "SPN"), ~{
  ggsave(plot = .x, filename = here("images", "paper_2", "average_waveforms", "all_blocks", paste0(.y, "_all.png")), device = "png", width = 8, height = 5, scale = 1.5)
})
# passive blocks
map2(plots_passive, c("LPP", "Late_LPP", "EPN", "N170", "SPN"), ~{
  ggsave(plot = .x, filename = here("images", "paper_2", "average_waveforms", "passive_blocks", paste0(.y, "_passive.png")), device = "png", width = 8, height = 5, scale = 1.5)
})
# positive blocks
map2(plots_positive, c("LPP", "Late_LPP", "EPN", "N170", "SPN"), ~{
  ggsave(plot = .x, filename = here("images", "paper_2", "average_waveforms", "positive_blocks", paste0(.y, "_positive.png")), device = "png", width = 8, height = 5, scale = 1.5)
})
# negative blocks
map2(plots_negative, c("LPP", "Late_LPP", "EPN", "N170", "SPN"), ~{
  ggsave(plot = .x, filename = here("images", "paper_2", "average_waveforms", "negative_blocks", paste0(.y, "_negative.png")), device = "png", width = 8, height = 5, scale = 1.5)
})


eeg_df_mast %>%
  select(all_of(lpp_elec),  pid:prop_trials) %>%
  filter(ms < 2000) %>%
  pivot_longer(., cols = all_of(lpp_elec), names_to = "electrode", values_to = "mv") %>%
  group_by(block, ms) %>%
  mutate(mean_mv = mean(mv, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(mean_mv, ms), color = "blue") +
  #geom_line(aes(mv, ms, group = pid), alpha = 0.3) +
  facet_wrap(~ block) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(400, 800), linetype = "solid", size = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Time (ms)",
       y = expression(paste("Amplitude (",mu,"V)")),
       title = paste("Average", "LPP", "Waveforms")) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2, "line"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 16))


eeg_df_avr %>%
  filter(ms < 500,
         pid == 20620181) %>%
  select(pid, block, ms, all_of(epn_elec)) %>%
  pivot_longer(., cols = all_of(epn_elec), names_to = "electrode", values_to = "mv") %>%
  group_by(pid, block, ms) %>%
  summarize(mv = mean(mv, na.rm = TRUE)) %>%
  group_by(block, ms) %>%
  mutate(avg_mv = mean(mv, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(ms, mv, group = pid), alpha = 0.3) +
  geom_line(aes(ms, avg_mv), color = "red", size = 1.2) +
  facet_wrap(~ block, ncol = 2, scales = "free_y") +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = 250, xmax = 375, ymin = -Inf, ymax = Inf, alpha = .15) +
  labs(x = "Time (ms)",
       y = expression(paste("Amplitude (",mu,"V)")),
       title = paste("Waveform Variability Plots")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2, "line"),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.position = "none")
