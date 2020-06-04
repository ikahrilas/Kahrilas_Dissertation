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
eeg_df_mast <- read_csv(here("data", "paper_two", "data_from_PER_R_project", "created_data", "erp_mast_no_lp.csv"))
eeg_df_avr <- read_csv(here("data", "paper_two", "data_from_PER_R_project", "created_data", "erp_avr_no_lp.csv"))
spn_df <- read_csv(here("data", "paper_two", "created_data", "spn.csv"))
#eeg_df_mast <- eeg_df_mast %>% filter(prop_trials > .50)
#eeg_df_avr <- eeg_df_avr %>% filter(prop_trials > .50)
#'
#' define clusters of electrodes and time windows for each component
#+ electrode clusters and time windows
# clusters
lpp_elec <- "B28" # look at 450 - 800 ms windows
epn_elec_right <- "B26" # 250 - 375 ms
epn_elec_left <- "A29" # 250 - 375 ms
epn_elec_avg <- c("A29", "B27")
N170_right <- "B26" # 160 - 200 ms
N170_left <- "A29" # 160 - 200 ms
N170_avg <- c("A29", "B27")
spn_elec <- "A5"
#'
#' Omit extreme cases as indicated by histograms and QQ plots in Univarite Exploration script
#+ Omit exteme cases
eeg_df_mast$A25[eeg_df_mast$pid %in% c(206201843, 206201831) & eeg_df_mast$block %in% c("Neg_Watch", "Pos_Watch") & between(eeg_df_mast$ms, 450, 800)] <- NA
eeg_df_mast$B28[eeg_df_mast$pid %in% c(206201843, 206201831) & eeg_df_mast$block %in% c("Neg_Watch", "Pos_Watch") & between(eeg_df_mast$ms, 450, 800)] <- NA
eeg_df_mast$B22[eeg_df_mast$pid %in% c(206201843, 206201831) & eeg_df_mast$block %in% c("Neg_Watch", "Pos_Watch") & between(eeg_df_mast$ms, 450, 800)] <- NA

#' Create plots for each component with all conditions
#+ plot creation
erp_plot_fun <- function(dat, cluster, comp_name, time_window_low, time_window_high) {
  dat %>%
    select(all_of(cluster),  block:prop_trials) %>%
    filter(ms < 1000) %>%
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
          title = element_text(size = 16)) +
    # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
    scale_color_viridis_d(name = "Group",
                          breaks = c("Neg_Dec", "Neg_Watch", "Neg_Inc",
                                     "Neu_Watch",
                                     "Pos_Dec", "Pos_Watch", "Pos_Inc"),
                          labels = c("Negative Decrease", "Negative Watch", "Negative Increase",
                                     "Neutral Watch",
                                     "Positive Decrease", "Positive Watch", "Positive Increase"))
}
#' plots for just passive conditions
#+ passive plot function
erp_plot_fun_passive <- function(dat, cluster, comp_name, time_window_low, time_window_high) {
  dat %>%
    filter(block %in% c("Neg_Watch", "Neu_Watch", "Pos_Watch")) %>%
    select(all_of(cluster),  block:prop_trials) %>%
    filter(ms < 1000) %>%
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
          title = element_text(size = 16)) +
    # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
    scale_color_viridis_d(name = "Group",
                          breaks = c("Neg_Watch", "Neu_Watch", "Pos_Watch"),
                          labels = c("Negative Watch", "Neutral Watch", "Positive Watch"))
}
#' plots for positive conditions
#+ positive plot function
erp_plot_fun_positive <- function(dat, cluster, comp_name, time_window_low, time_window_high) {
  dat %>%
    filter(block %in% c("Pos_Dec", "Pos_Watch", "Pos_Inc")) %>%
    select(all_of(cluster),  block:prop_trials) %>%
    filter(ms < 1000) %>%
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
          title = element_text(size = 16)) +
    # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
    scale_color_viridis_d(name = "Group",
                          breaks = c("Pos_Dec", "Pos_Watch", "Pos_Inc"),
                          labels = c("Positive Decrease", "Positive Watch", "Positive Increase"))
}
#' plot function for negative condition
#+ negative plot function
erp_plot_fun_negative <- function(dat, cluster, comp_name, time_window_low, time_window_high) {
  dat %>%
    filter(block %in% c("Neg_Dec", "Neg_Watch", "Neg_Inc")) %>%
    select(all_of(cluster),  block:prop_trials) %>%
    filter(ms < 1000) %>%
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
          title = element_text(size = 16)) +
    # scale_color_manual(name = "Group", values = c("green", "blue", "red"))
    scale_color_viridis_d(name = "Group",
                          breaks = c("Neg_Dec", "Neg_Watch", "Neg_Inc"),
                          labels = c("Negative Decrease", "Negative Watch", "Negative Increase"))
}
#'
#' Use pmap to iterate plotting function over list of parameters.
#+ iterate and plot
# all
plots_all <- pmap(list(dat = list(eeg_df_mast,
                              eeg_df_avr,
                              eeg_df_avr,
                              eeg_df_avr,
                              eeg_df_avr,
                              eeg_df_avr,
                              eeg_df_avr),
                   cluster = list(lpp_elec,
                               epn_elec_right,
                               epn_elec_left,
                               epn_elec_avg,
                               N170_left,
                               N170_right,
                               N170_avg),
                   comp_name = c("LPP",
                                 "Right EPN",
                                 "Left EPN",
                                 "Bilateral EPN",
                                 "Left N170",
                                 "Right N170",
                                 "Bilateral N170"),
                   time_window_low = c(450,
                                       225,
                                       225,
                                       225,
                                       150,
                                       150,
                                       150),
                   time_window_high = c(800,
                                        375,
                                        375,
                                        375,
                                        190,
                                        190,
                                        190)),
                   .f = erp_plot_fun)
# passive blocks
plots_passive <- pmap(list(dat = list(eeg_df_mast,
                                      eeg_df_avr,
                                      eeg_df_avr,
                                      eeg_df_avr,
                                      eeg_df_avr,
                                      eeg_df_avr,
                                      eeg_df_avr),
                           cluster = list(lpp_elec,
                                          epn_elec_right,
                                          epn_elec_left,
                                          epn_elec_avg,
                                          N170_left,
                                          N170_right,
                                          N170_avg),
                           comp_name = c("LPP",
                                         "Right EPN",
                                         "Left EPN",
                                         "Bilateral EPN",
                                         "Left N170",
                                         "Right N170",
                                         "Bilateral N170"),
                           time_window_low = c(450,
                                               225,
                                               225,
                                               225,
                                               150,
                                               150,
                                               150),
                           time_window_high = c(800,
                                                375,
                                                375,
                                                375,
                                                190,
                                                190,
                                                190)),
                      .f = erp_plot_fun_passive)
# positive blocks
plots_positive <- pmap(list(dat = list(eeg_df_mast,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       eeg_df_avr),
                            cluster = list(lpp_elec,
                                           epn_elec_right,
                                           epn_elec_left,
                                           epn_elec_avg,
                                           N170_left,
                                           N170_right,
                                           N170_avg),
                            comp_name = c("LPP",
                                          "Right EPN",
                                          "Left EPN",
                                          "Bilateral EPN",
                                          "Left N170",
                                          "Right N170",
                                          "Bilateral N170"),
                            time_window_low = c(450,
                                                225,
                                                225,
                                                225,
                                                150,
                                                150,
                                                150),
                            time_window_high = c(800,
                                                 375,
                                                 375,
                                                 375,
                                                 190,
                                                 190,
                                                 190)),
                       .f = erp_plot_fun_positive)
# negative blocks
plots_negative <- pmap(list(dat = list(eeg_df_mast,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       eeg_df_avr,
                                       eeg_df_avr),
                            cluster = list(lpp_elec,
                                           epn_elec_right,
                                           epn_elec_left,
                                           epn_elec_avg,
                                           N170_left,
                                           N170_right,
                                           N170_avg),
                            comp_name = c("LPP",
                                          "Right EPN",
                                          "Left EPN",
                                          "Bilateral EPN",
                                          "Left N170",
                                          "Right N170",
                                          "Bilateral N170"),
                            time_window_low = c(450,
                                                225,
                                                225,
                                                225,
                                                150,
                                                150,
                                                150),
                            time_window_high = c(800,
                                                 375,
                                                 375,
                                                 375,
                                                 190,
                                                 190,
                                                 190)),
                       .f = erp_plot_fun_negative)
#'
#' save images to workspace
#+ save the images
# all
map2(plots_all, c("LPP", "EPN_right", "EPN_left", "EPN_bilateral", "N170_left", "N170_right", "N170_bilateral"), ~{
  ggsave(plot = .x, filename = here("Images", "average_waveforms", "all_blocks", paste0(.y, "_all.png")), device = "png", width = 8, height = 5, scale = 1.5)
})
# passive blocks
map2(plots_passive, c("LPP", "EPN_right", "EPN_left", "EPN_bilateral", "N170_left", "N170_right", "N170_bilateral"), ~{
  ggsave(plot = .x, filename = here("Images", "average_waveforms", "passive_blocks", paste0(.y, "_passive.png")), device = "png", width = 8, height = 5, scale = 1.5)
})
# positive blocks
map2(plots_positive, c("LPP", "EPN_right", "EPN_left", "EPN_bilateral", "N170_left", "N170_right", "N170_bilateral"), ~{
  ggsave(plot = .x, filename = here("Images", "average_waveforms", "positive_blocks", paste0(.y, "_positive.png")), device = "png", width = 8, height = 5, scale = 1.5)
})
# negative blocks
map2(plots_negative, c("LPP", "EPN_right", "EPN_left", "EPN_bilateral", "N170_left", "N170_right", "N170_bilateral"), ~{
  ggsave(plot = .x, filename = here("Images", "average_waveforms", "negative_blocks", paste0(.y, "_negative.png")), device = "png", width = 8, height = 5, scale = 1.5)
})
