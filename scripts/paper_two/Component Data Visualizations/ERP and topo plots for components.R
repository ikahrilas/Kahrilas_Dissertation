#' ---
#' title: "ERP Component Plots"
#' author: "Ian J. Kahrilas"
#' date: "2021/3/2"
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
dat <- read_csv(here("data", "paper_two", "temp_fac_score_erp.csv"))

# derive ERP plots in mV units for each temporal component
## define list of electrodes of interest for each component based on visual
## inspection of topo plots

###################################################
# temporal components that change across conditions
## RC2 @ c(A29, B26)
## RC3 @ c(A29, B26, A26, B23, B28)
## RC5 @ c(A29, B26)
## RC11 @ c(A29, B26)
## RC12 @ c(B21, B28) for negativity
## RC12 @ c(A29, A31, B30, B26) for positivity

# order electrode sites chronologically
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

# titles for each of the plots
titles <- list("124 ms Component Waveforms",
               "162 ms Component Waveforms",
               "Negative 259 ms Component Waveforms",
               "Positive 259 ms Component Waveforms",
               "381 ms Component Waveforms",
               "740 ms Component Waveforms")


# change block variable to factor, reorder, and rename
dat$block <- factor(dat$block, levels = c("Neg_Inc",
                                          "Neg_Watch",
                                          "Neg_Dec",
                                          "Pos_Inc",
                                          "Pos_Watch",
                                          "Pos_Dec",
                                          "Neu_Watch"))

levels(dat$block) <- c("Negative Increase",
                       "Negative Watch",
                       "Negative Decrease",
                       "Positive Increase",
                       "Positive Watch",
                       "Positive Decrease",
                       "Neutral Watch")

# iterate over each component and electrode selection to create ERP plots
watch_plots <-
  pmap(list(component_list,
          elec_selections,
          titles), ~ {
  # watch plots
  dat %>%
       filter(elec %in% ..2,
              comp == ..1,
              block %in% c("Positive Watch", "Neutral Watch", "Negative Watch")) %>%
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
             plot.title = element_text(hjust = 0.5),
             title = element_text(size = 16)) +
       scale_color_manual(breaks = c("Negative Watch", "Positive Watch", "Neutral Watch"),
                          values=c("red", "blue", "gray"))
})

pos_plots <-
  pmap(list(component_list,
            elec_selections,
            titles), ~ {
  # positive conditions
  dat %>%
    filter(elec %in% ..2,
           comp == ..1,
           block %in% c("Positive Watch", "Positive Increase", "Positive Decrease")) %>%
    group_by(block, ms) %>%
    mutate(ms = as.numeric(ms)) %>%
    summarise(mv = mean(mv)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
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
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16)) +
    scale_color_manual(breaks = c("Positive Watch", "Positive Increase", "Positive Decrease"),
                       values=c("blue", "springgreen", "cadetblue1"))
})

neg_plots <-
  pmap(list(component_list,
            elec_selections,
            titles), ~ {
  # negative conditions
  dat %>%
    filter(elec %in% ..2,
           comp == ..1,
           block %in% c("Negative Watch", "Negative Increase", "Negative Decrease")) %>%
    group_by(block, ms) %>%
    mutate(ms = as.numeric(ms)) %>%
    summarise(mv = mean(mv)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
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
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16)) +
    scale_color_manual(breaks = c("Negative Watch", "Negative Increase", "Negative Decrease"),
                       values=c("red", "magenta", "coral"))
})

reg_plots <-
  pmap(list(component_list,
            elec_selections,
            titles), ~ {
  # Inc and Dec conditions
  dat %>%
    filter(elec %in% ..2,
           comp == ..1,
           block %in% c("Negative Increase", "Negative Decrease", "Positive Decrease", "Positive Increase")) %>%
    group_by(block, ms) %>%
    mutate(ms = as.numeric(ms)) %>%
    summarise(mv = mean(mv)) %>%
    ggplot(., aes(ms, mv, color = block)) +
    geom_line(size = 1.1) +
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
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 16)) +
    scale_color_manual(breaks = c("Negative Increase", "Negative Decrease", "Positive Increase", "Positive Decrease"),
                       values=c("magenta", "coral", "springgreen", "cadetblue1"))
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

temp_dat_fac_scores <- read_csv(here("data", "paper_two", "temp_fac_score_dat.csv")) %>%
  rename("electrode" = "elec") %>%
  full_join(elec_loc, by = c("electrode" = "channel")) %>%
  filter(electrode != "EXG1",
         electrode != "EXG2") # not used in plotting, so omit

# create variables representing contrasts among conditions
temp_dat_fac_scores <-
  temp_dat_fac_scores %>%
  pivot_wider(names_from = block,
              values_from = c(n_trials, prop_trials, RC2, RC3, RC5, RC11, RC12)) %>%
  mutate(RC5_neg_watch_neu_watch = RC5_Neg_Watch - RC5_Neu_Watch,
         RC5_pos_watch_neu_watch = RC5_Pos_Watch - RC5_Neu_Watch,
         RC5_neg_watch_pos_watch = RC5_Neg_Watch - RC5_Pos_Watch,
         RC5_neg_dec_neg_watch = RC5_Neg_Dec - RC5_Neg_Watch,
         RC5_neg_inc_neg_watch = RC5_Neg_Inc - RC5_Neg_Watch,
         RC5_neg_inc_neg_dec = RC5_Neg_Inc - RC5_Neg_Dec,
         RC5_pos_dec_pos_watch = RC5_Pos_Dec - RC5_Pos_Watch,
         RC5_pos_inc_pos_watch = RC5_Pos_Inc - RC5_Pos_Watch,
         RC5_pos_inc_pos_dec = RC5_Pos_Inc - RC5_Pos_Dec,
         RC5_neg_inc_pos_inc = RC5_Neg_Inc - RC5_Pos_Inc,
         RC5_neg_dec_pos_dec = RC5_Neg_Dec - RC5_Pos_Dec,
         RC11_neg_watch_neu_watch = RC11_Neg_Watch - RC11_Neu_Watch,
         RC11_pos_watch_neu_watch = RC11_Pos_Watch - RC11_Neu_Watch,
         RC11_neg_watch_pos_watch = RC11_Neg_Watch - RC11_Pos_Watch,
         RC11_neg_dec_neg_watch = RC11_Neg_Dec - RC11_Neg_Watch,
         RC11_neg_inc_neg_watch = RC11_Neg_Inc - RC11_Neg_Watch,
         RC11_neg_inc_neg_dec = RC11_Neg_Inc - RC11_Neg_Dec,
         RC11_pos_dec_pos_watch = RC11_Pos_Dec - RC11_Pos_Watch,
         RC11_pos_inc_pos_watch = RC11_Pos_Inc - RC11_Pos_Watch,
         RC11_pos_inc_pos_dec = RC11_Pos_Inc - RC11_Pos_Dec,
         RC11_neg_inc_pos_inc = RC11_Neg_Inc - RC11_Pos_Inc,
         RC11_neg_dec_pos_dec = RC11_Neg_Dec - RC11_Pos_Dec,
         RC12_neg_watch_neu_watch = RC12_Neg_Watch - RC12_Neu_Watch,
         RC12_pos_watch_neu_watch = RC12_Pos_Watch - RC12_Neu_Watch,
         RC12_neg_watch_pos_watch = RC12_Neg_Watch - RC12_Pos_Watch,
         RC12_neg_dec_neg_watch = RC12_Neg_Dec - RC12_Neg_Watch,
         RC12_neg_inc_neg_watch = RC12_Neg_Inc - RC12_Neg_Watch,
         RC12_neg_inc_neg_dec = RC12_Neg_Inc - RC12_Neg_Dec,
         RC12_pos_dec_pos_watch = RC12_Pos_Dec - RC12_Pos_Watch,
         RC12_pos_inc_pos_watch = RC12_Pos_Inc - RC12_Pos_Watch,
         RC12_pos_inc_pos_dec = RC12_Pos_Inc - RC12_Pos_Dec,
         RC12_neg_inc_pos_inc = RC12_Neg_Inc - RC12_Pos_Inc,
         RC12_neg_dec_pos_dec = RC12_Neg_Dec - RC12_Pos_Dec,
         RC2_neg_watch_neu_watch = RC2_Neg_Watch - RC2_Neu_Watch,
         RC2_pos_watch_neu_watch = RC2_Pos_Watch - RC2_Neu_Watch,
         RC2_neg_watch_pos_watch = RC2_Neg_Watch - RC2_Pos_Watch,
         RC2_neg_dec_neg_watch = RC2_Neg_Dec - RC2_Neg_Watch,
         RC2_neg_inc_neg_watch = RC2_Neg_Inc - RC2_Neg_Watch,
         RC2_neg_inc_neg_dec = RC2_Neg_Inc - RC2_Neg_Dec,
         RC2_pos_dec_pos_watch = RC2_Pos_Dec - RC2_Pos_Watch,
         RC2_pos_inc_pos_watch = RC2_Pos_Inc - RC2_Pos_Watch,
         RC2_pos_inc_pos_dec = RC2_Pos_Inc - RC2_Pos_Dec,
         RC2_neg_inc_pos_inc = RC2_Neg_Inc - RC2_Pos_Inc,
         RC2_neg_dec_pos_dec = RC2_Neg_Dec - RC2_Pos_Dec,
         RC3_neg_watch_neu_watch = RC3_Neg_Watch - RC3_Neu_Watch,
         RC3_pos_watch_neu_watch = RC3_Pos_Watch - RC3_Neu_Watch,
         RC3_neg_watch_pos_watch = RC3_Neg_Watch - RC3_Pos_Watch,
         RC3_neg_dec_neg_watch = RC3_Neg_Dec - RC3_Neg_Watch,
         RC3_neg_inc_neg_watch = RC3_Neg_Inc - RC3_Neg_Watch,
         RC3_neg_inc_neg_dec = RC3_Neg_Inc - RC3_Neg_Dec,
         RC3_pos_dec_pos_watch = RC3_Pos_Dec - RC3_Pos_Watch,
         RC3_pos_inc_pos_watch = RC3_Pos_Inc - RC3_Pos_Watch,
         RC3_pos_inc_pos_dec = RC3_Pos_Inc - RC3_Pos_Dec,
         RC3_neg_inc_pos_inc = RC3_Neg_Inc - RC3_Pos_Inc,
         RC3_neg_dec_pos_dec = RC3_Neg_Dec - RC3_Pos_Dec
         )

# define vectors for iteration
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

watch_contrasts <- list(c("RC5_neg_watch_neu_watch",
                        "RC5_pos_watch_neu_watch",
                        "RC5_neg_watch_pos_watch"),
                        c("RC11_neg_watch_neu_watch",
                        "RC11_pos_watch_neu_watch",
                        "RC11_neg_watch_pos_watch"),
                        c("RC12_neg_watch_neu_watch",
                        "RC12_pos_watch_neu_watch",
                        "RC12_neg_watch_pos_watch"),
                        c("RC12_neg_watch_neu_watch",
                        "RC12_pos_watch_neu_watch",
                        "RC12_neg_watch_pos_watch"),
                        c("RC2_neg_watch_neu_watch",
                        "RC2_pos_watch_neu_watch",
                        "RC2_neg_watch_pos_watch"),
                        c("RC3_neg_watch_neu_watch",
                        "RC3_pos_watch_neu_watch",
                        "RC3_neg_watch_pos_watch"))

pos_contrasts <- list(c("RC5_pos_inc_pos_watch",
                        "RC5_pos_dec_pos_watch",
                        "RC5_pos_inc_pos_dec"),
                      c("RC11_pos_inc_pos_watch",
                        "RC11_pos_dec_pos_watch",
                        "RC11_pos_inc_pos_dec"),
                      c("RC12_pos_inc_pos_watch",
                        "RC12_pos_dec_pos_watch",
                        "RC12_pos_inc_pos_dec"),
                      c("RC12_pos_inc_pos_watch",
                        "RC12_pos_dec_pos_watch",
                        "RC12_pos_inc_pos_dec"),
                      c("RC2_pos_inc_pos_watch",
                        "RC2_pos_dec_pos_watch",
                        "RC2_pos_inc_pos_dec"),
                      c("RC3_pos_inc_pos_watch",
                        "RC3_pos_dec_pos_watch",
                        "RC3_pos_inc_pos_dec"))

neg_contrasts <- list(c("RC5_neg_inc_neg_watch",
                        "RC5_neg_dec_neg_watch",
                        "RC5_neg_inc_neg_dec"),
                      c("RC11_neg_inc_neg_watch",
                        "RC11_neg_dec_neg_watch",
                        "RC11_neg_inc_neg_dec"),
                      c("RC12_neg_inc_neg_watch",
                        "RC12_neg_dec_neg_watch",
                        "RC12_neg_inc_neg_dec"),
                      c("RC12_neg_inc_neg_watch",
                        "RC12_neg_dec_neg_watch",
                        "RC12_neg_inc_neg_dec"),
                      c("RC2_neg_inc_neg_watch",
                        "RC2_neg_dec_neg_watch",
                        "RC2_neg_inc_neg_dec"),
                      c("RC3_neg_inc_neg_watch",
                        "RC3_neg_dec_neg_watch",
                        "RC3_neg_inc_neg_dec"))

reg_contrasts <- list(c("RC5_neg_inc_pos_inc",
                        "RC5_neg_dec_pos_dec"),
                      c("RC11_neg_inc_pos_inc",
                        "RC11_neg_dec_pos_dec"),
                      c("RC12_neg_inc_pos_inc",
                        "RC12_neg_dec_pos_dec"),
                      c("RC12_neg_inc_pos_inc",
                        "RC12_neg_dec_pos_dec"),
                      c("RC2_neg_inc_pos_inc",
                        "RC2_neg_dec_pos_dec"),
                      c("RC3_neg_inc_pos_inc",
                        "RC3_neg_dec_pos_dec"))

highlights <- list(c("A29", "B26"),
                   c("A29", "B26"),
                   c("B21", "B28"),
                   c("A29", "B26"),
                   c("A29", "B26"),
                   c("A29", "B26", "A26", "B23","B28", "A30", "B27", "A25", "B22"))

# define functions for topo plots
plot_fun_watch <- function(contrasts, highlights, erp){
  values <-
    temp_dat_fac_scores %>%
    select(contrasts, electrode) %>%
    filter(electrode %in% highlights) %>%
    group_by(electrode) %>%
    summarise(across(starts_with("RC"), ~mean(.x, na.rm = TRUE))) %>%
    select(-electrode) %>%
    pivot_longer(cols = everything()) %>%
    select(value) %>%
    pull()

  maximum <- max(values)
  minimum <- min(values)

  neg_watch_neu_watch <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[1]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Neg Watch - Neu Watch") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  pos_watch_neu_watch <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[2]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Pos Watch - Neu Watch") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  neg_watch_pos_watch <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[3]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Neg Watch - Pos Watch") +
    theme(plot.title = element_text(hjust = 0.5))

list(neg_watch_neu_watch, pos_watch_neu_watch, neg_watch_pos_watch)

  # (neg_watch_neu_watch | plot_spacer() | pos_watch_neu_watch | plot_spacer() | neg_watch_pos_watch) / erp +
  #   plot_layout(heights = c(1, 1.4))
}

plot_fun_pos <- function(contrasts, highlights, erp){
  values <-
    temp_dat_fac_scores %>%
    select(contrasts, electrode) %>%
    filter(electrode %in% highlights) %>%
    group_by(electrode) %>%
    summarise(across(starts_with("RC"), ~mean(.x, na.rm = TRUE))) %>%
    select(-electrode) %>%
    pivot_longer(cols = everything()) %>%
    select(value) %>%
    pull()

  maximum <- max(values)
  minimum <- min(values)

  pos_inc_pos_watch <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[1]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Pos Inc - Pos Watch") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  pos_dec_pos_watch <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[2]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Pos Dec - Pos Watch") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  pos_inc_pos_dec <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[3]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Pos Inc - Pos Dec") +
    theme(plot.title = element_text(hjust = 0.5))

  list(pos_inc_pos_watch, pos_dec_pos_watch, pos_inc_pos_dec)
  # (pos_inc_pos_watch | plot_spacer() | pos_dec_pos_watch | plot_spacer() | pos_inc_pos_dec) / erp +
  #   plot_layout(heights = c(1, 1.4))
}

plot_fun_neg <- function(contrasts, highlights, erp){
  values <-
    temp_dat_fac_scores %>%
    select(contrasts, electrode) %>%
    filter(electrode %in% highlights) %>%
    group_by(electrode) %>%
    summarise(across(starts_with("RC"), ~mean(.x, na.rm = TRUE))) %>%
    select(-electrode) %>%
    pivot_longer(cols = everything()) %>%
    select(value) %>%
    pull()

  maximum <- max(values)
  minimum <- min(values)

  neg_inc_neg_watch <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[1]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Neg Inc - Neg Watch") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  neg_dec_neg_watch <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[2]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Neg Dec - Neg Watch") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  neg_inc_neg_dec <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[3]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Neg Inc - Neg Dec") +
    theme(plot.title = element_text(hjust = 0.5))

  list(neg_inc_neg_watch, neg_dec_neg_watch, neg_inc_neg_dec)
  # (neg_inc_neg_watch | plot_spacer() | neg_dec_neg_watch | plot_spacer() | neg_inc_neg_dec) / erp +
  #   plot_layout(heights = c(1, 1.4))
}

plot_fun_reg <- function(contrasts, highlights, erp){
  values <-
    temp_dat_fac_scores %>%
    select(contrasts, electrode) %>%
    filter(electrode %in% highlights) %>%
    group_by(electrode) %>%
    summarise(across(starts_with("RC"), ~mean(.x, na.rm = TRUE))) %>%
    select(-electrode) %>%
    pivot_longer(cols = everything()) %>%
    select(value) %>%
    pull()

  maximum <- max(values)
  minimum <- min(values)

  neg_inc_pos_inc <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[1]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Neg Inc - Pos Inc") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")

  neg_dec_pos_dec <-
    temp_dat_fac_scores %>%
    rename("amplitude" = contrasts[[2]]) %>%
    topoplot(interp_limit = "head",
             highlights = highlights,
             limits = c(minimum - .3, maximum + .3),
             scaling = 0.8) +
    ggtitle("Neg Dec - Pos Dec") +
    theme(plot.title = element_text(hjust = 0.5))

list(neg_inc_pos_inc, neg_dec_pos_dec)
}

# iterate and create topo plots
all_watch_topos <- pmap(list(watch_contrasts, highlights, watch_plots), ~ plot_fun_watch(..1, ..2, ..3))

all_pos_topos <- pmap(list(pos_contrasts, highlights, pos_plots), ~ plot_fun_pos(..1, ..2, ..3))

all_neg_topos <- pmap(list(neg_contrasts, highlights, neg_plots), ~ plot_fun_neg(..1, ..2, ..3))

all_reg_topos <- pmap(list(reg_contrasts, highlights, reg_plots), ~ plot_fun_reg(..1, ..2, ..3))

# compose the final images

layout <- '
AAABBBCCC
DDDDDDDDD
EEEFFFGGG
HHHHHHHHH
IIIJJJKKK
LLLLLLLLL
MMM###NNN
OOOOOOOOO
'

rc5_component_plots <-
  all_watch_topos[[1]][[1]] +
    all_watch_topos[[1]][[2]] +
    all_watch_topos[[1]][[3]] +
    watch_plots[[1]] +
    all_pos_topos[[1]][[1]] +
    all_pos_topos[[1]][[2]] +
    all_pos_topos[[1]][[3]] +
    pos_plots[[1]] +
    all_neg_topos[[1]][[1]] +
    all_neg_topos[[1]][[2]] +
    all_neg_topos[[1]][[3]] +
    neg_plots[[1]] +
    all_reg_topos[[1]][[1]] +
    all_reg_topos[[1]][[2]] +
    reg_plots[[1]] +
    plot_layout(design = layout,
                heights = c(1, 1.4, 1, 1.4, 1, 1.4, 1, 1.4),
                widths = 2,
                guides = "auto") +
    plot_annotation(title = "124 ms Component",
                    theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

rc11_component_plots <-
  all_watch_topos[[2]][[1]] +
  all_watch_topos[[2]][[2]] +
  all_watch_topos[[2]][[3]] +
  watch_plots[[2]] +
  all_pos_topos[[2]][[1]] +
  all_pos_topos[[2]][[2]] +
  all_pos_topos[[2]][[3]] +
  pos_plots[[2]] +
  all_neg_topos[[2]][[1]] +
  all_neg_topos[[2]][[2]] +
  all_neg_topos[[2]][[3]] +
  neg_plots[[2]] +
  all_reg_topos[[2]][[1]] +
  all_reg_topos[[2]][[2]] +
  reg_plots[[2]] +
  plot_layout(design = layout,
              heights = c(1, 1.4, 1, 1.4, 1, 1.4, 1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "162 ms Component",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

neg_rc12_component_plots <-
  all_watch_topos[[3]][[1]] +
  all_watch_topos[[3]][[2]] +
  all_watch_topos[[3]][[3]] +
  watch_plots[[3]] +
  all_pos_topos[[3]][[1]] +
  all_pos_topos[[3]][[2]] +
  all_pos_topos[[3]][[3]] +
  pos_plots[[3]] +
  all_neg_topos[[3]][[1]] +
  all_neg_topos[[3]][[2]] +
  all_neg_topos[[3]][[3]] +
  neg_plots[[3]] +
  all_reg_topos[[3]][[1]] +
  all_reg_topos[[3]][[2]] +
  reg_plots[[3]] +
  plot_layout(design = layout,
              heights = c(1, 1.4, 1, 1.4, 1, 1.4, 1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "Negative 259 ms Component",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

pos_rc12_component_plots <-
  all_watch_topos[[4]][[1]] +
  all_watch_topos[[4]][[2]] +
  all_watch_topos[[4]][[3]] +
  watch_plots[[4]] +
  all_pos_topos[[4]][[1]] +
  all_pos_topos[[4]][[2]] +
  all_pos_topos[[4]][[3]] +
  pos_plots[[4]] +
  all_neg_topos[[4]][[1]] +
  all_neg_topos[[4]][[2]] +
  all_neg_topos[[4]][[3]] +
  neg_plots[[4]] +
  all_reg_topos[[4]][[1]] +
  all_reg_topos[[4]][[2]] +
  reg_plots[[4]] +
  plot_layout(design = layout,
              heights = c(1, 1.4, 1, 1.4, 1, 1.4, 1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "Positive 259 ms Component",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

rc2_component_plots <-
  all_watch_topos[[5]][[1]] +
  all_watch_topos[[5]][[2]] +
  all_watch_topos[[5]][[3]] +
  watch_plots[[5]] +
  all_pos_topos[[5]][[1]] +
  all_pos_topos[[5]][[2]] +
  all_pos_topos[[5]][[3]] +
  pos_plots[[5]] +
  all_neg_topos[[5]][[1]] +
  all_neg_topos[[5]][[2]] +
  all_neg_topos[[5]][[3]] +
  neg_plots[[5]] +
  all_reg_topos[[5]][[1]] +
  all_reg_topos[[5]][[2]] +
  reg_plots[[5]] +
  plot_layout(design = layout,
              heights = c(1, 1.4, 1, 1.4, 1, 1.4, 1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "381 ms Component",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

rc3_component_plots <-
  all_watch_topos[[6]][[1]] +
  all_watch_topos[[6]][[2]] +
  all_watch_topos[[6]][[3]] +
  watch_plots[[6]] +
  all_pos_topos[[6]][[1]] +
  all_pos_topos[[6]][[2]] +
  all_pos_topos[[6]][[3]] +
  pos_plots[[6]] +
  all_neg_topos[[6]][[1]] +
  all_neg_topos[[6]][[2]] +
  all_neg_topos[[6]][[3]] +
  neg_plots[[6]] +
  all_reg_topos[[6]][[1]] +
  all_reg_topos[[6]][[2]] +
  reg_plots[[6]] +
  plot_layout(design = layout,
              heights = c(1, 1.4, 1, 1.4, 1, 1.4, 1, 1.4),
              widths = 2,
              guides = "auto") +
  plot_annotation(title = "740 ms Component",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))
