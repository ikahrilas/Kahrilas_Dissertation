# load packages
library(tidyverse)
library(here)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(gridExtra)
library(patchwork)
library(ggpubr)

# read in factor score data
temp_score_dat <- read_csv(here("data", "paper_three", "hs_per_fac_score_dat.csv")) %>%
  arrange(pid)
temp_score_dat$pid[temp_score_dat$pid == 22585479] <- 22586479

# derive factor scores using electrode sites of interest
## electrodes for each component
rc2_elec <- c("A29", "B26")
rc3_elec <- c("A29", "B26", "A26", "B23",
              "B28",
              "A30", "B27", "A25", "B22")
rc5_elec <- c("A29", "B26")
rc7_elec <- c("A29", "B26", "A31", "B30")
neg_rc8_elec <- c("B21", "B28")
pos_rc8_elec <- c("A29", "B26", "A31", "B30")
rc17_elec <- c("A29", "B26")

## code for creating component variables
temp_score_dat <-
  temp_score_dat %>%
  filter(elec %in% rc2_elec) %>%
  group_by(pid, block, group) %>%
  summarise(RC2 = mean(RC2, na.rm = TRUE)) %>%
  full_join(
    temp_score_dat %>%
      filter(elec %in% rc3_elec) %>%
      group_by(pid, block, group) %>%
      summarise(RC3 = mean(RC3, na.rm = TRUE)),
    by = c("pid", "block", "group")) %>%
  full_join(
    temp_score_dat %>%
      filter(elec %in% rc5_elec) %>%
      group_by(pid, block, group) %>%
      summarise(RC5 = mean(RC5, na.rm = TRUE)),
    by = c("pid", "block", "group")) %>%
  full_join(
    temp_score_dat %>%
      filter(elec %in% rc7_elec) %>%
      group_by(pid, block, group) %>%
      summarise(RC7 = mean(RC7, na.rm = TRUE)),
    by = c("pid", "block", "group")) %>%
  full_join(
    temp_score_dat %>%
      filter(elec %in% neg_rc8_elec) %>%
      group_by(pid, block, group) %>%
      summarise(nRC8 = mean(RC8, na.rm = TRUE)),
    by = c("pid", "block", "group")) %>%
  full_join(
    temp_score_dat %>%
      filter(elec %in% pos_rc8_elec) %>%
      group_by(pid, block, group) %>%
      summarise(pRC8 = mean(RC8, na.rm = TRUE)),
    by = c("pid", "block", "group")) %>%
  full_join(
    temp_score_dat %>%
      filter(elec %in% rc17_elec) %>%
      group_by(pid, block, group) %>%
      summarise(RC17 = mean(RC17, na.rm = TRUE)),
    by = c("pid", "block", "group"))

# read in questionnaire data
per_hs_ques_dat <- read_csv(here("data", "paper_three", "total_questionnaire_data.csv")) %>%
  arrange(pid)
summary(per_hs_ques_dat)
unique(temp_score_dat$pid[!(temp_score_dat$pid %in% per_hs_ques_dat$pid)])

# merge eeg and questionnaire data
dat <- full_join(temp_score_dat, per_hs_ques_dat, by = c("pid", "group"))

# define block as a factor variable
dat$block <- factor(dat$block, levels = c("Neg_Watch",
                                          "Neu_Watch",
                                          "Pos_Watch"))
dat$block <- droplevels(dat$block)
levels(dat$block)
# observe histograms for each of the factors
comps <- c("RC2", "RC3", "RC5", "RC7", "nRC8", "pRC8", "RC17")

# observe histograms of component data
map(comps, ~ {
dat %>%
  filter(!is.na(RC2)) %>%
  ggplot(., aes(.data[[.x]])) +
  geom_density(color = "black", aes(fill = group), alpha = 0.6) +
  facet_wrap(~ block, nrow = 3) +
  theme_classic() +
  ggtitle(.x)
})

# observe means and standard deviations
map(comps, ~ {
  dat %>%
    group_by(block) %>%
    summarise(mean(.data[[.x]], na.rm = TRUE),
              sd(.data[[.x]], na.rm = TRUE))
})

# boxplots
map(comps, ~ {
  dat %>%
    ggplot(aes(block, .data[[.x]], fill = group)) +
    geom_boxplot() +
    theme_classic()
})


# now inspect univariate distributions of questionnaire data
vars <- c("phq_total",
          "masq_pa",
          "masq_na",
          "masq_aa",
          "masq_ad",
          "pswq_total",
          "panas_pos_total",
          "panas_neg_total",
          "pss_total")

map(vars, ~ {
ggplot(dat, aes(.data[[.x]])) +
  geom_density(color = "black", fill = "red", alpha = 0.6) +
    ggtitle(.x) +
  theme_classic()
})

# bivariate relations
ggscatterhist(dat, x = "masq_pa", y = "masq_na",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_pa", y = "masq_aa",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_na", y = "masq_aa",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_na", y = "phq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_pa", y = "phq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_aa", y = "phq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_ad", y = "phq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_ad", y = "pswq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_na", y = "pswq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_pa", y = "pswq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_ad", y = "pswq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "masq_aa", y = "pswq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "panas_pos_total", y = "panas_neg_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "panas_pos_total", y = "phq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "panas_neg_total", y = "phq_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "panas_pos_total", y = "masq_aa",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

## consider polynomial term if including panas pos and masq_aa in cfa
ggplot(dat, aes(panas_pos_total, masq_aa)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, fill=NA,
              formula=y ~ poly(x, 2), colour="red")

ggscatterhist(dat, x = "panas_pos_total", y = "pss_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat, x = "panas_neg_total", y = "pss_total",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

## -- bivariate relations with EEG components
dat_wide <-
  dat %>%
  pivot_wider(names_from = block,
              values_from = RC2:RC17) %>%
  relocate(pid:race, RC2_Neg_Watch:RC17_NA) %>%
  select(-c(RC2_NA, RC3_NA, RC5_NA, RC7_NA, nRC8_NA, pRC8_NA, RC17_NA))

## RC2
ggscatterhist(dat_wide, x = "RC2_Pos_Watch", y = "RC2_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC2_Pos_Watch", y = "RC2_Neu_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC2_Neu_Watch", y = "RC2_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

## RC3
ggscatterhist(dat_wide, x = "RC3_Pos_Watch", y = "RC3_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC3_Pos_Watch", y = "RC3_Neu_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC3_Neu_Watch", y = "RC3_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

## RC5
ggscatterhist(dat_wide, x = "RC5_Pos_Watch", y = "RC5_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC5_Pos_Watch", y = "RC5_Neu_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC5_Neu_Watch", y = "RC5_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

## RC7
ggscatterhist(dat_wide, x = "RC7_Pos_Watch", y = "RC7_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC7_Pos_Watch", y = "RC7_Neu_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC7_Neu_Watch", y = "RC7_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

## RC8
ggscatterhist(dat_wide, x = "RC8_Pos_Watch", y = "RC8_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC8_Pos_Watch", y = "RC8_Neu_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC8_Neu_Watch", y = "RC8_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

## RC17
ggscatterhist(dat_wide, x = "RC17_Pos_Watch", y = "RC17_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC17_Pos_Watch", y = "RC17_Neu_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

ggscatterhist(dat_wide, x = "RC17_Neu_Watch", y = "RC17_Neg_Watch",
              color = "#00AFBB",
              margin.params = list(fill = "red"))

# save the file
write_csv(x = dat,
          file = paste0("data/paper_three/dat_for_analyses_", Sys.Date(), ".csv"))

