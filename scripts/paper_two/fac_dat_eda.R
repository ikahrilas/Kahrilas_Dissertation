# load packages
library(tidyverse)
library(here)
library(lmerTest)

# read in factor score data
temp_score_dat <- read_csv(here("data", "paper_two", "temp_fac_score_dat.csv"))
temp_spat_score_dat <- read_csv(here("data", "paper_two", "temp_spat_fac_score_dat.csv"))

# read in quesionnaire data
per_sr_dat <- read_csv(here("data", "paper_two", "archive", "data_from_PER_R_project", "created_data", "per_measures.csv"))

# merge the data sets and omit variables
dat <- full_join(temp_score_dat, per_sr_dat, by = c("pid", "block")) %>%
  select(pid, age, sex, Race, block:valence, anticipating:masq_aa, tmms_repair:depression)

# define block as a factor variable
dat$block <- factor(dat$block, levels = c("Neg_Inc", "Neg_Dec", "Neg_Watch",
                                          "Neu_Watch",
                                          "Pos_Watch", "Pos_Dec", "Pos_Inc"))

# observe histograms for each of the factors
comps <- c("RC2", "RC3", "RC5", "RC11", "RC12", "RC12")

# define electrodes
rc2_elec <- c("A29", "B26")
rc3_elec <- c("A29", "B26", "A26", "B23",
              "B28", "A30", "B27", "A25", "B22")
rc5_elec <- c("A29", "B26")
rc11_elec <- c("A29", "B26")
rc12_elec <- c("B21", "B28")
pos_rc12_elec <- c("A29", "B26")

## list for function
elec_list <- list(rc2_elec,
                  rc3_elec,
                  rc5_elec,
                  rc11_elec,
                  rc12_elec,
                  pos_rc12_elec)

map2(comps, elec_list, ~ {
  dat %>%
    filter(elec %in% .y) %>%
  ggplot(., aes(.data[[.x]])) +
    geom_histogram(color = "black", fill = "white") +
    facet_wrap(~ block) +
    theme_classic()
})

# observe means and standard deviations
map2(comps, elec_list, ~ {
dat %>%
  filter(elec %in% .y) %>%
  group_by(block) %>%
  summarise(mean(.data[[.x]], na.rm = TRUE),
            sd(.data[[.x]], na.rm = TRUE))
})

# boxplots
map2(comps, elec_list, ~ {
dat %>%
  filter(elec %in% .y) %>%
ggplot(aes(block, .data[[.x]])) +
  geom_boxplot() +
  theme_classic()
})

# all looks fine, so create data set that has factor scores calculated based on appropriate electrodes
# and save the file

dat %>%
  filter(elec %in% rc2_elec) %>%
  group_by(pid, block) %>%
  summarise(RC2 = mean(RC2, na.rm = TRUE)) %>%
  bind_cols(
    dat %>%
      filter(elec %in% rc3_elec) %>%
      group_by(pid, block) %>%
      summarise(RC3 = mean(RC3, na.rm = TRUE))
            ) %>%
  bind_cols(
    dat %>%
      filter(elec %in% rc5_elec) %>%
      group_by(pid, block) %>%
      summarise(RC5 = mean(RC5, na.rm = TRUE))
            ) %>%
  bind_cols(
    dat %>%
      filter(elec %in% rc11_elec) %>%
      group_by(pid, block) %>%
      summarise(RC11 = mean(RC11, na.rm = TRUE))
      ) %>%
  bind_cols(
    dat %>%
      filter(elec %in% rc12_elec) %>%
      group_by(pid, block) %>%
      summarise(RC12 = mean(RC12, na.rm = TRUE))
  ) %>%
  bind_cols(
    dat %>%
      filter(elec %in% pos_rc12_elec) %>%
      group_by(pid, block) %>%
      summarise(pos_RC12 = mean(RC12, na.rm = TRUE))
  ) %>%
  select(pid...1, block...2, RC2, RC3, RC5, RC11, RC12, pos_RC12) %>%
  rename("pid" = pid...1,
         "block" = block...2) %>%
  full_join(per_sr_dat, by = c("pid", "block")) %>%
  select(pid:valence, ethnicity, Race, anticipating:masq_aa, tmms_repair:depression) %>%
  write_csv(file = paste0("data/paper_two/created_data/temp_fac_score_dat_analyses", Sys.Date(), ".csv"))
