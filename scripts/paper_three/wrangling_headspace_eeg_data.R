## reading in headspace t1 EEG data

# load packages
library(tidyverse)
library(here)

# vector of all files in working directory partitioned by avg montage and evt files
files <- list.files(here("data", "paper_three", "headspace_mul_files"))
files_evt <- list.files(here("data", "paper_three", "headspace_evt_files"))

# remove 22585472, 22585472b, 22585452, 22585473, and 22585488 as these files will have to be handled separately
files <- files[!str_detect(files, c("22585472_T1_a_av-export.mul"))]
files <- files[!str_detect(files, c("22585472_T1_av-export.mul"))]
files <- files[!str_detect(files, c("22585452_T1_av-export.mul"))]
files <- files[!str_detect(files, c("22585473_T1_av-export.mul"))]
files <- files[!str_detect(files, c("22585488_T1_av-export.mul"))]

# Read in .mul and .evt files
## preallocate space
avg_dat <- as_tibble(matrix(data = NA_real_, nrow = 580454, ncol = 80))
evt <- tibble(
  Tmu = rep(NA, times = 51*7),
  block = rep(NA, times = 51*7),
  n_trials = rep(NA, times = 51*7),
  pid = rep(NA, times = 51*7),
  prop_trials = rep(NA, times = 51*7)
)

# make list of block names
block_names <- c("Pos_Inc",
                 "Pos_Dec",
                 "Pos_Watch",
                 "Neg_Inc",
                 "Neg_Dec",
                 "Neg_Watch",
                 "Neu_Watch")

# read in evt files, make columns for block, number of trials, and pid. This cleans names up, too.
evt <- map_df(files_evt, ~ {
  read_table(here("data", "paper_three", "headspace_evt_files", .x)) %>%
    select(-Tmu) %>%
    rename("block:n_trials" = `Code\tTriNo\tComnt`) %>%
    mutate(`block:n_trials` = str_remove(`block:n_trials`, "42\t200000\t")) %>%
    separate(`block:n_trials`, into = c("block", "n_trials"), sep = ":") %>%
    mutate(n_trials = str_remove(n_trials, "avs"),
           n_trials = as.numeric(n_trials),
           pid = as.numeric(str_extract(.x, "[0-9]+")))
})

# this reads all average-referenced files into a list and generates, pid, block, and ms variables
eeg_dat <- map_df(here("data", "paper_three", "headspace_mul_files",  files), ~ {
  read_table2(.x, skip = 1) %>%
    mutate(
      pid = as.numeric(str_extract(.x, "[0-9]{7,}")),
      block = rep(block_names, each = nrow(.) / 7),
      ms = rep(seq(from = -200, to = 3000,
                   by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
               times = 7)) %>%
    relocate(pid, block, ms, everything())
}
)
length(unique(eeg_dat$pid))

# Merge lingering cases with the rest of the data
## r take care of split case
split_file <- c("22585472_T1_a_av-export.mul", "22585472_T1_av-export.mul")

# average reference
eeg_dat <-
  bind_rows(
    read_table2(here("data", "paper_three", "headspace_mul_files", split_file[1]), skip = 1),
    read_table2(here("data", "paper_three", "headspace_mul_files", split_file[2]), skip = 1)) %>%
  mutate(
    pid = as.numeric(str_extract(split_file[1], "[0-9]{7,}")),
    block = rep(block_names, each = nrow(.) / 7),
    ms = rep(seq(from = -200, to = 3000,
                 by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
             times = 7)) %>%
  bind_rows(eeg_dat)

# handle cases with only six blocks
files_with_six_blocks <- c("22585452_T1_av-export.mul",
                           "22585473_T1_av-export.mul",
                           "22585488_T1_av-export.mul")

eeg_dat <-
  map_df(here("data", "paper_three", "headspace_mul_files",  files_with_six_blocks), ~ {
  read_table2(.x, skip = 1) %>%
  mutate(
    pid = as.numeric(str_extract(.x, "[0-9]{7,}")),
    block = rep(evt %>%
                  filter(pid == as.numeric(str_extract(.x, "[0-9]{7,}"))) %>%
                  select(block) %>%
                  pull(),
                each = nrow(.) / 6),
    ms = rep(seq(from = -200, to = 3000,
                 by = ((3200 + (3200 / (nrow(.)/6))) / (nrow(.)/6))),
             times = 6)) %>%
    relocate(pid, block, ms, everything())
}) %>%
  bind_rows(eeg_dat)

# clean up electrode names
names(eeg_dat) <- gsub("_.*", "", names(eeg_dat))

# merge eeg and evt data
eeg_evt_dat <- full_join(eeg_dat, evt, by = c("pid", "block")) %>%
  relocate(pid, block, n_trials)

write_csv(eeg_evt_dat, file = here("data", "paper_three", "erp_dat.csv"))
