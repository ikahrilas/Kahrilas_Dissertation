#' ---
#' title: "Wrangle N170/EPN/LPP Data"
#' author: "Ian J. Kahrilas"
#' date: "2020/6/12"
#' output: "html_document"
#' ---
#' Load in packages
#+ set up
library(tidyverse)
library(here)
library(glue)
#'
#' Create a vectors with .mul and .evt names
#+ files names
LPP_files <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/LPP mul/", list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/LPP mul/"))
n170_epn_files <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N170_EPN mul/", list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N170_EPN mul/"))
evt_files <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N170_EPN_LPP evt/", list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N170_EPN_LPP evt/"))
## get rid of bad files
LPP_files <- str_subset(LPP_files, "206201831", negate = TRUE)
LPP_files <- str_subset(LPP_files, "206201827", negate = TRUE)
LPP_files <- str_subset(LPP_files, "20620183_av", negate = TRUE)
n170_epn_files <- str_subset(n170_epn_files, "206201831", negate = TRUE)
n170_epn_files <- str_subset(n170_epn_files, "206201827", negate = TRUE)
n170_epn_files <- str_subset(n170_epn_files, "20620183_av", negate = TRUE)
evt_files <- str_subset(evt_files, "206201827", negate = TRUE)
evt_files <- str_subset(evt_files, "206201831", negate = TRUE)
evt_files <- str_subset(evt_files, "20620183_av", negate = TRUE)
#'
#' Read in .mul and .evt files
#+ read data in
# make list of block names
block_names <- c("Pos_Inc",
                 "Pos_Dec",
                 "Pos_Watch",
                 "Neg_Inc",
                 "Neg_Dec",
                 "Neg_Watch",
                 "Neu_Watch")

# read in .evt files, make columns for block, number of trials, and pid. This cleans names up, too.
evt <- map_dfr(evt_files, ~{
  read_table(.x) %>%
    separate(`Code	TriNo	Comnt`, into = c("block", "n_trials"), sep = " ") %>%
    mutate(block = str_extract(block, block_names),
           n_trials = as.numeric(str_extract(n_trials, "[0-9]{2}")),
           pid = as.numeric(str_extract(.x, "[0-9]+")),
           prop_trials = n_trials / 40)
})

# this reads all mastoid-referenced files into a list and generates pid, block, and ms variables
mastoid <- map_dfr(LPP_files, ~ {
  read_table2(.x, skip = 1) %>%
    mutate(
      pid = as.numeric(str_extract(.x, "[0-9]{7,}")),
      block = rep(block_names, each = nrow(.) / 7),
      ms = rep(seq(from = -200, to = 3000,
                   by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
               times = 7),
    )
}
)

# drop extraneous columns
mastoid <- mastoid %>% select(-c(E67:X74, X71:X72))

# this reads all average-referenced files into a list and generates, pid, block, and ms variables
avg_lp <- map_df(here("data", "paper_two", "data_from_PER_R_project",  "avr_lp", files_avg_lp), ~ {
  read_table2(.x, skip = 1) %>%
    mutate(
      pid = as.numeric(str_extract(.x, "[0-9]{7,}")),
      block = rep(block_names, each = nrow(.) / 7),
      ms = rep(seq(from = -200, to = 3000,
                   by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
               times = 7))
}
)
avg_no_lp <- map_df(here("data", "paper_two", "data_from_PER_R_project", "avr_no_lp", files_avg_no_lp), ~ {
  read_table2(.x, skip = 1) %>%
    mutate(
      pid = as.numeric(str_extract(.x, "[0-9]{7,}")),
      block = rep(block_names, each = nrow(.) / 7),
      ms = rep(seq(from = -200, to = 3000,
                   by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
               times = 7))
}
)

pre_block <- c("Pre_Pos_Inc",
               "Pre_Pos_Dec",
               "Pre_Pos_Watch",
               "Pre_Neu_Watch",
               "Pre_Neg_Inc",
               "Pre_Neg_Dec",
               "Pre_Neg_Watch")

spn <- map_df(here("data", "paper_two", "SPN_mul", files_spn), ~ {
  read_table2(.x, skip = 1) %>%
    mutate(
      pid = as.numeric(str_extract(.x, "[0-9]{7,}")),
      block = rep(pre_block, each = nrow(.) / 7),
      ms = rep(seq(from = -200, to = 2000,
                   by = ((2200 + (2200 / (nrow(.)/7))) / (nrow(.)/7))),
               times = 7))
}
)

spn_evt <- map_df(files_spn_evt, ~{
  read_table(here("data", "paper_two", "SPN_evt", .x)) %>%
    separate(`Code	TriNo	Comnt`, into = c("block", "n_trials"), sep = " ") %>%
    mutate(block = str_extract(block, pre_block),
           n_trials = as.numeric(str_extract(n_trials, "[0-9]{2}")),
           pid = as.numeric(str_extract(.x, "[0-9]+")),
           prop_trials = n_trials / 40)
})
```

Merge split case with rest of the data

```{r take care of split case}
# mastoid reference
split_case_mr_lp <- bind_rows(read_table2(here("data","paper_two", "data_from_PER_R_project", "mast_lp", split_case_mast_lp[1]), skip = 1), read_table2(here("data", "paper_two", "data_from_PER_R_project", "mast_lp", split_case_mast_lp[2]), skip = 1)) %>%
  mutate(
    pid = as.numeric(str_extract(split_case_mast_lp[1], "[0-9]{7,}")),
    block = rep(block_names, each = nrow(.) / 7),
    ms = rep(seq(from = -200, to = 3000,
                 by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
             times = 7)) %>%
  select(-c(E67:X74))

mastoid_lp <- bind_rows(mastoid_lp, split_case_mr_lp)

split_case_mr_no_lp <- bind_rows(read_table2(here("data","paper_two", "data_from_PER_R_project",  "mast_no_lp", split_case_mast_no_lp[1]), skip = 1), read_table2(here("data", "paper_two", "data_from_PER_R_project",  "mast_no_lp", split_case_mast_no_lp[2]), skip = 1)) %>%
  mutate(
    pid = as.numeric(str_extract(split_case_mast_no_lp[1], "[0-9]{7,}")),
    block = rep(block_names, each = nrow(.) / 7),
    ms = rep(seq(from = -200, to = 3000,
                 by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
             times = 7)) %>%
  select(-c(E67:X74))

mastoid_no_lp <- bind_rows(mastoid_no_lp, split_case_mr_no_lp)

# average reference
split_case_ar_lp <- bind_rows(read_table2(here("data", "paper_two", "data_from_PER_R_project",  "avr_lp", split_case_avr_lp[1]), skip = 1), read_table2(here("data", "paper_two", "data_from_PER_R_project",  "avr_lp", split_case_avr_lp[2]), skip = 1)) %>%
  mutate(
    pid = as.numeric(str_extract(split_case_avr_lp[1], "[0-9]{7,}")),
    block = rep(block_names, each = nrow(.) / 7),
    ms = rep(seq(from = -200, to = 3000,
                 by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
             times = 7))

avg_lp <- bind_rows(avg_lp, split_case_ar_lp)

split_case_ar_no_lp <- bind_rows(read_table2(here("data", "paper_two", "data_from_PER_R_project",  "avr_no_lp", split_case_avr_no_lp[1]), skip = 1), read_table2(here("data", "paper_two", "data_from_PER_R_project", "avr_no_lp", split_case_avr_no_lp[2]), skip = 1)) %>%
  mutate(
    pid = as.numeric(str_extract(split_case_avr_no_lp[1], "[0-9]{7,}")),
    block = rep(block_names, each = nrow(.) / 7),
    ms = rep(seq(from = -200, to = 3000,
                 by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
             times = 7))

avg_no_lp <- bind_rows(avg_no_lp, split_case_ar_no_lp)

# evt data
split_case_evt <- bind_rows(read_table(here("data", "paper_two", "data_from_PER_R_project",  "mast_evt", split_case_evt[1])),
                            read_table(here("data", "paper_two", "data_from_PER_R_project",  "mast_evt", split_case_evt[2]))) %>%
  separate(`Code	TriNo	Comnt`, into = c("block", "n_trials"), sep = " ") %>%
  mutate(block = c("Pos_Dec",
                   "Neg_Inc",
                   "Neu_Watch",
                   "Pos_Inc",
                   "Pos_Watch",
                   "Neg_Dec",
                   "Neg_Watch"),
         n_trials = as.numeric(n_trials),
         pid = as.numeric(str_extract(split_case_evt[1], "[0-9]+")),
         prop_trials = n_trials / 40)
evt <- bind_rows(evt, split_case_evt)

# bind_rows(read_table2(here("data", "paper_two", "SPN_mul", split_case_spn[1]), skip = 1), read_table2(here("data", "paper_two", "SPN_mul", split_case_spn[2]), skip = 1)) %>%
#   nrow()
#   mutate(
#     pid = as.numeric(str_extract(split_case_spn[1], "[0-9]{7,}")),
#     block = rep(pre_block, each = nrow(.) / 7),
#     ms = rep(seq(from = -200, to = 2000,
#                    by = ((2200 + (2200 / (nrow(.)/7))) / (nrow(.)/7))),
#                times = 7))
```

```{r clean up variable/condition names}
# clean up electrode names
names(avg_lp) <- gsub("_.*", "", names(avg_lp))
names(avg_no_lp) <- gsub("_.*", "", names(avg_no_lp))
names(mastoid_lp) <- gsub("_.*", "", names(mastoid_lp))
names(mastoid_no_lp) <- gsub("_.*", "", names(mastoid_no_lp))
names(spn) <- gsub("_.*", "", names(spn))
```

```{r merge mul and evt dataframes}
# preallocate space
erp_mast_lp <- as_tibble(matrix(data = NA_real_, nrow = nrow(mastoid_lp), ncol = ncol(mastoid_lp) + ncol(evt)))
erp_mast_no_lp <- as_tibble(matrix(data = NA_real_, nrow = nrow(mastoid_no_lp), ncol = ncol(mastoid_no_lp) + ncol(evt)))
erp_avr_lp <- as_tibble(matrix(data = NA_real_, nrow = nrow(avg_lp), ncol = ncol(avg_lp) + ncol(evt)))
erp_avr_no_lp <- as_tibble(matrix(data = NA_real_, nrow = nrow(avg_no_lp), ncol = ncol(avg_no_lp) + ncol(evt)))
spn_dat <- as_tibble(matrix(data = NA_real_, nrow = nrow(spn), ncol = ncol(spn) + ncol(spn_evt)))
# merge
erp_mast_lp <- inner_join(mastoid_lp, evt, by = c("pid", "block"))
erp_mast_no_lp <- inner_join(mastoid_no_lp, evt, by = c("pid", "block"))
erp_avr_lp <- inner_join(avg_lp, evt, by = c("pid", "block"))
erp_avr_no_lp <- inner_join(avg_no_lp, evt, by = c("pid", "block"))
spn_dat <- inner_join(spn, spn_evt, by = c("pid", "block"))
# fix incorrect pid
erp_mast_lp$pid[erp_mast_lp$pid == 201206832] <- 206201832
erp_mast_no_lp$pid[erp_mast_no_lp$pid == 201206832] <- 206201832
erp_avr_lp$pid[erp_avr_lp$pid == 201206832] <- 206201832
erp_avr_no_lp$pid[erp_avr_no_lp$pid == 201206832] <- 206201832
spn_dat$pid[spn_dat$pi == 201206832] <- 206201832
# set pid as character class
erp_mast_lp$pid <- as.character(erp_mast_lp$pid)
erp_mast_no_lp$pid <- as.character(erp_mast_no_lp$pid)
erp_avr_lp$pid <- as.character(erp_avr_lp$pid)
erp_avr_no_lp$pid <- as.character(erp_avr_no_lp$pid)
spn_dat$pid <- as.character(spn_dat$pid)
# write files
write_csv(erp_mast_lp, here("data", "created_data", "erp_mast_lp.csv"))
write_csv(erp_mast_no_lp, here("data", "created_data", "erp_mast_no_lp.csv"))
write_csv(erp_avr_lp, here("data", "created_data", "erp_avr_lp.csv"))
write_csv(erp_avr_no_lp, here("data", "created_data", "erp_avr_no_lp.csv"))
write_csv(spn_dat, here("data", "paper_two", "created_data", "spn.csv"))
```
