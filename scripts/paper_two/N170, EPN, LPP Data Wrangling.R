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
n170_epn_files <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N_EPN mul/", list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N_EPN mul/"))
evt_files <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N_EPN_LPP evt/", list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N_EPN_LPP evt/"))
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
           n_trials = as.numeric(n_trials),
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
) %>%
  select(-c(E67:X74, X71:X72)) # drop extraneous columns

# this reads all average-referenced files into a list and generates, pid, block, and ms variables
avr <- map_df(n170_epn_files, ~ {
  read_table2(.x, skip = 1) %>%
    mutate(
      pid = as.numeric(str_extract(.x, "[0-9]{7,}")),
      block = rep(block_names, each = nrow(.) / 7),
      ms = rep(seq(from = -200, to = 3000,
                   by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
               times = 7))
}
)
#'
#' Merge split case (206201821b) with rest of the data
#+ merging split case
split_case_mr <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/LPP mul/", str_subset(list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/LPP mul/"), "206201831"))
split_case_avr <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N_EPN mul/", str_subset(list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N_EPN mul/"), "206201831"))

split_case_mr <- bind_rows(read_table2(split_case_mr[1], skip = 1), read_table2(split_case_mr[2], skip = 1)) %>%
  mutate(
    pid = as.numeric(str_extract(split_case_mr[1], "[0-9]{7,}")),
    block = rep(block_names, each = nrow(.) / 7),
    ms = rep(seq(from = -200, to = 3000,
                 by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
             times = 7)) %>%
  select(-c(E67:X74))

mastoid <- bind_rows(mastoid, split_case_mr)

# average reference
split_case_avr <- bind_rows(read_table2(split_case_avr[1], skip = 1), read_table2(split_case_avr[2], skip = 1)) %>%
  mutate(
    pid = as.numeric(str_extract(split_case_avr[1], "[0-9]{7,}")),
    block = rep(block_names, each = nrow(.) / 7),
    ms = rep(seq(from = -200, to = 3000,
                 by = ((3200 + (3200 / (nrow(.)/7))) / (nrow(.)/7))),
             times = 7))

avr <- bind_rows(avr, split_case_avr)

# evt data
split_case_evt <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N_EPN_LPP evt/", str_subset(list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/N_EPN_LPP evt/"), "206201831"))

split_case_evt <- bind_rows(read_table(split_case_evt[1]),
                            read_table(split_case_evt[2])) %>%
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
#'
#' Clean up variable names
#+ cleaning
# clean up electrode names
names(avr) <- gsub("_.*", "", names(avr))
names(mastoid) <- gsub("_.*", "", names(mastoid))
#'
#' merge all data together and write to work space
#+ merge .mul and .evt data frames
erp_mast <- inner_join(mastoid, evt, by = c("pid", "block"))
erp_avr <- inner_join(avr, evt, by = c("pid", "block"))

# fix incorrect pid
erp_mast$pid[erp_mast$pid == 201206832] <- 206201832
erp_avr$pid[erp_avr$pid == 201206832] <- 206201832

# set pid as character class
erp_mast$pid <- as.character(erp_mast$pid)
erp_avr$pid <- as.character(erp_avr$pid)

# write files
write_csv(erp_mast, here("data", "paper_two", "created_data", "erp_mast.csv"))
write_csv(erp_avr, here("data", "paper_two", "created_data", "erp_avr.csv"))
