#' ---
#' title: "Wrangle SPN Data"
#' author: "Ian J. Kahrilas"
#' date: "2020/6/4"
#' output: "html_document"
#' ---
# load packages
library(tidyverse)
library(here)
library(glue)

# file names
files <- str_subset(list.files(path = "data/paper_two/SPN mul"), "206201831", negate = TRUE)
files <- glue("data/paper_two/SPN mul/{files}")
# evt file names
files_evt <- str_subset(list.files(path = "data/paper_two/SPN evt"), "206201831", negate = TRUE)
files_evt <- glue("data/paper_two/SPN evt/{files_evt}")

# block names
block_names <- c("Pre_Pos_Inc",
                 "Pre_Pos_Dec",
                 "Pre_Pos_Watch",
                 "Pre_Neu_Watch",
                 "Pre_Neg_Inc",
                 "Pre_Neg_Dec",
                 "Pre_Neg_Watch")

# read in data
spn_dat <- files %>%
  map_dfr(~{
    read_table2(.x, skip = 1) %>%
      mutate(pid = as.numeric(str_extract(.x, "[0-9]{7,}")),
             block = rep(block_names, each = nrow(.) / 7),
             ms = rep(seq(from = -200, to = 2000,
                          by = ((2200 + (2200 / (nrow(.)/7))) / (nrow(.)/7))),
                      times = 7)
             )
  }) %>%
  select(-X71, -X72, -X73, -X74)

# read in evt files
spn_evt <- files_evt %>%
  map_dfr(~{
    read_delim(.x, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      separate(Comnt, into = c("block", "n_trials"), sep = " ") %>%
      mutate(block = str_extract(block, block_names),
             n_trials = as.numeric(n_trials),
             pid = as.numeric(str_extract(.x, "[0-9]+")),
             prop_trials = n_trials / 40)
  })

# handle split file - 206201831
tmp <- bind_rows(read_table2("data/paper_two/SPN mul/206201831_av-export.mul", skip = 1),
          read_table2("data/paper_two/SPN mul/206201831b_av-export.mul", skip = 1)) %>%
  mutate(pid = 206201831,
         block = rep(c("Pos_Dec",
                       "Neg_Inc",
                       "Neu_Watch",
                       "Pos_Inc",
                       "Pos_Watch",
                       "Neg_Dec",
                       "Neg_Watch"), each = nrow(.) / 7),
         ms = rep(seq(from = -200, to = 2000,
                      by = ((2200 + (2200 / (nrow(.)/7))) / (nrow(.)/7))),
                  times = 7)
  ) %>%
  select(-X74)

# split evt files
split_case_evt <- paste0("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/SPN evt/", str_subset(list.files("/Users/ian/tmp/Kahrilas_Dissertation/data/paper_two/SPN evt/"), "206201831"))

split_case_evt <- bind_rows(read_delim(split_case_evt[1], "\t", escape_double = FALSE, trim_ws = TRUE),
                            read_delim(split_case_evt[2], "\t", escape_double = FALSE, trim_ws = TRUE)) %>%
  separate(Comnt, into = c("block", "n_trials"), sep = " ") %>%
  mutate(block = c("Pos_Dec",
                   "Neg_Inc",
                   "Neu_Watch",
                   "Pos_Inc",
                   "Pos_Watch",
                   "Neg_Dec",
                   "Neg_Watch"),
         pid = as.numeric(str_extract(split_case_evt[1], "[0-9]+")),
         n_trials = as.numeric(n_trials),
         prop_trials = n_trials / 40)

evt <- bind_rows(evt, split_case_evt)

# merge split file with rest of spn data
spn_dat <- bind_rows(spn_dat, tmp)
# clean up names
names(spn_dat) <- gsub("_.*", "", names(spn_dat))

# merge mul and evt files
spn <- full_join(spn_dat, spn_evt, by = c("pid", "block")) %>%
  select(-Tmu, -Code, -TriNo)

# write file
write_csv(spn, here("data", "paper_two", "created_data", "spn.csv"))
