# load packages

library(tidyverse)
library(here)

# read in factor score data
temp_score_dat <- read_csv(here("data", "paper_two", "temp_fac_score_dat.csv"))
temp_spat_score_dat <- read_csv(here("data", "paper_two", "temp_spat_fac_score_dat.csv"))

# read in quesionnaire data
per_sr_dat <- read_csv(here("data", "paper_two", "archive", "data_from_PER_R_project", "created_data", "per_measures.csv"))
