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
files <- str_subset(list.files(path = "data/paper_two/SPN_mul"), "206201831", negate = TRUE)
files <- glue("data/paper_two/SPN_mul/{files}")

dat <- files %>%
  map_dfr(~ {

  })
