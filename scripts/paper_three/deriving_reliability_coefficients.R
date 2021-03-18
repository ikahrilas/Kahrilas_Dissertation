# deriving reliability coefficients for variables
library(here)
library(MBESS)
library(tidyverse)

# read in data
dat <- read_csv(here("data", "paper_three", "total_questionnaire_data.csv"))

## -- Code for deriving omega reliability coefficients
# make list of all items for sub scales
items = list(masq_pa = c("masq2", "masq4", "masq5", "masq7", "masq11", "masq14",
                         "masq19", "masq23", "masq26", "masq28", "masq32", "masq34",
                         "masq36", "masq37"),
             masq_na = c("masq9", "masq13", "masq17", "masq21", "masq29",
                         "masq30", "masq35", "masq38"),
             masq_aa = c("masq1", "masq3", "masq6", "masq8", "masq10",
                         "masq12", "masq15", "masq16", "masq18", "masq20",
                         "masq22", "masq24", "masq25", "masq27", "masq31", "masq33", "masq39"),
             masq_ad = c("masq2r", "masq4r", "masq5r", "masq7r", "masq11r", "masq14r",
                         "masq19r", "masq23r", "masq26r", "masq28r", "masq32r", "masq34r",
                         "masq36r", "masq37r", "masq9", "masq13", "masq17", "masq21", "masq29",
                         "masq30", "masq35", "masq38"),
             phq_9 = c(paste0("phq", 1:9)),
             pswq = c("pswq1r", "pswq2", "pswq3r", "pswq4", "pswq5", "pswq6", "pswq7",
                      "pswq8r", "pswq9", "pswq10r", "pswq11r", "pswq12", "pswq13", "pswq14",
                      "pswq15", "pswq16"),
             pss = c("pss1", "pss2r", "pss3r", "pss4"),
             panas_pos = c("panas1", "panas3", "panas5", "panas9", "panas10",
                           "panas12", "panas14", "panas16", "panas17", "panas19"),
             panas_neg = c("panas2", "panas4", "panas6", "panas7", "panas8", "panas11",
                           "panas13", "panas15", "panas18", "panas20")
             )

# define function that finds omega for each set of items and prints the coefficient
reliability_fun <- function(items) {
  tmp <- ci.reliability(dat %>% select(all_of(items)), type = "omega")
  round(tmp$est, digits = 2)
}

# loop over all items to derive the coefficients
map_dbl(items, ~ reliability_fun(.x))

