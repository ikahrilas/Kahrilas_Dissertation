## Deriving subscales
## Author: Ian J. Kahrirlas
## This script calculates subscales and reliability coeffcients (McDonald's Omega) for analyses. 
## This script should be run after the "Loading and Cleaning data.R" script.

# load packages
library(MBESS) # install.packages("MBESS")
library(tidyverse) # install.packages("tidyverse")

# load dataset 
sad_data <- read_csv("total_data.csv")

# reverse score SBI itmes. Items are a 7-item likert scale, so subtract each score from 8.
sbi_rev_items <- paste("SBI_", seq(2, 24, by = 2), sep = "") #variable containing SBI items to be reversed, which is every other item from 2 to 24
# create names of reversed items to be used in for loop
sbi_rev_names <- paste0(sbi_rev_items, "_r")
# define reverse scoring function that can be used for any measure
reverse <- function(item, subtraction) {
  subtraction - item
} 
# for loop that derives reverse scored SBI items
for (i in seq_along(sbi_rev_items)) {
  sad_data[, sbi_rev_names[i]] <- reverse(sad_data[, sbi_rev_items[i]], 8)
}

# derive SBI subscales of anticipating, savoring the moment, and anticipating by averaging the appropriate items.
#establish variable containing the number of items in each average, which is 8.
sbi_sub_n_items <- 8

# change variable names to lower case so the below code works
names(sad_data) <- tolower(names(sad_data))

# the following code derives the subscales and adds these variables to the per dataset. 
# Subscales are derived by averaging select items.
sad_data <- sad_data %>%
  mutate(anticipating = (sbi_1 + sbi_7 + sbi_13 + sbi_19 + sbi_4_r + sbi_10_r + sbi_16_r + sbi_22_r),
         savoring_moment = (sbi_5 + sbi_11 + sbi_17 + sbi_23 + sbi_2_r + sbi_8_r + sbi_14_r + sbi_20_r),
         reminiscing  = (sbi_3 + sbi_9 + sbi_15 + sbi_21 + sbi_6_r + sbi_12_r + sbi_18_r + sbi_24_r),
         sbi_tot = (sbi_1 + sbi_7 + sbi_13 + sbi_19 + sbi_4_r + sbi_10_r + sbi_16_r + sbi_22_r + sbi_5 + sbi_11 + sbi_17 + sbi_23 + sbi_2_r + sbi_8_r + sbi_14_r + sbi_20_r + sbi_3 + sbi_9 + sbi_15 + sbi_21 + sbi_6_r + sbi_12_r + sbi_18_r + sbi_24_r) 
  )

# derive McDonald's Omega for each subscale
anticipating_items <- c("sbi_1", "sbi_7", "sbi_13", "sbi_19", "sbi_4_r", "sbi_10_r", "sbi_16_r", "sbi_22_r")
moment_items <- c("sbi_5", "sbi_11", "sbi_17", "sbi_23", "sbi_2_r", "sbi_8_r", "sbi_14_r", "sbi_20_r")
reminiscing_items <- c("sbi_3", "sbi_9", "sbi_15", "sbi_21", "sbi_6_r", "sbi_12_r", "sbi_18_r", "sbi_24_r")

items <- list(anticipating = select(sad_data, anticipating_items), 
              savoring_moment = select(sad_data, moment_items), 
              reminiscing = select(sad_data, reminiscing_items)
)

map(items, ~{
  ci.reliability(.x, type = "omega")
})

# derive positive affectivity, negative affectivity, and anxious arousal MASQ subscales, 
# which are calculated by simply summing select items
sad_data <- sad_data %>%
  mutate(masq_pa = masq_2 + masq_4 + masq_5 + masq_7 + masq_11 + masq_14 + masq_19 + masq_23 + masq_26 + masq_28 + masq_32 + masq_34 + masq_36 + masq_37,
         masq_na = masq_9 + masq_13 + masq_17 + masq_21 + masq_29 + masq_30 + masq_35 + masq_38,
         masq_aa = masq_1 + masq_3 + masq_6 + masq_8 + masq_10 + masq_12 + masq_15 + masq_16 + masq_18 + masq_20 + masq_22 + masq_24 + masq_25 + masq_27 + masq_31 + masq_33 + masq_39)

# derive McDonald's Omega for each subscale
masq_pa_items <- c("masq_2", "masq_4", "masq_5", "masq_7", "masq_11", "masq_14", "masq_19", "masq_23", "masq_26", "masq_28", "masq_32", "masq_34", "masq_36", "masq_37")
masq_na_items <- c("masq_9", "masq_13", "masq_17", "masq_21", "masq_29", "masq_30", "masq_35", "masq_38")
masq_aa_items <- c("masq_1", "masq_3", "masq_6", "masq_8", "masq_10", "masq_12", "masq_15", "masq_16", "masq_18", "masq_20", "masq_22", "masq_24", "masq_25", "masq_27", "masq_31", "masq_33", "masq_39")

items_masq <- list(pa = select(sad_data, masq_pa_items),
                   na = select(sad_data, masq_na_items),
                   aa = select(sad_data, masq_aa_items)
)

map(items_masq, ~{
  ci.reliability(.x, type = "omega")
})

# derive the total for the PHQ-9
sad_data <- sad_data %>%
  mutate(phq_total = phq_1 + phq_2 + phq_3 + phq_4 + phq_5 + phq_6 + phq_7 + phq_8 + phq_9)

# McDonald's omega
phq_items <- c("phq_1", "phq_2", "phq_3", "phq_4", "phq_5", "phq_6", "phq_7", "phq_8", "phq_9")
ci.reliability(select(sad_data, phq_items), type = "omega")

# derive total for PSWQ
sad_data <- sad_data %>%
  mutate(pswq_total = pswq_1 + pswq_2 + pswq_3 + pswq_4 + pswq_5 + pswq_6 + pswq_7 + pswq_8
         + pswq_9 + pswq_10 + pswq_11 + pswq_12 + pswq_13 + pswq_14 + pswq_15 + pswq_16)

# McDonalds omega
pswq_items <- paste0("pswq_", 1:16)
ci.reliability(select(sad_data, pswq_items), type = "omega")
 
# save data
write_csv(sad_data, "total_data.csv")