# Repeated measures ANOVA

## load packages
library(tidyverse)
library(lmerTest)
library(here)
library(afex)
library(papaja)

## read in data
dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-18.csv")

########## ANALYSES ##########

dv <- c("RC5", "RC11", "RC12", "pos_RC12", "RC2", "RC3", "valence", "arousal", "difficulty")

aov_lst <- map(dv, ~ {
afex::aov_ez(
  data = dat
  , dv = .x
  , id = "pid"
  , within = "block"
)
})
## EVERYTHING'S SIGNIFICANT
# 124 ms component
rc5_results <- apa_print(aov_lst[[1]])$full_result$block

# 162 ms component
rc11_results <- apa_print(aov_lst[[2]])$full_result$block

# negative 259 component
neg_rc12_results <- apa_print(aov_lst[[3]])$full_result$block

# positive 259 component
pos_rc12_results <- apa_print(aov_lst[[4]])$full_result$block

# 381 component
pos_rc2_results <- apa_print(aov_lst[[5]])$full_result$block

# 740 component
pos_rc3_results <- apa_print(aov_lst[[6]])$full_result$block

# valence
val_results <- apa_print(aov_lst[[7]])$full_result$block

# arousal
ar_results <- apa_print(aov_lst[[8]])$full_result$block

# difficulty
diff_results <- apa_print(aov_lst[[9]])$full_result$block

# save the results to be accessed in the paper
save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_repeated_measures_ANOVA_analysis-data", ".RData"))
