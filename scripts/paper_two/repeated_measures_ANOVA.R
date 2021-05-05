# Repeated measures ANOVA

## load packages
library(tidyverse)
library(lmerTest)
library(here)
library(afex)
library(papaja)
library(kableExtra)

## read in data
dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-18.csv")

########## ANALYSES ##########

dv <- c("RC5", "RC11", "RC12", "RC2", "RC3", "valence", "arousal", "difficulty")

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

apa_print(aov_lst[[1]])$table %>% tibble()

# 162 ms component
rc11_results <- apa_print(aov_lst[[2]])$full_result$block

# negative 259 component
neg_rc12_results <- apa_print(aov_lst[[3]])$full_result$block

# 381 component
pos_rc2_results <- apa_print(aov_lst[[4]])$full_result$block

# 740 component
pos_rc3_results <- apa_print(aov_lst[[5]])$full_result$block

# valence
val_results <- apa_print(aov_lst[[6]])$full_result$block

# arousal
ar_results <- apa_print(aov_lst[[7]])$full_result$block

# difficulty
diff_results <- apa_print(aov_lst[[8]])$full_result$block

# save the results to be accessed in the paper
save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_repeated_measures_ANOVA_analysis-data", ".RData"))

# create table for paper
aov_tab <-
  map_df(aov_lst, ~ {
    apa_print(.x)$table
  }) %>%
    mutate(DV = c("124 ms Component",
                  "162 ms Component",
                  "259 ms Component",
                  "381 ms Component",
                  "740 ms Component",
                  "Valence Ratings",
                  "Arousal Ratings",
                  "Difficulty Ratings")) %>%
    select(-Effect) %>%
    relocate(DV)

names(aov_tab) <- c("DV",
                    "$F$",
                    "$df_1^{GG}$",
                    "$df_2^{GG}$",
                    "$MSE$",
                    "$p$",
                    "$\\hat{\\eta}^2_G$")

aov_tab %>%
  kable("latex", escape = FALSE, booktabs = TRUE, align = c("r", "r", "c", "c", "c", "r", "c"), linesep = "", caption = "(ref:anova-table)") %>%
  row_spec(0, align = "c") %>%
  footnote(escape = FALSE,
           footnote_as_chunk = TRUE,
           general_title = "Note.",
           general = "DV = dependent variable, $F$ = F statistic with Greenhouse-Geisser corrected degrees of freedom, $df_1^{GG}$ = Greenhouse-Geisser corrected numerator degrees of freedom,
           $df_2^{GG}$ = Greenhouse-Geisser corrected denominatr degrees of freedom,
           $MSE$ = mean square error,
           $p$ = p value,
           $\\\\hat{\\\\eta}^2_G$ = partial eta squared as measure of effect size.
           The above table summarized findings from repeated measures ANOVA with block as the independent variable and each of the ERP components and behavioral ratings as
           dependent variables (as indiciated in the DV column).",
           threeparttable = TRUE) %>%
  save_kable(file = "anova_table.pdf")


