# positive affectivity moderation analyses with factor scores

# load pacakges
library(tidyverse)
library(lmerTest)
library(emmeans)
library(here)
library(performance)
library(effectsize)

# read in data
dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-11.csv")

# relevel block variable
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")

# center savoring variable
dat$savoring_moment <- scale(dat$savoring_moment, center = TRUE, scale = FALSE)

# RC2 as outcome
RC2_stm_mod <- lmer(RC2 ~ block*savoring_moment + (1|pid), data = dat)
## check assumptions
check_model(RC2_stm_mod)
## model results
summary(RC2_stm_mod)
### nothing
## extract model stats
RC2_stm_mod_dat <- data.frame(coef(summary(RC2_stm_mod)))
RC2_stm_mod_beta <- RC2_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
RC2_stm_mod_std_beta <- standardize_parameters(RC2_stm_mod) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
RC2_stm_mod_se <- RC2_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
RC2_stm_mod_df <- RC2_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
RC2_stm_mod_t <- RC2_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
RC2_stm_mod_p <- RC2_stm_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# frontal LPP as outcome
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
RC3_stm_mod <- lmer(RC3 ~ block*savoring_moment + (1|pid), data = dat)
## check assumptions
check_model(RC3_stm_mod)
## model results
summary(RC3_stm_mod)
### nothing
## extract info
RC3_stm_mod_dat <- data.frame(coef(summary(RC3_stm_mod)))
RC3_stm_mod_beta <- RC3_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
RC3_stm_mod_std_beta <- standardize_parameters(RC3_stm_mod) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
RC3_stm_mod_se <- RC3_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
RC3_stm_mod_df <- RC3_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
RC3_stm_mod_t <- RC3_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
RC3_stm_mod_p <- RC3_stm_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# RC5 as outcome
RC5_stm_mod <- lmer(RC5 ~ block*savoring_moment + (1|pid), data = dat)
## check assumptions
check_model(RC5_stm_mod)
## model results
summary(RC5_stm_mod)
### nothing
## extract info
RC5_stm_mod_dat <- data.frame(coef(summary(RC5_stm_mod)))
RC5_stm_mod_beta <- RC5_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
RC5_stm_mod_std_beta <- standardize_parameters(RC5_stm_mod) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
RC5_stm_mod_se <- RC5_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
RC5_stm_mod_df <- RC5_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
RC5_stm_mod_t <- RC5_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
RC5_stm_mod_p <- RC5_stm_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# RC11 as outcome
RC11_stm_mod <- lmer(RC11 ~ block*savoring_moment + (1|pid), data = dat)
## check assumptions
check_model(RC11_stm_mod)
## model results
summary(RC11_stm_mod)
### nothing
## extract information
RC11_stm_mod_dat <- data.frame(coef(summary(RC11_stm_mod)))
RC11_stm_mod_beta <- RC11_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
RC11_stm_mod_std_beta <- standardize_parameters(RC11_stm_mod) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
RC11_stm_mod_se <- RC11_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
RC11_stm_mod_df <- RC11_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
RC11_stm_mod_t <- RC11_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
RC11_stm_mod_p <- RC11_stm_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# RC12 as outcome
RC12_stm_mod <- lmer(RC12 ~ block*savoring_moment + (1|pid), data = dat)
## check assumptions
check_model(RC12_stm_mod)
## model results
summary(RC12_stm_mod)
### nothing
## extract information
RC12_stm_mod_dat <- data.frame(coef(summary(RC12_stm_mod)))
RC12_stm_mod_beta <- RC12_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
RC12_stm_mod_std_beta <- standardize_parameters(RC12_stm_mod) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
RC12_stm_mod_se <- RC12_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
RC12_stm_mod_df <- RC12_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
RC12_stm_mod_t <- RC12_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
RC12_stm_mod_p <- RC12_stm_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# positive RC12 as outcome
pos_RC12_stm_mod <- lmer(pos_RC12 ~ block*savoring_moment + (1|pid), data = dat)
## check assumptions
check_model(pos_RC12_stm_mod)
## model results
summary(pos_RC12_stm_mod)
### nothing
## extract information
pos_RC12_stm_mod_dat <- data.frame(coef(summary(pos_RC12_stm_mod)))
pos_RC12_stm_mod_beta <- pos_RC12_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
pos_RC12_stm_mod_std_beta <- standardize_parameters(pos_RC12_stm_mod) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
pos_RC12_stm_mod_se <- pos_RC12_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
pos_RC12_stm_mod_df <- pos_RC12_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
pos_RC12_stm_mod_t <- pos_RC12_stm_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
pos_RC12_stm_mod_p <- pos_RC12_stm_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# arousal as outcome
ar_stm <- lmer(arousal ~ block*savoring_moment + (1|pid), data = dat)
summary(ar_stm)
ar_stm_dat <- data.frame(coef(summary(ar_stm)))
ar_stm_beta <- ar_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
ar_stm_std_beta <- standardize_parameters(ar_stm) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
ar_stm_se <- ar_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
ar_stm_df <- ar_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
ar_stm_t <- ar_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
ar_stm_p <- ar_stm_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# valence as outcome
val_stm <- lmer(valence ~ block*savoring_moment + (1|pid), data = dat)
summary(val_stm)
val_stm_dat <- data.frame(coef(summary(val_stm)))
val_stm_beta <- val_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
val_stm_std_beta <- standardize_parameters(val_stm) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
val_stm_se <- val_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
val_stm_df <- val_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
val_stm_t <- val_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
val_stm_p <- val_stm_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# save image
save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_stm_factor_score_moderation-data", ".RData"))
