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

# RC2 as outcome
## relevel block variable
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
## center positive affectivity
dat$pos_affectivity <- scale(dat$pos_affectivity, center = TRUE, scale = FALSE)
rc2_pa_mod <- lmer(RC2 ~ block*pos_affectivity + (1|pid), data = dat,)
## check assumptions
check_model(rc2_pa_mod)
## results
summary(rc2_pa_mod)
### There is no block by PA moderation effect in predicting RC2
## pull stats for paper
rc2_pa_dat <- data.frame(coef(summary(rc2_pa_mod)))
rc2_pa_beta <- rc2_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
rc2_pa_std_beta <- standardize_parameters(rc2_pa_mod) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
rc2_pa_se <- rc2_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
rc2_pa_df <- rc2_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
rc2_pa_t <- rc2_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
rc2_pa_p <- rc2_pa_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# RC3 as outcome
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
rc3_pa_mod <- lmer(RC3 ~ block*pos_affectivity + (1|pid), data = dat)
## check assumptions
check_model(rc3_pa_mod)
## model results
summary(rc3_pa_mod)
### no significant moderation effect, p = .09
## extract model statistcs
rc3_pa_mod_dat <- data.frame(coef(summary(rc3_pa_mod)))
rc3_pa_beta <- rc3_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
rc3_pa_std_beta <- standardize_parameters(rc3_pa_mod) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
rc3_pa_se <- rc3_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
rc3_pa_df <- rc3_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
rc3_pa_t <- rc3_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
rc3_pa_p <- rc3_pa_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# RC5 as outcome
rc5_pa_mod <- lmer(RC5 ~ block*pos_affectivity + (1|pid), data = dat)
## check assumptions
check_model(rc5_pa_mod)
## model results
summary(rc5_pa_mod)
### nothing
## extract model stats
rc5_pa_mod_dat <- data.frame(coef(summary(rc5_pa_mod)))
rc5_pa_beta <- rc5_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
rc5_pa_std_beta <- standardize_parameters(rc5_pa_mod) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
rc5_pa_se <- rc5_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
rc5_pa_df <- rc5_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
rc5_pa_t <- rc5_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
rc5_pa_p <- rc5_pa_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# RC11 as outcome
rc11_pa_mod <- lmer(RC11 ~ block*pos_affectivity + (1|pid), data = dat)
## check model
check_model(rc11_pa_mod)
## model results
summary(rc11_pa_mod)
### nothing
## extract info
rc11_pa_mod_dat <- data.frame(coef(summary(rc11_pa_mod)))
rc11_pa_beta <- rc11_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
rc11_pa_std_beta <- standardize_parameters(rc11_pa_mod) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
rc11_pa_se <- rc11_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
rc11_pa_df <- rc11_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
rc11_pa_t <- rc11_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
rc11_pa_p <- rc11_pa_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# RC12 as outcome
rc12_pa_mod <- lmer(RC12 ~ block*pos_affectivity + (1|pid), data = dat)
## check model
check_model(rc12_pa_mod)
## model results
summary(rc12_pa_mod)
### nothing
## extract info
rc12_pa_mod_dat <- data.frame(coef(summary(rc12_pa_mod)))
rc12_pa_beta <- rc12_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
rc12_pa_std_beta <- standardize_parameters(rc12_pa_mod) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
rc12_pa_se <- rc12_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
rc12_pa_df <- rc12_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
rc12_pa_t <- rc12_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
rc12_pa_p <- rc12_pa_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# pos_RC12 as outcome
pos_RC12_pa_mod <- lmer(pos_RC12 ~ block*pos_affectivity + (1|pid), data = dat)
## check model
check_model(pos_RC12_pa_mod)
## model results
summary(pos_RC12_pa_mod)
### nothing
## extract info
pos_RC12_pa_mod_dat <- data.frame(coef(summary(pos_RC12_pa_mod)))
pos_RC12_pa_beta <- pos_RC12_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
pos_RC12_pa_std_beta <- standardize_parameters(pos_RC12_pa_mod) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
pos_RC12_pa_se <- pos_RC12_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
pos_RC12_pa_df <- pos_RC12_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
pos_RC12_pa_t <- pos_RC12_pa_mod_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
pos_RC12_pa_p <- pos_RC12_pa_mod_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()


# arousal as outcome
ar_pa <- lmer(arousal ~ block*pos_affectivity + (1|pid), data = dat)
summary(ar_pa)
ar_pa_dat <- data.frame(coef(summary(ar_pa)))
ar_pa_beta <- ar_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
ar_pa_std_beta <- standardize_parameters(ar_pa) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
ar_pa_se <- ar_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
ar_pa_df <- ar_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
ar_pa_t <- ar_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
ar_pa_p <- ar_pa_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# valence as outcome
val_pa <- lmer(valence ~ block*pos_affectivity + (1|pid), data = dat)
summary(val_pa)
val_pa_dat <- data.frame(coef(summary(val_pa)))
val_pa_beta <- val_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
val_pa_std_beta <- standardize_parameters(val_pa) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
val_pa_se <- val_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
val_pa_df <- val_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
val_pa_t <- val_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
val_pa_p <- val_pa_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# difficulty as outcome
diff_pa <- lmer(difficulty ~ block*pos_affectivity + (1|pid), data = dat)
summary(diff_pa)
diff_pa_dat <- data.frame(coef(summary(diff_pa)))
diff_pa_beta <- diff_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
diff_pa_std_beta <- standardize_parameters(diff_pa) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
diff_pa_se <- diff_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
diff_pa_df <- diff_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
diff_pa_t <- diff_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
diff_pa_p <- diff_pa_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_pa_factor_score_moderation-data", ".RData"))
