# stm mlm analyses

library(tidyverse)
library(lmerTest)
library(effectsize)

dat <- read_csv("data/paper_two/created_data/per_data_analyses_2020_5_19.csv") %>%
  mutate(savoring_moment = scale(savoring_moment, scale = FALSE, center = TRUE)) %>%
  select(pid:LPP, N170, EPN, everything())

# # SPN as outcome
# dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
# spn_stm <- lmer(SPN ~ block*anticipating + (1|pid), data = dat)
# summary(spn_stm)
# spn_stm_dat <- data.frame(coef(summary(spn_stm)))
# spn_stm_beta <- spn_stm_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(Estimate = sprintf("%.2f", Estimate)) %>%
#   filter(var == "blockPos_Inc:ancitipating") %>%
#   select(Estimate) %>%
#   pull()
# spn_stm_std_beta <- standardize_parameters(spn_stm) %>%
#   filter(Parameter == "blockPos_Inc:ancitipating") %>%
#   mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
#   select(`Std_Coefficient`) %>%
#   pull()
# spn_stm_se <- spn_stm_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
#   filter(var == "blockPos_Inc:ancitipating") %>%
#   select(Std..Error) %>%
#   pull()
# spn_stm_df <- spn_stm_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(df = sprintf("%.2f", df)) %>%
#   filter(var == "blockPos_Inc:ancitipating") %>%
#   select(df) %>%
#   pull()
# spn_stm_t <- spn_stm_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(t.value = sprintf("%.2f", t.value)) %>%
#   filter(var == "blockPos_Inc:ancitipating") %>%
#   select(t.value) %>%
#   pull()
# spn_stm_p <- spn_stm_dat %>%
#   mutate(var = row.names(.),
#          Pr...t.. = sprintf("%.3f", Pr...t..),
#          Pr...t.. = str_remove(Pr...t.., "^0+"),
#          Pr...t.. = paste("=", Pr...t..),
#          Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
#   filter(var == "blockPos_Inc:ancitipating") %>%
#   select(Pr...t..) %>%
#   pull()

# # SPN as outcome reminiscing moderation
# dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
# spn_rem <- lmer(SPN ~ block*reminiscing + (1|pid), data = dat)
# summary(spn_rem)
# spn_rem_dat <- data.frame(coef(summary(spn_rem)))
# spn_rem_beta <- spn_rem_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(Estimate = sprintf("%.2f", Estimate)) %>%
#   filter(var == "blockPos_Inc:reminiscing") %>%
#   select(Estimate) %>%
#   pull()
# spn_rem_std_beta <- standardize_parameters(spn_rem) %>%
#   filter(Parameter == "blockPos_Inc:reminiscing") %>%
#   mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
#   select(`Std_Coefficient`) %>%
#   pull()
# spn_rem_se <- spn_rem_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
#   filter(var == "blockPos_Inc:reminiscing") %>%
#   select(Std..Error) %>%
#   pull()
# spn_rem_df <- spn_rem_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(df = sprintf("%.2f", df)) %>%
#   filter(var == "blockPos_Inc:reminiscing") %>%
#   select(df) %>%
#   pull()
# spn_rem_t <- spn_rem_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(t.value = sprintf("%.2f", t.value)) %>%
#   filter(var == "blockPos_Inc:reminiscing") %>%
#   select(t.value) %>%
#   pull()
# spn_rem_p <- spn_rem_dat %>%
#   mutate(var = row.names(.),
#          Pr...t.. = sprintf("%.3f", Pr...t..),
#          Pr...t.. = str_remove(Pr...t.., "^0+"),
#          Pr...t.. = paste("=", Pr...t..),
#          Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
#   filter(var == "blockPos_Inc:reminiscing") %>%
#   select(Pr...t..) %>%
#   pull()

# LPP as outcome
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
lpp_stm <- lmer(LPP ~ block*savoring_moment + (1|pid), data = dat)
summary(lpp_stm)
lpp_stm_dat <- data.frame(coef(summary(lpp_stm)))
lpp_stm_beta <- lpp_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
lpp_stm_std_beta <- standardize_parameters(lpp_stm) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
lpp_stm_se <- lpp_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
lpp_stm_df <- lpp_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
lpp_stm_t <- lpp_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
lpp_stm_p <- lpp_stm_dat %>%
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
front_lpp_stm <- lmer(LPP_front ~ block*savoring_moment + (1|pid), data = dat)
summary(front_lpp_stm)
front_lpp_stm_dat <- data.frame(coef(summary(front_lpp_stm)))
front_lpp_stm_beta <- front_lpp_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
front_lpp_stm_std_beta <- standardize_parameters(front_lpp_stm) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
front_lpp_stm_se <- front_lpp_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
front_lpp_stm_df <- front_lpp_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
front_lpp_stm_t <- front_lpp_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
front_lpp_stm_p <- front_lpp_stm_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# EPN as outcome
epn_stm <- lmer(EPN ~ block*savoring_moment + (1|pid), data = dat)
summary(epn_stm)
epn_stm_dat <- data.frame(coef(summary(epn_stm)))
epn_stm_beta <- epn_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
epn_stm_std_beta <- standardize_parameters(epn_stm) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
epn_stm_se <- epn_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
epn_stm_df <- epn_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
epn_stm_t <- epn_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
epn_stm_p <- epn_stm_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Pr...t..) %>%
  pull()

# n170 as outcome
n170_stm <- lmer(N170 ~ block*savoring_moment + (1|pid), data = dat)
summary(n170_stm)
n170_stm_dat <- data.frame(coef(summary(n170_stm)))
n170_stm_beta <- n170_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Estimate) %>%
  pull()
n170_stm_std_beta <- standardize_parameters(n170_stm) %>%
  filter(Parameter == "blockPos_Inc:savoring_moment") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
n170_stm_se <- n170_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(Std..Error) %>%
  pull()
n170_stm_df <- n170_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(df) %>%
  pull()
n170_stm_t <- n170_stm_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Inc:savoring_moment") %>%
  select(t.value) %>%
  pull()
n170_stm_p <- n170_stm_dat %>%
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

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_stm_moderation-data", ".RData"))
