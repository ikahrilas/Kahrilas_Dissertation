# pa mlm analyses

library(tidyverse)
library(lmerTest)
library(effectsize)

dat <- read_csv("data/paper_two/created_data/per_data_analyses_2020_5_19.csv") %>%
  mutate(pos_affectivity = scale(pos_affectivity, scale = FALSE, center = TRUE)) %>%
  select(pid:LPP, N170, EPN, everything())

# # SPN as outcome
# dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
# spn_pa <- lmer(SPN ~ block*pos_affectivity + (1|pid), data = dat)
# summary(spn_pa)
# spn_pa_dat <- data.frame(coef(summary(spn_pa)))
# spn_pa_beta <- spn_pa_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(Estimate = sprintf("%.2f", Estimate)) %>%
#   filter(var == "blockPos_Watch:pos_affectivity") %>%
#   select(Estimate) %>%
#   pull()
# spn_pa_std_beta <- standardize_parameters(spn_pa) %>%
#   filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
#   mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
#   select(`Std_Coefficient`) %>%
#   pull()
# spn_pa_se <- spn_pa_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
#   filter(var == "blockPos_Watch:pos_affectivity") %>%
#   select(Std..Error) %>%
#   pull()
# spn_pa_df <- spn_pa_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(df = sprintf("%.2f", df)) %>%
#   filter(var == "blockPos_Watch:pos_affectivity") %>%
#   select(df) %>%
#   pull()
# spn_pa_t <- spn_pa_dat %>%
#   mutate(var = row.names(.)) %>%
#   mutate(t.value = sprintf("%.2f", t.value)) %>%
#   filter(var == "blockPos_Watch:pos_affectivity") %>%
#   select(t.value) %>%
#   pull()
# spn_pa_p <- spn_pa_dat %>%
#   mutate(var = row.names(.),
#          Pr...t.. = sprintf("%.3f", Pr...t..),
#          Pr...t.. = str_remove(Pr...t.., "^0+"),
#          Pr...t.. = paste("=", Pr...t..),
#          Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
#   filter(var == "blockPos_Watch:pos_affectivity") %>%
#   select(Pr...t..) %>%
#   pull()

# LPP as outcome
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
lpp_pa <- lmer(LPP ~ block*pos_affectivity + (1|pid), data = dat)
summary(lpp_pa)
lpp_pa_dat <- data.frame(coef(summary(lpp_pa)))
lpp_pa_beta <- lpp_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
lpp_pa_std_beta <- standardize_parameters(lpp_pa) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
lpp_pa_se <- lpp_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
lpp_pa_df <- lpp_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
lpp_pa_t <- lpp_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
lpp_pa_p <- lpp_pa_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# frontal LPP as outcome
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
front_lpp_pa <- lmer(LPP_front ~ block*pos_affectivity + (1|pid), data = dat)
summary(front_lpp_pa)
front_lpp_pa_dat <- data.frame(coef(summary(front_lpp_pa)))
front_lpp_pa_beta <- front_lpp_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
front_lpp_pa_std_beta <- standardize_parameters(front_lpp_pa) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
front_lpp_pa_se <- front_lpp_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
front_lpp_pa_df <- front_lpp_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
front_lpp_pa_t <- front_lpp_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
front_lpp_pa_p <- front_lpp_pa_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# EPN as outcome
epn_pa <- lmer(EPN ~ block*pos_affectivity + (1|pid), data = dat)
summary(epn_pa)
epn_pa_dat <- data.frame(coef(summary(epn_pa)))
epn_pa_beta <- epn_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
epn_pa_std_beta <- standardize_parameters(epn_pa) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
epn_pa_se <- epn_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
epn_pa_df <- epn_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
epn_pa_t <- epn_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
epn_pa_p <- epn_pa_dat %>%
  mutate(var = row.names(.),
         Pr...t.. = sprintf("%.3f", Pr...t..),
         Pr...t.. = str_remove(Pr...t.., "^0+"),
         Pr...t.. = paste("=", Pr...t..),
         Pr...t.. = if_else(Pr...t.. == "= .000", "< .001", Pr...t..)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Pr...t..) %>%
  pull()

# n170 as outcome
n170_pa <- lmer(N170 ~ block*pos_affectivity + (1|pid), data = dat)
summary(n170_pa)
n170_pa_dat <- data.frame(coef(summary(n170_pa)))
n170_pa_beta <- n170_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Estimate = sprintf("%.2f", Estimate)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Estimate) %>%
  pull()
n170_pa_std_beta <- standardize_parameters(n170_pa) %>%
  filter(Parameter == "blockPos_Watch:pos_affectivity") %>%
  mutate(Std_Coefficient = sprintf("%.2f", Std_Coefficient)) %>%
  select(`Std_Coefficient`) %>%
  pull()
n170_pa_se <- n170_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(Std..Error = sprintf("%.2f", Std..Error)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(Std..Error) %>%
  pull()
n170_pa_df <- n170_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(df = sprintf("%.2f", df)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(df) %>%
  pull()
n170_pa_t <- n170_pa_dat %>%
  mutate(var = row.names(.)) %>%
  mutate(t.value = sprintf("%.2f", t.value)) %>%
  filter(var == "blockPos_Watch:pos_affectivity") %>%
  select(t.value) %>%
  pull()
n170_pa_p <- n170_pa_dat %>%
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

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_pa_moderation-data", ".RData"))
