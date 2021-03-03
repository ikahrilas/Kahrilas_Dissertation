# pairwise regulation table
library(tidyverse)
library(emmeans)
library(lmerTest)
library(r2glmm)
library(effectsize)
library(kableExtra)

# read in data and make variables for valence and regulation conditions
dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-11.csv")

per_dat_cond <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-11.csv") %>%
  separate(block, c("valence_cond", "regulation_cond"), "_") %>%
  mutate(valence_condition = if_else(valence_cond == "Neg", "Negative",
                                     if_else(valence_cond == "Pos", "Positive", "Neutral")),
         regulation_condition = if_else(regulation_cond == "Dec", "Decrease",
                                        if_else(regulation_cond == "Inc", "Increase", "Watch")))

## -- differences between increase and decrease conditions not included in table yet
# rc5
dat$block <- relevel(as.factor(dat$block), "Pos_Inc")
standardize_parameters(lmer(RC5 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Inc")
dat$block <- relevel(as.factor(dat$block), "Pos_Dec")
standardize_parameters(lmer(RC5 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Dec")

# rc11
dat$block <- relevel(as.factor(dat$block), "Pos_Inc")
standardize_parameters(lmer(RC11 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Inc")
dat$block <- relevel(as.factor(dat$block), "Pos_Dec")
standardize_parameters(lmer(RC11 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Dec")

# rc12
dat$block <- relevel(as.factor(dat$block), "Pos_Inc")
standardize_parameters(lmer(RC12 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Inc")
dat$block <- relevel(as.factor(dat$block), "Pos_Dec")
standardize_parameters(lmer(RC12 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Dec")

# pos rc12
dat$block <- relevel(as.factor(dat$block), "Pos_Inc")
standardize_parameters(lmer(pos_RC12 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Inc")
dat$block <- relevel(as.factor(dat$block), "Pos_Dec")
standardize_parameters(lmer(pos_RC12 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Dec")

# rc2
dat$block <- relevel(as.factor(dat$block), "Pos_Inc")
standardize_parameters(lmer(RC2 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Inc")
dat$block <- relevel(as.factor(dat$block), "Pos_Dec")
standardize_parameters(lmer(RC2 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Dec")

# rc3
dat$block <- relevel(as.factor(dat$block), "Pos_Inc")
standardize_parameters(lmer(RC3 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Inc")
dat$block <- relevel(as.factor(dat$block), "Pos_Dec")
standardize_parameters(lmer(RC3 ~ block + (1|pid), data = dat)) %>%
  as_tibble() %>%
  mutate(label = interpret_d(Std_Coefficient, rules = "cohen1988")) %>%
  filter(Parameter == "blockNeg_Dec")



# RC2 regulation comparisons
reg_mod_RC2 <- lmer(RC2 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
RC2_reg <- emmeans(reg_mod_RC2, pairwise ~ regulation_condition | valence_condition)

RC2_reg_confint <- data.frame(confint(RC2_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

RC2_reg_tab <- data.frame(RC2_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., RC2_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
RC2_std_beta_pos <- standardize_parameters(lmer(RC2 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(RC2 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
RC2_std_beta_neg <- standardize_parameters(lmer(RC2 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(RC2 ~ block + (1|pid), data = dat))

RC2_std_beta_pos <- RC2_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Positive")

RC2_std_beta_neg <- RC2_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Negative")

RC2_std_beta <- bind_rows(RC2_std_beta_pos, RC2_std_beta_neg)

RC2_reg_tab <- full_join(RC2_reg_tab, RC2_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "First LPP Component") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# RC3 regulation comparisons
reg_mod_RC3 <- lmer(RC3 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
RC3_reg <- emmeans(reg_mod_RC3, pairwise ~ regulation_condition | valence_condition)

RC3_reg_confint <- data.frame(confint(RC3_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

RC3_reg_tab <- data.frame(RC3_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., RC3_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
RC3_std_beta_pos <- standardize_parameters(lmer(RC3 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(RC3 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
RC3_std_beta_neg <- standardize_parameters(lmer(RC3 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(RC3 ~ block + (1|pid), data = dat))

RC3_std_beta_pos <- RC3_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Positive")

RC3_std_beta_neg <- RC3_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Negative")

RC3_std_beta <- bind_rows(RC3_std_beta_pos, RC3_std_beta_neg)

RC3_reg_tab <- full_join(RC3_reg_tab, RC3_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "Second LPP Component") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# RC5 regulation comparisons
reg_mod_RC5 <- lmer(RC5 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
RC5_reg <- emmeans(reg_mod_RC5, pairwise ~ regulation_condition | valence_condition)

RC5_reg_confint <- data.frame(confint(RC5_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

RC5_reg_tab <- data.frame(RC5_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., RC5_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
RC5_std_beta_pos <- standardize_parameters(lmer(RC5 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(RC5 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
RC5_std_beta_neg <- standardize_parameters(lmer(RC5 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(RC5 ~ block + (1|pid), data = dat))

RC5_std_beta_pos <- RC5_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Positive")

RC5_std_beta_neg <- RC5_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Negative")

RC5_std_beta <- bind_rows(RC5_std_beta_pos, RC5_std_beta_neg)

RC5_reg_tab <- full_join(RC5_reg_tab, RC5_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "P100 Component") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# RC11 regulation comparisons
reg_mod_RC11 <- lmer(RC11 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
RC11_reg <- emmeans(reg_mod_RC11, pairwise ~ regulation_condition | valence_condition)

RC11_reg_confint <- data.frame(confint(RC11_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

RC11_reg_tab <- data.frame(RC11_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., RC11_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
RC11_std_beta_pos <- standardize_parameters(lmer(RC11 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(RC11 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
RC11_std_beta_neg <- standardize_parameters(lmer(RC11 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(RC11 ~ block + (1|pid), data = dat))

RC11_std_beta_pos <- RC11_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Positive")

RC11_std_beta_neg <- RC11_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Negative")

RC11_std_beta <- bind_rows(RC11_std_beta_pos, RC11_std_beta_neg)

RC11_reg_tab <- full_join(RC11_reg_tab, RC11_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "N170 Component") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# RC12 regulation comparisons
reg_mod_RC12 <- lmer(RC12 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
RC12_reg <- emmeans(reg_mod_RC12, pairwise ~ regulation_condition | valence_condition)

RC12_reg_confint <- data.frame(confint(RC12_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

RC12_reg_tab <- data.frame(RC12_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., RC12_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
RC12_std_beta_pos <- standardize_parameters(lmer(RC12 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(RC12 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
RC12_std_beta_neg <- standardize_parameters(lmer(RC12 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(RC12 ~ block + (1|pid), data = dat))

RC12_std_beta_pos <- RC12_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Positive")

RC12_std_beta_neg <- RC12_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Negative")

RC12_std_beta <- bind_rows(RC12_std_beta_pos, RC12_std_beta_neg)

RC12_reg_tab <- full_join(RC12_reg_tab, RC12_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "EPN Component") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# pos_RC12 regulation comparisons
reg_mod_pos_RC12 <- lmer(pos_RC12 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
pos_RC12_reg <- emmeans(reg_mod_pos_RC12, pairwise ~ regulation_condition | valence_condition)

pos_RC12_reg_confint <- data.frame(confint(pos_RC12_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

pos_RC12_reg_tab <- data.frame(pos_RC12_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., pos_RC12_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
pos_RC12_std_beta_pos <- standardize_parameters(lmer(pos_RC12 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(pos_RC12 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
pos_RC12_std_beta_neg <- standardize_parameters(lmer(pos_RC12 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(pos_RC12 ~ block + (1|pid), data = dat))

pos_RC12_std_beta_pos <- pos_RC12_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Positive")

pos_RC12_std_beta_neg <- pos_RC12_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Negative")

pos_RC12_std_beta <- bind_rows(pos_RC12_std_beta_pos, pos_RC12_std_beta_neg)

pos_RC12_reg_tab <- full_join(pos_RC12_reg_tab, pos_RC12_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "EPP Component") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)


# Arousal regulation comparisons
reg_mod_ar <- lmer(arousal ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
ar_reg <- emmeans(reg_mod_ar, pairwise ~ regulation_condition | valence_condition)

ar_reg_confint <- data.frame(confint(ar_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

ar_reg_tab <- data.frame(ar_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., ar_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
ar_std_beta_pos <- standardize_parameters(lmer(arousal ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(arousal ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
ar_std_beta_neg <- standardize_parameters(lmer(arousal ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(arousal ~ block + (1|pid), data = dat))

ar_std_beta_pos <- ar_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Positive")

ar_std_beta_neg <- ar_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Negative")

ar_std_beta <- bind_rows(ar_std_beta_pos, ar_std_beta_neg)

ar_reg_tab <- full_join(ar_reg_tab, ar_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "arousal") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# Valence comparisons
reg_mod_val <- lmer(valence ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
val_reg <- emmeans(reg_mod_val, pairwise ~ regulation_condition | valence_condition)

val_reg_confint <- data.frame(confint(val_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

val_reg_tab <- data.frame(val_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., val_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
val_std_beta_pos <- standardize_parameters(lmer(valence ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(valence ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
val_std_beta_neg <- standardize_parameters(lmer(valence ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(valence ~ block + (1|pid), data = dat))

val_std_beta_pos <- val_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Positive")

val_std_beta_neg <- val_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "cohen1988")),
         valence_condition = "Negative")

val_std_beta <- bind_rows(val_std_beta_pos, ar_std_beta_neg)

val_reg_tab <- full_join(val_reg_tab, val_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "valence") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# collate components and behavioral ratings
reg_tab <- bind_rows(RC2_reg_tab, RC3_reg_tab, RC5_reg_tab, RC11_reg_tab, RC12_reg_tab, pos_RC12_reg_tab, ar_reg_tab, val_reg_tab) %>%
  mutate(Std_Coefficient = abs(Std_Coefficient))
names(reg_tab) <- c("comp", "valence_condition", "Contrast", "Estimate", "lower.CL", "upper.CL", "Std. Beta", "interpretation", "Sig.")
reg_tab <- reg_tab %>%
  mutate(Sig. = sprintf("%.3f", reg_tab$Sig.),
         Sig. = str_remove(Sig., "^0+"),
         Sig. = if_else(Sig. == ".000", "<.001", Sig.))
reg_tab <- map_df(reg_tab, ~ {
  if(class(.x) == "numeric") {
    sprintf("%.2f", .x)
  } else {
    .x
  }
})

reg_tab <- reg_tab %>%
  mutate(Estimate = if_else(Estimate == "-0.00", "0.00", Estimate),
         lower.CL = if_else(lower.CL == "-0.00", "0.00", lower.CL),
         upper.CL = if_else(upper.CL == "-0.00", "0.00", upper.CL),
         `Std. Beta` = if_else(`Std. Beta` == "-0.00", "0.00", `Std. Beta`),
         CI = paste0(lower.CL, ", ", upper.CL)) %>%
  select(comp:Estimate, "CI", everything()) %>%
  select(-lower.CL, -upper.CL) %>%
  mutate("Estimate (95\\% CI)" = paste0(Estimate, " (", CI, ")"),
         "Std. Beta (Label)" = paste0(`Std. Beta`, " (", interpretation, ")")) %>%
  select(comp, valence_condition, Contrast, "Estimate (95\\% CI)", "Std. Beta (Label)", Sig.)

reg_tab_pos <- reg_tab %>%
  filter(valence_condition == "Positive") %>%
  select(-valence_condition)

reg_tab_neg <- reg_tab %>%
  filter(valence_condition == "Negative") %>%
  select(`Estimate (95\\% CI)`:Sig.)

reg_tab_wide <- bind_cols(reg_tab_pos, reg_tab_neg)

names(reg_tab_wide) <- c(names(reg_tab_pos), names(reg_tab_pos)[3:5])

reg_tab_wide[-1] %>%
  kable(., escape = FALSE, booktabs = TRUE, align = c("l", rep("c", times = 6)), linesep = "", caption = "(ref:pairwise-reg-comparison-table)") %>%
  kable_styling(latex_options = "scale_down") %>%
  add_header_above(c(" ", "Positive Images" = 3, "Negative Images" = 3), bold = TRUE, italic = TRUE) %>%
  row_spec(0, align = "c") %>%
  pack_rows("Early LPP Component", 1, 3) %>%
  pack_rows("Late LPP Component", 4, 6) %>%
  pack_rows("P125 Component", 7, 9) %>%
  pack_rows("N170 Component", 10, 12) %>%
  pack_rows("EPN Component", 13, 15) %>%
  pack_rows("EPP Component", 16, 18) %>%
  pack_rows("Arousal Ratings", 19, 21) %>%
  pack_rows("Valence Ratings", 22, 24) %>%
  footnote(escape = FALSE,
           general_title = "Note.",
           general = "Std. Beta = Absolute value of standardized beta coefficient as measure of effect size derived by fitting model to standardized dataset with effect size label as per Cohen's (1988) recommendations, Sig. = $p$ value. $P$ values and confidence intervals adjusted using the Tukey method for comparing a family of three estimates.",
           threeparttable = TRUE,
           footnote_as_chunk = TRUE)

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "factor_score_regulations_table-data", ".RData"))
