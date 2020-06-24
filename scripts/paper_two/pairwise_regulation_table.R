# pairwise regulation table
library(tidyverse)
library(emmeans)
library(lmerTest)
library(r2glmm)
library(effectsize)


# read in data and make variables for valence and regulation conditions
dat <- read_csv("data/paper_two/created_data/per_data_analyses_2020_5_19.csv") %>%
  select(pid:LPP, N170, EPN, everything())

per_dat_cond <- read_csv("data/paper_two/created_data/per_data_analyses_2020_5_19.csv") %>%
  select(pid:LPP, N170, EPN, everything()) %>%
  separate(block, c("valence_cond", "regulation_cond"), "_") %>%
  mutate(valence_condition = if_else(valence_cond == "Neg", "Negative",
                                     if_else(valence_cond == "Pos", "Positive", "Neutral")),
         regulation_condition = if_else(regulation_cond == "Dec", "Decrease",
                                        if_else(regulation_cond == "Inc", "Increase", "Watch")))

# # SPN regulation comparisons
# reg_mod_spn <- lmer(SPN ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
# spn_reg <- emmeans(reg_mod_spn, pairwise ~ regulation_condition | valence_condition)
#
# spn_reg_confint <- data.frame(confint(spn_reg)$contrasts) %>%
#   filter(valence_condition %in% c("Negative", "Positive")) %>%
#   select(contrast, valence_condition, lower.CL, upper.CL)
#
# spn_reg_tab <- data.frame(spn_reg$contrasts) %>%
#   filter(valence_condition %in% c("Negative", "Positive")) %>%
#   left_join(., spn_reg_confint, by = c("contrast", "valence_condition"))
#
# # derive standardized beta for effect size
# dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
# spn_std_beta_pos <- standardize_parameters(lmer(SPN ~ block + (1|pid), data = dat))
# dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
# tmp_pos <- standardize_parameters(lmer(SPN ~ block + (1|pid), data = dat))
# dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
# spn_std_beta_neg <- standardize_parameters(lmer(SPN ~ block + (1|pid), data = dat))
# dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
# tmp_neg <- standardize_parameters(lmer(SPN ~ block + (1|pid), data = dat))
#
# spn_std_beta_pos <- spn_std_beta_pos %>%
#   filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
#   mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
#   bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
#   mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
#          interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
#          valence_condition = "Positive")
#
# spn_std_beta_neg <- spn_std_beta_neg %>%
#   filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
#   mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
#   bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
#   mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
#          interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
#          valence_condition = "Negative")
#
# spn_std_beta <- bind_rows(spn_std_beta_pos, spn_std_beta_neg)
#
# spn_reg_tab <- full_join(spn_reg_tab, spn_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
#   mutate(comp = "SPN") %>%
#   select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# Frontal LPP regulation comparisons
reg_mod_lpp_front <- lmer(LPP_front ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
lpp_front_reg <- emmeans(reg_mod_lpp_front, pairwise ~ regulation_condition | valence_condition)

lpp_front_reg_confint <- data.frame(confint(lpp_front_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

lpp_front_reg_tab <- data.frame(lpp_front_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., lpp_front_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
lpp_front_std_beta_pos <- standardize_parameters(lmer(LPP_front ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(LPP_front ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
lpp_front_std_beta_neg <- standardize_parameters(lmer(LPP_front ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(LPP_front ~ block + (1|pid), data = dat))

lpp_front_std_beta_pos <- lpp_front_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Positive")

lpp_front_std_beta_neg <- lpp_front_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Negative")

lpp_front_std_beta <- bind_rows(lpp_front_std_beta_pos, lpp_front_std_beta_neg)

lpp_front_reg_tab <- full_join(lpp_front_reg_tab, lpp_front_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "Frontal LPP") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# LPP regulation comparisons
reg_mod_lpp <- lmer(LPP ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
lpp_reg <- emmeans(reg_mod_lpp, pairwise ~ regulation_condition | valence_condition)

lpp_reg_confint <- data.frame(confint(lpp_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

lpp_reg_tab <- data.frame(lpp_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., lpp_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
lpp_std_beta_pos <- standardize_parameters(lmer(LPP ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(LPP ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
lpp_std_beta_neg <- standardize_parameters(lmer(LPP ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(LPP ~ block + (1|pid), data = dat))

lpp_std_beta_pos <- lpp_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Positive")

lpp_std_beta_neg <- lpp_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Negative")

lpp_std_beta <- bind_rows(lpp_std_beta_pos, lpp_std_beta_neg)

lpp_reg_tab <- full_join(lpp_reg_tab, lpp_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "LPP") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# EPN regulation comparisons
reg_mod_epn <- lmer(EPN ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
epn_reg <- emmeans(reg_mod_epn, pairwise ~ regulation_condition | valence_condition)

epn_reg_confint <- data.frame(confint(epn_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

epn_reg_tab <- data.frame(epn_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., epn_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
epn_std_beta_pos <- standardize_parameters(lmer(EPN ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(EPN ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
epn_std_beta_neg <- standardize_parameters(lmer(EPN ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(EPN ~ block + (1|pid), data = dat))

epn_std_beta_pos <- epn_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Positive")

epn_std_beta_neg <- epn_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Negative")

epn_std_beta <- bind_rows(epn_std_beta_pos, epn_std_beta_neg)

epn_reg_tab <- full_join(epn_reg_tab, epn_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "EPN") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# N170 regulation comparisons
reg_mod_n170 <- lmer(N170 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
n170_reg <- emmeans(reg_mod_n170, pairwise ~ regulation_condition | valence_condition)

n170_reg_confint <- data.frame(confint(n170_reg)$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  select(contrast, valence_condition, lower.CL, upper.CL)

n170_reg_tab <- data.frame(n170_reg$contrasts) %>%
  filter(valence_condition %in% c("Negative", "Positive")) %>%
  left_join(., n170_reg_confint, by = c("contrast", "valence_condition"))

# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
n170_std_beta_pos <- standardize_parameters(lmer(N170 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Pos_Inc")
tmp_pos <- standardize_parameters(lmer(N170 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Watch")
n170_std_beta_neg <- standardize_parameters(lmer(N170 ~ block + (1|pid), data = dat))
dat$block <- relevel(factor(dat$block), ref = "Neg_Inc")
tmp_neg <- standardize_parameters(lmer(N170 ~ block + (1|pid), data = dat))

n170_std_beta_pos <- n170_std_beta_pos %>%
  filter(Parameter %in% c("blockPos_Dec", "blockPos_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_pos, Parameter == "blockPos_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockPos_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Positive")

n170_std_beta_neg <- n170_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Negative")

n170_std_beta <- bind_rows(n170_std_beta_pos, n170_std_beta_neg)

n170_reg_tab <- full_join(n170_reg_tab, n170_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "N170") %>%
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
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Positive")

ar_std_beta_neg <- ar_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
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
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Positive")

val_std_beta_neg <- val_std_beta_neg %>%
  filter(Parameter %in% c("blockNeg_Dec", "blockNeg_Inc")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Watch", "Increase - Watch")) %>%
  bind_rows(., filter(tmp_neg, Parameter == "blockNeg_Dec")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Dec", "Decrease - Increase", Parameter),
         interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "funder2019")),
         valence_condition = "Negative")

val_std_beta <- bind_rows(val_std_beta_pos, ar_std_beta_neg)

val_reg_tab <- full_join(val_reg_tab, val_std_beta, by = c("contrast" = "Parameter", "valence_condition")) %>%
  mutate(comp = "valence") %>%
  select(comp, valence_condition, contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, p.value)

# collate components and behavioral ratings
reg_tab <- bind_rows(n170_reg_tab, epn_reg_tab, lpp_reg_tab, lpp_front_reg_tab, ar_reg_tab, val_reg_tab) %>%
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

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_regulations_table-data", ".RData"))

# reg_tab_wide[-1] %>%
#   kable(., escape = FALSE, booktabs = TRUE, align = c("l", rep("c", times = 6)), linesep = "", caption = "(ref:pairwise-reg-comparison-table)") %>%
#   kable_styling(latex_options = "scale_down") %>%
#   add_header_above(c(" ", "Positive Images" = 3, "Negative Images" = 3), bold = TRUE, italic = TRUE) %>%
#   row_spec(0, align = "c") %>%
#   pack_rows("N170", 1, 3) %>%
#   pack_rows("EPN", 4, 6) %>%
#   pack_rows("LPP", 7, 9) %>%
#   pack_rows("Frontal LPP", 10, 12) %>%
#   pack_rows("Arousal Ratings", 13, 15) %>%
#   pack_rows("Valence Ratings", 16, 18) %>%
#   footnote(escape = FALSE,
#            general_title = "Note.",
#            general = "Std. Beta = Absolute value of standardized beta coefficient as measure of effect size derived by fitting model to standardized dataset with effect size label as per Funder's (2019) recommendations, Sig. = $p$ value. $P$ values and confidence intervals adjusted using the Tukey method for comparing a family of three estimates.",
#            threeparttable = TRUE,
#            footnote_as_chunk = TRUE)
