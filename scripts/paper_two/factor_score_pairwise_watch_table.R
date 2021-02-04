# Pairwise watch comparison table for factor score data

# load packages
library(tidyverse)
library(emmeans)
library(lmerTest)
library(r2glmm)
library(effectsize)
library(kableExtra)

# read in data and make variables for valence and regulation conditions
dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-02.csv")

per_dat_cond <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-02.csv") %>%
  separate(block, c("valence_cond", "regulation_cond"), "_") %>%
  mutate(valence_condition = if_else(valence_cond == "Neg", "Negative",
                                     if_else(valence_cond == "Pos", "Positive", "Neutral")),
         regulation_condition = if_else(regulation_cond == "Dec", "Decrease",
                                        if_else(regulation_cond == "Inc", "Increase", "Watch")))

# RC2 watch comparisons
reg_mod_rc2 <- lmer(RC2 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
rc2_watch <- emmeans(reg_mod_rc2, pairwise ~ valence_condition | regulation_condition)
rc2_watch_confint <- data.frame(confint(rc2_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
rc2_watch_tab <- data.frame(rc2_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., rc2_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
rc2_std_beta <- standardize_parameters(lmer(RC2 ~ block + (1|pid), data = dat))
rc2_std_beta <- rc2_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(RC2 ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
rc2_std_beta <- bind_rows(rc2_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
rc2_watch_tab <- full_join(rc2_watch_tab, rc2_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
rc2_watch_tab <- rc2_watch_tab %>%
  mutate(comp = "First LPP Component") %>%
  select(comp, everything())

# RC3 watch comparisons
reg_mod_rc3 <- lmer(RC3 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
rc3_watch <- emmeans(reg_mod_rc3, pairwise ~ valence_condition | regulation_condition)
rc3_watch_confint <- data.frame(confint(rc3_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
rc3_watch_tab <- data.frame(rc3_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., rc3_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
rc3_std_beta <- standardize_parameters(lmer(RC3 ~ block + (1|pid), data = dat))
rc3_std_beta <- rc3_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(RC3 ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
rc3_std_beta <- bind_rows(rc3_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
rc3_watch_tab <- full_join(rc3_watch_tab, rc3_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
rc3_watch_tab <- rc3_watch_tab %>%
  mutate(comp = "Second LPP Component") %>%
  select(comp, everything())

# RC5 watch comparisons
reg_mod_RC5 <- lmer(RC5 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
RC5_watch <- emmeans(reg_mod_RC5, pairwise ~ valence_condition | regulation_condition)
RC5_watch_confint <- data.frame(confint(RC5_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
RC5_watch_tab <- data.frame(RC5_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., RC5_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
RC5_std_beta <- standardize_parameters(lmer(RC5 ~ block + (1|pid), data = dat))
RC5_std_beta <- RC5_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(RC5 ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
RC5_std_beta <- bind_rows(RC5_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
RC5_watch_tab <- full_join(RC5_watch_tab, RC5_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
RC5_watch_tab <- RC5_watch_tab %>%
  mutate(comp = "P100 Component") %>%
  select(comp, everything())

# RC11 watch comparisons
reg_mod_RC11 <- lmer(RC11 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
RC11_watch <- emmeans(reg_mod_RC11, pairwise ~ valence_condition | regulation_condition)
RC11_watch_confint <- data.frame(confint(RC11_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
RC11_watch_tab <- data.frame(RC11_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., RC11_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
RC11_std_beta <- standardize_parameters(lmer(RC11 ~ block + (1|pid), data = dat))
RC11_std_beta <- RC11_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(RC11 ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
RC11_std_beta <- bind_rows(RC11_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
RC11_watch_tab <- full_join(RC11_watch_tab, RC11_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
RC11_watch_tab <- RC11_watch_tab %>%
  mutate(comp = "N170 Component") %>%
  select(comp, everything())

# RC12 watch comparisons
reg_mod_RC12 <- lmer(RC12 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
RC12_watch <- emmeans(reg_mod_RC12, pairwise ~ valence_condition | regulation_condition)
RC12_watch_confint <- data.frame(confint(RC12_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
RC12_watch_tab <- data.frame(RC12_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., RC12_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
RC12_std_beta <- standardize_parameters(lmer(RC12 ~ block + (1|pid), data = dat))
RC12_std_beta <- RC12_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(RC12 ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
RC12_std_beta <- bind_rows(RC12_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
RC12_watch_tab <- full_join(RC12_watch_tab, RC12_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
RC12_watch_tab <- RC12_watch_tab %>%
  mutate(comp = "EPN Component") %>%
  select(comp, everything())

# Arousal comparisons
reg_mod_ar <- lmer(arousal ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
ar_watch <- emmeans(reg_mod_ar, pairwise ~ valence_condition | regulation_condition)
ar_watch_confint <- data.frame(confint(ar_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
ar_watch_tab <- data.frame(ar_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., ar_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
ar_std_beta <- standardize_parameters(lmer(arousal ~ block + (1|pid), data = dat))
ar_std_beta <- ar_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(arousal ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
ar_std_beta <- bind_rows(ar_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
ar_watch_tab <- full_join(ar_watch_tab, ar_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
ar_watch_tab <- ar_watch_tab %>%
  mutate(comp = "arousal_rating") %>%
  select(comp, everything())

# Valence comparisons
reg_mod_val <- lmer(valence ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
val_watch <- emmeans(reg_mod_val, pairwise ~ valence_condition | regulation_condition)
val_watch_confint <- data.frame(confint(val_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
val_watch_tab <- data.frame(val_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., ar_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
val_std_beta <- standardize_parameters(lmer(valence ~ block + (1|pid), data = dat))
val_std_beta <- val_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(valence ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
val_std_beta <- bind_rows(val_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
val_watch_tab <- full_join(val_watch_tab, val_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
val_watch_tab <- val_watch_tab %>%
  mutate(comp = "valence_rating") %>%
  select(comp, everything())

watch_tab <- bind_rows(rc2_watch_tab, rc3_watch_tab, RC5_watch_tab, RC11_watch_tab, RC12_watch_tab, ar_watch_tab, val_watch_tab) %>%
  mutate(Std_Coefficient = abs(Std_Coefficient))
names(watch_tab) <- c("comp", "Contrast", "Estimate", "lower.CL", "upper.CL", "Std. Beta", "interpretation", "$t$ ratio", "Sig.")
watch_tab <- watch_tab %>%
  mutate(Sig. = sprintf("%.3f", watch_tab$Sig.),
         Sig. = str_remove(Sig., "^0+"),
         Sig. = if_else(Sig. == ".000", "<.001", Sig.))
watch_tab <- map_df(watch_tab, ~ {
  if(class(.x) == "numeric") {
    sprintf("%.2f", .x)
  } else {
    .x
  }
})
watch_tab <- watch_tab %>%
  mutate(Estimate = if_else(Estimate == "-0.00", "0.00", Estimate),
         lower.CL = if_else(lower.CL == "-0.00", "0.00", lower.CL),
         upper.CL = if_else(upper.CL == "-0.00", "0.00", upper.CL),
         `Std. Beta` = if_else(`Std. Beta` == "-0.00", "0.00", `Std. Beta`),
         CI = paste0(lower.CL, ", ", upper.CL)) %>%
  select(comp:Estimate, "CI", everything()) %>%
  select(-lower.CL, -upper.CL, -`$t$ ratio`) %>%
  mutate("Estimate (95\\% CI)" = paste0(watch_tab$Estimate, " (", CI, ")"),
         "Std. Beta (Label)" = paste0(watch_tab$`Std. Beta`, " (", interpretation, ")")) %>%
  select(comp, Contrast, "Estimate (95\\% CI)", "Std. Beta (Label)", Sig.)

watch_tab %>%
  select(-comp) %>%
  kable(., escape = FALSE, booktabs = TRUE, align = c("l", "c", "c", "c"), linesep = "", caption = "(ref:pairwise-watch-comparison-table)") %>%
  row_spec(0, align = "c") %>%
  pack_rows("First LPP Component", 1, 3) %>%
  pack_rows("Second LPP Component", 4, 6) %>%
  pack_rows("P100 Component", 7, 9) %>%
  pack_rows("N170 Component", 10, 12) %>%
  pack_rows("EPN Component", 13, 15) %>%
  pack_rows("Arousal Ratings", 16, 18) %>%
  pack_rows("Valence Ratings", 19, 21) %>%
  footnote(escape = FALSE,
           footnote_as_chunk = TRUE,
           general_title = "Note.",
           general = "Std. Beta (Label) = Absolute value of standardized
beta coefficient as measure of effect size derived by fitting model to s
tandardized dataset with effect size label as per Gignac's (2016) recommendations,
Sig. = $p$ value. $P$ values and confidence intervals adjusted using the Tukey method
for comparing a family of three estimates.",
           threeparttable = TRUE)

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "factor_score_pairwise_watch_table-data", ".RData"))
