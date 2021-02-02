# Pairwise watch comparison table

library(tidyverse)
library(emmeans)
library(lmerTest)
library(r2glmm)
library(effectsize)
library(kableExtra)

# read in data and make variables for valence and regulation conditions
dat <- read_csv("data/paper_two/created_data/per_data_analyses_2020_5_19.csv")

per_dat_cond <- read_csv("data/paper_two/created_data/per_data_analyses_2020_5_19.csv") %>%
  select(pid:LPP, N170, EPN, everything()) %>%
  separate(block, c("valence_cond", "regulation_cond"), "_") %>%
  mutate(valence_condition = if_else(valence_cond == "Neg", "Negative",
                                     if_else(valence_cond == "Pos", "Positive", "Neutral")),
         regulation_condition = if_else(regulation_cond == "Dec", "Decrease",
                                        if_else(regulation_cond == "Inc", "Increase", "Watch")))

# LPP watch comparisons
reg_mod_lpp <- lmer(LPP ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
lpp_watch <- emmeans(reg_mod_lpp, pairwise ~ valence_condition | regulation_condition)
lpp_watch_confint <- data.frame(confint(lpp_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
lpp_watch_tab <- data.frame(lpp_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., lpp_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
lpp_std_beta <- standardize_parameters(lmer(LPP ~ block + (1|pid), data = dat))
lpp_std_beta <- lpp_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(LPP ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
lpp_std_beta <- bind_rows(lpp_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
lpp_watch_tab <- full_join(lpp_watch_tab, lpp_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
lpp_watch_tab <- lpp_watch_tab %>%
  mutate(comp = "LPP") %>%
  select(comp, everything())

# frontal LPP watch comparisons
reg_mod_lpp_front <- lmer(LPP_front ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
lpp_front_watch <- emmeans(reg_mod_lpp_front, pairwise ~ valence_condition | regulation_condition)
lpp_front_watch_confint <- data.frame(confint(lpp_front_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
lpp_front_watch_tab <- data.frame(lpp_front_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., lpp_front_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
lpp_front_std_beta <- standardize_parameters(lmer(LPP_front ~ block + (1|pid), data = dat))
lpp_front_std_beta <- lpp_front_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(LPP_front ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
lpp_front_std_beta <- bind_rows(lpp_front_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
lpp_front_watch_tab <- full_join(lpp_front_watch_tab, lpp_front_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
lpp_front_watch_tab <- lpp_front_watch_tab %>%
  mutate(comp = "Frontal LPP") %>%
  select(comp, everything())

# EPN comparisons
reg_mod_epn <- lmer(EPN ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
epn_watch <- emmeans(reg_mod_epn, pairwise ~ valence_condition | regulation_condition)
epn_watch_confint <- data.frame(confint(epn_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
epn_watch_tab <- data.frame(epn_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., epn_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
epn_std_beta <- standardize_parameters(lmer(EPN ~ block + (1|pid), data = dat))
epn_std_beta <- epn_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(EPN ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
epn_std_beta <- bind_rows(epn_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
epn_watch_tab <- full_join(epn_watch_tab, epn_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
epn_watch_tab <- epn_watch_tab %>%
  mutate(comp = "EPN") %>%
  select(comp, everything())

# N170 comparisons
reg_mod_n170 <- lmer(N170 ~ valence_condition * regulation_condition + (1|pid), data = per_dat_cond)
n170_watch <- emmeans(reg_mod_n170, pairwise ~ valence_condition | regulation_condition)
n170_watch_confint <- data.frame(confint(n170_watch)$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  select(contrast, lower.CL, upper.CL)
n170_watch_tab <- data.frame(n170_watch$contrasts) %>%
  filter(regulation_condition == "Watch") %>%
  left_join(., n170_watch_confint, by = "contrast")
# derive standardized beta for effect size
dat$block <- relevel(factor(dat$block), ref = "Neu_Watch")
n170_std_beta <- standardize_parameters(lmer(N170 ~ block + (1|pid), data = dat))
n170_std_beta <- n170_std_beta %>%
  filter(Parameter %in% c("blockNeg_Watch", "blockPos_Watch")) %>%
  mutate(Parameter = if_else(Parameter == "blockNeg_Watch", "Negative - Neutral", "Neutral - Positive"))
dat$block <- relevel(factor(dat$block), ref = "Pos_Watch")
tmp <- standardize_parameters(lmer(N170 ~ block + (1|pid), data = dat))
tmp <- tmp %>%
  filter(Parameter == "blockNeg_Watch") %>%
  mutate(Parameter = "Negative - Positive")
n170_std_beta <- bind_rows(n170_std_beta, tmp) %>%
  mutate(interpretation = tools::toTitleCase(interpret_d(Std_Coefficient, rules = "gignac2016")))
n170_watch_tab <- full_join(n170_watch_tab, n170_std_beta, by = c("contrast" = "Parameter")) %>%
  select(contrast, estimate, lower.CL, upper.CL, Std_Coefficient, interpretation, t.ratio, p.value)
n170_watch_tab <- n170_watch_tab %>%
  mutate(comp = "N170") %>%
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

watch_tab <- bind_rows(n170_watch_tab, epn_watch_tab, lpp_watch_tab, lpp_front_watch_tab, ar_watch_tab, val_watch_tab) %>%
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
  pack_rows("N170", 1, 3) %>%
  pack_rows("EPN", 4, 6) %>%
  pack_rows("LPP", 7, 9) %>%
  pack_rows("Frontal LPP", 10, 12) %>%
  pack_rows("Arousal Ratings", 13, 15) %>%
  pack_rows("Valence Ratings", 16, 18) %>%
  footnote(escape = FALSE,
           footnote_as_chunk = TRUE,
           general_title = "Note.",
           general = "Std. Beta (Label) = Absolute value of standardized
beta coefficient as measure of effect size derived by fitting model to s
tandardized dataset with effect size label as per Funder's (2019) recommendations,
Sig. = $p$ value. $P$ values and confidence intervals adjusted using the Tukey method
for comparing a family of three estimates.",
           threeparttable = TRUE)

 save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_pairwise_watch_table-data", ".RData"))
