# block analyses for factor scores

# load pacakges
library(tidyverse)
library(lmerTest)
library(emmeans)
library(here)
library(performance)

# read in data
dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-02.csv")

# RC2 analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
RC2_block_mod <- lmer(RC2 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC2_block_mod)
## check results
summary(RC2_block_mod)
## all contrasts
emmeans(RC2_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "none")
### significant differences in pos/neg watch and neutral watch
### no regulation effects for negative stimuli
### significant difference between positive watch and positive inc, no positive decrease

RC2_block_aov <- anova(RC2_block_mod)
df_num_rc2 <- RC2_block_aov$NumDF
df_den_rc2 <- sprintf("%.2f", RC2_block_aov$DenDF)
f_rc2 <- sprintf("%.2f", RC2_block_aov[["F value"]])


# RC3 analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
RC3_block_mod <- lmer(RC3 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC3_block_mod)
## check results
summary(RC3_block_mod)
## all contrasts
emmeans(RC3_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "none")
### significant differences in pos/neg watch and neutral watch
### no regulation effects for negative stimuli, but trending for negative decrease
### trending difference between positive watch and positive inc, no positive decrease

RC3_block_aov <- anova(RC3_block_mod)
df_num_rc3 <- RC3_block_aov$NumDF
df_den_rc3 <- sprintf("%.2f", RC3_block_aov$DenDF)
f_rc3 <- sprintf("%.2f", RC2_block_aov[["F value"]])


# RC5 analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
RC5_block_mod <- lmer(RC5 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC5_block_mod)
## check results
summary(RC5_block_mod)
## all contrasts
emmeans(RC5_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "none")
### total dud - nothing interesting

RC5_block_aov <- anova(RC5_block_mod)
df_num_rc5 <- RC5_block_aov$NumDF
df_den_rc5 <- sprintf("%.2f", RC5_block_aov$DenDF)
f_rc5 <- sprintf("%.2f", RC5_block_aov[["F value"]])


# RC11 analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
RC11_block_mod <- lmer(RC11 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC11_block_mod)
## check results
summary(RC11_block_mod)
## all contrasts
emmeans(RC11_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "none")
### arousal effects
### no negative regulatory effects
### no positive regulatory effects

RC11_block_aov <- anova(RC11_block_mod)
df_num_rc11 <- RC11_block_aov$NumDF
df_den_rc11 <- sprintf("%.2f", RC11_block_aov$DenDF)
f_rc11 <- sprintf("%.2f", RC11_block_aov[["F value"]])


# RC12 analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
RC12_block_mod <- lmer(RC12 ~ block + (1|pid), data = dat)
## check assumptions
check_model(RC12_block_mod)
## check results
summary(RC12_block_mod)
## all contrasts
emmeans(RC12_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "none")
### arousal effect for negative stimuli, not so with positive
### negative increase effect
### positive increase effect (!)

RC12_block_aov <- anova(RC12_block_mod)
df_num_rc12 <- RC12_block_aov$NumDF
df_den_rc12 <- sprintf("%.2f", RC12_block_aov$DenDF)
f_rc12 <- sprintf("%.2f", RC12_block_aov[["F value"]])


# arousal analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
arousal_block_mod <- lmer(arousal ~ block + (1|pid), data = dat)
## check assumptions
check_model(arousal_block_mod)
## check results
summary(arousal_block_mod)
## all contrasts
emmeans(arousal_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "none")
### arousal effect
### negative increase effect
### positive increase effect

arousal_block_aov <- anova(arousal_block_mod)
df_num_arousal <- arousal_block_aov$NumDF
df_den_arousal <- sprintf("%.2f", arousal_block_aov$DenDF)
f_arousal <- sprintf("%.2f", arousal_block_aov[["F value"]])


# valence analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
valence_block_mod <- lmer(valence ~ block + (1|pid), data = dat)
## check assumptions
check_model(valence_block_mod)
## check results
summary(valence_block_mod)
## all contrasts
emmeans(valence_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "none")
### arousal effect
### negative increase effect
### positive regulatory effects

valence_block_aov <- anova(valence_block_mod)
df_num_valence <- valence_block_aov$NumDF
df_den_valence <- sprintf("%.2f", valence_block_aov$DenDF)
f_valence <- sprintf("%.2f", valence_block_aov[["F value"]])

# difficulty analyses
## relevel block factor so that neutral watch is the reference
dat$block <- relevel(as.factor(dat$block), ref = "Neu_Watch")
difficulty_block_mod <- lmer(difficulty ~ block + (1|pid), data = dat)
## check assumptions
check_model(difficulty_block_mod)
## check results
summary(difficulty_block_mod)
## all contrasts
emmeans(difficulty_block_mod, data = dat, ~block) %>%
  contrast("pairwise", adjust = "none")
### no arousal differences
### negative regulatory differences
### positive regulatory differences

difficulty_block_aov <- anova(difficulty_block_mod)
df_num_difficulty <- difficulty_block_aov$NumDF
df_den_difficulty <- sprintf("%.2f", difficulty_block_aov$DenDF)
f_difficulty <- sprintf("%.2f", difficulty_block_aov[["F value"]])

save.image(file = paste0("data/paper_two/analyses/", Sys.Date(), "_repeated_measures_ANOVA_analysis-data", ".RData"))
