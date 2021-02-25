# Pairwise comparisons for positive/negative images within increase/decrease conditions

# load packages
library(tidyverse)
library(emmeans)
library(lmerTest)
library(r2glmm)
library(effectsize)
library(kableExtra)

dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-11.csv")

per_dat_cond <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-11.csv") %>%
  separate(block, c("valence_cond", "regulation_cond"), "_") %>%
  mutate(valence_condition = if_else(valence_cond == "Neg", "Negative",
                                     if_else(valence_cond == "Pos", "Positive", "Neutral")),
         regulation_condition = if_else(regulation_cond == "Dec", "Decrease",
                                        if_else(regulation_cond == "Inc", "Increase", "Watch")))

# change block to factor type
dat$block <- as.factor(dat$block)
dat$block <- relevel(dat$block, "Neg_Inc")

# run analyses of increase conditions
## RC2
inc_mod_rc2 <- lmer(RC2 ~ block + (1|pid), data = dat)
summary(inc_mod_rc2)
anova(inc_mod_rc2)
### no differences

## RC3
inc_mod_rc3 <- lmer(RC3 ~ block + (1|pid), data = dat)
summary(inc_mod_rc3)
### sig difference, negative increase greater than positive increase

## RC5
inc_mod_rc5 <- lmer(RC5 ~ block + (1|pid), data = dat)
summary(inc_mod_rc5)
### no sig difference

## RC11
inc_mod_rc11 <- lmer(RC11 ~ block + (1|pid), data = dat)
summary(inc_mod_rc11)
### no sig difference

## RC12
inc_mod_rc12 <- lmer(RC12 ~ block + (1|pid), data = dat)
summary(inc_mod_rc12)
### sig difference, positive increase greater than negative increase

# run analyses of decrease conditions
dat$block <- relevel(dat$block, "Neg_Dec")

## RC2
dec_mod_rc2 <- lmer(RC2 ~ block + (1|pid), data = dat)
summary(dec_mod_rc2)
### sig difference, negative decrease greater than positive decrease

## RC3
dec_mod_rc3 <- lmer(RC3 ~ block + (1|pid), data = dat)
summary(dec_mod_rc3)
### sig difference, negative decrease greater than positive decrease

## RC5
dec_mod_rc5 <- lmer(RC5 ~ block + (1|pid), data = dat)
summary(dec_mod_rc5)
### sig difference, negative decrease greater than positive decrease

## RC11
dec_mod_rc11 <- lmer(RC11 ~ block + (1|pid), data = dat)
summary(dec_mod_rc11)
### no sig difference

## RC12
dec_mod_rc12 <- lmer(RC12 ~ block + (1|pid), data = dat)
summary(dec_mod_rc12)
### sig difference, positive decrease greater than negative decrease


