# Pairwise comparisons for positive/negative images within increase/decrease conditions

# load packages
library(tidyverse)
library(emmeans)
library(lmerTest)
library(r2glmm)
library(effectsize)
library(kableExtra)

dat <- read_csv("data/paper_two/created_data/temp_fac_score_dat_analyses2021-02-11.csv")

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
effectsize(inc_mod_rc3)
interpret_d(0.23, "cohen1988")
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
effectsize(inc_mod_rc12)
interpret_d(0.36, rules = "cohen1988")
### sig difference, positive increase greater than negative increase

## positive RC12
inc_mod_pos_rc12 <- lmer(pos_RC12 ~ block + (1|pid), data = dat)
summary(inc_mod_pos_rc12)
### sig difference, positive increase greater than negative increase

# run analyses of decrease conditions
dat$block <- relevel(dat$block, "Neg_Dec")

## RC2
dec_mod_rc2 <- lmer(RC2 ~ block + (1|pid), data = dat)
summary(dec_mod_rc2)
effectsize(dec_mod_rc2)
interpret_d(-0.23, "cohen1988")
### sig difference, negative decrease greater than positive decrease

## RC3
dec_mod_rc3 <- lmer(RC3 ~ block + (1|pid), data = dat)
summary(dec_mod_rc3)
effectsize(dec_mod_rc3)
interpret_d(-0.51, "cohen1988")
### sig difference, negative decrease greater than positive decrease

## RC5
dec_mod_rc5 <- lmer(RC5 ~ block + (1|pid), data = dat)
summary(dec_mod_rc5)
### no sig difference

## RC11
dec_mod_rc11 <- lmer(RC11 ~ block + (1|pid), data = dat)
summary(dec_mod_rc11)
### no sig difference

## RC12
dec_mod_rc12 <- lmer(RC12 ~ block + (1|pid), data = dat)
summary(dec_mod_rc12)
### sig difference, positive decrease greater than negative decrease

## positive RC12
dec_mod_pos_rc12 <- lmer(pos_RC12 ~ block + (1|pid), data = dat)
summary(dec_mod_pos_rc12)
### no sig difference
