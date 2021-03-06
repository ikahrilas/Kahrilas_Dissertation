## Path analyses
## Author: Ian J. Kahrilas
## This script runs the path analysis in the manuscript. The scripts "Loading and Cleaning Data.R" and "Deriving subscales.R"
## should be run first for the following code to work properly.
  
# load packages
library(readr)
library(tidyverse)
library(lavaan) # install.packages("lavaan")

# load in data
sad_data <- read_csv("total_data.csv")

# construct model
# We have two predictors, positive affect (MASQAD14) and negative affect (MASQAD8) and an outcome of depression severity (PHQDEP). 
# This relationship is mediated by the three temporal domains of savoring: reminiscing (SBI_REM), savoring the moment (SAV_MOM), 
# and anticipating (SBI_ANT). We also include the covariates of sex, worry, and anxious arousal.

sad_model <- '# covariates
             reminiscing ~ sex + pswq_total + masq_aa
             savoring_moment ~ sex + pswq_total + masq_aa
             anticipating ~ sex + pswq_total + masq_aa
             phq_total ~ sex + pswq_total + masq_aa
             masq_pa ~ sex + pswq_total + masq_aa
             masq_na ~ sex + pswq_total + masq_aa
           # direct effects
             phq_total ~ c1*masq_pa
             phq_total ~ ci*masq_na
           # mediators
             reminiscing ~ a1*masq_pa
             savoring_moment ~ a2*masq_pa
             anticipating ~ a3*masq_pa
             phq_total ~ b1*reminiscing
             phq_total ~ b2*savoring_moment
             phq_total ~ b3*anticipating
             reminiscing ~ ai*masq_na
             savoring_moment ~ aii*masq_na
             anticipating ~ aiii*masq_na
           # indirect effects (a*b)
             pa_rem := a1*b1
             pa_mom := a2*b2
             pa_ant := a3*b3
             na_rem := ai*b1
             na_mom := aii*b2
             na_ant := aiii*b3
           # total effects
             pa_total:= c1 + (a1*b1) + (a2*b2) + (a3*b3)
             na_total := ci + (ai*b1) + (aii*b2) + (aii*b3)
           # exogenous variable correlations
             masq_pa ~~ masq_na
           # residual correlations
             anticipating ~~ savoring_moment + reminiscing
             savoring_moment ~~ reminiscing
'

set.seed(123)
fit <- sem(sad_model, data = sad_data, se = "bootstrap", bootstrap = 10000, meanstructure = TRUE)
summary(fit, fit.measures = TRUE, ci = TRUE, standardized = TRUE, rsquare = TRUE)
parms <- parameterEstimates(fit, boot.ci.type = "bca.simple", standardized = TRUE)
sad_rsquare <- lavInspect(fit, "rsquare")
# save parameter dataframe to working directory so that it can be accessed in the table compiler script.
write_csv(parms, "parms.csv")

# multivariate regressions
# the following package adds standardized regression coefficients to lm output
library(lm.beta) # install.packages("lm.beta")

# multiple regression to establish individual relations between PA, NA, and AA with depression
multi_mod <- lm(phq_total ~ masq_pa + masq_na + masq_aa, data = sad_data)
summary(multi_mod)
lm.beta(multi_mod)
