# analyses

## load packages
library(tidyverse)
library(gt)
library(broom)
library(lavaan)
library(here)
library(semPlot)
library(tidySEM)
library(performance)

## read in data and pivot to wide format
dat <- read_csv(here("data", "paper_three", "dat_for_analyses_2021-04-02.csv")) %>%
  pivot_wider(names_from = block,
              values_from = RC2:RC17) %>%
  relocate(pid:race, RC2_Neg_Watch:RC17_NA) %>%
  select(-c(RC2_NA, RC3_NA, RC5_NA, RC7_NA, nRC8_NA, pRC8_NA, RC17_NA)) %>%
  filter(!is.na(RC2_Pos_Watch), # only participants with EEG data for all conditions
         pid != 22585512,       # dropped out of study
         pid != 22585452)       # only participants with EEG data for all conditions

################################
###INTERNALIZING FACTOR MODEL###
################################
internalizing_mod <- '
# internalizing measurement model
  INT =~ masq_aa + masq_pa + pswq_total
'

fit_int <- cfa(internalizing_mod,
               data = dat,
               estimator = "MLR",
               missing = "ML",
               std.lv = TRUE)

summary(fit_int, fit.measures = TRUE, standardized = TRUE)

model_performance(fit_int, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

int_params <- tidy(fit_int)

###################
## -- neg RC8 -- ##
###################
# define measurement model with neg RC8
rc8_meas_mod <- '
# neural responsivity factor
  NR =~ nRC8_Neu_Watch + nRC8_Pos_Watch + nRC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
'

mod_meas_rc8 <- cfa(rc8_meas_mod,
                    data = dat, estimator = "MLR",
                    missing = "ML",
                    orthogonal = TRUE,
                    std.lv = TRUE)

summary(mod_meas_rc8, fit.measures = TRUE, standardized = TRUE)

rc8_meas_params <- tidy(mod_meas_rc8)

rc8_meas_fit <- model_performance(mod_meas_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(mod_meas_rc8,
         what = "diagram",
         style = "lisrel",
         whatLabels = "est",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

# define eci model with neg RC8
rc8_eci_mod <- '
# neural responsivity factor
  NR =~ nRC8_Neu_Watch + nRC8_Pos_Watch + nRC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ pswq_total + masq_pa + masq_aa
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ nRC8_Pos_Watch + nRC8_Neg_Watch
'

fit_eci_rc8 <- cfa(rc8_eci_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_eci_rc8, fit.measures = TRUE, standardized = TRUE)

rc8_eci_params <- tidy(fit_eci_rc8)

rc8_eci_fit <- model_performance(fit_eci_rc8,
                                 metrics = c("RMSEA", "SRMR", "CFI",
                                             "NNFI", "AIC", "BIC"))

semPaths(fit_eci_rc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_eci_meas <- lavTestLRT(mod_meas_rc8, fit_eci_rc8)

# define neural reactivity model with neg RC8
rc8_nr_mod <- '
# neural responsivity factor
  NR =~ nRC8_Neu_Watch + nRC8_Pos_Watch + nRC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no association between nr and int factory but parameter estimate between nr and pa
  NR ~~ 0*INT + masq_pa
'

fit_nr_rc8 <- cfa(rc8_nr_mod,
                  data = dat,
                  estimator = "MLR",
                  missing = "ML",
                  std.lv = TRUE)

summary(fit_nr_rc8, fit.measures = TRUE, standardized = TRUE)

rc8_nr_params <- tidy(fit_nr_rc8)

rc8_nr_fit <- model_performance(fit_nr_rc8,
                                metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_nr_rc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "PA", "AA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_nr_meas <- lavTestLRT(mod_meas_rc8, fit_nr_rc8)

# define internalizing model with RC8
rc8_int_mod <- '
# neural responsivity factor
  NR =~ nRC8_Neu_Watch + nRC8_Pos_Watch + nRC8_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
  INT ~~ nRC8_Pos_Watch + nRC8_Neg_Watch
'

fit_int_rc8 <- cfa(rc8_int_mod,     ##### HEYWOOD CASE
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_int_rc8, fit.measures = TRUE, standardized = TRUE)

nrc8_int_params <- tidy(fit_int_rc8)

nrc8_int_fit <- model_performance(fit_int_rc8,
                                  metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_int_rc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "PA", "AA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_int_meas <- lavTestLRT(mod_meas_rc8, fit_int_rc8)

# define anxiety model with neg RC8
rc8_anx_mod <- '
# neural responsivity factor
  NR =~ nRC8_Neu_Watch + nRC8_Pos_Watch + nRC8_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
masq_aa ~~ nRC8_Pos_Watch + nRC8_Neg_Watch
'

fit_anx_rc8 <- cfa(rc8_anx_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_anx_rc8, fit.measures = TRUE, standardized = TRUE)

rc8_anx_params <- tidy(fit_anx_rc8)

rc8_anx_fit <- model_performance(fit_anx_rc8,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anx_rc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_anx_meas <- lavTestLRT(mod_meas_rc8, fit_anx_rc8)
### comparison model does NOT fit the data significantly better

# define anxiety model with neg RC8
rc8_anxapp_mod <- '
# neural responsivity factor
  NR =~ nRC8_Neu_Watch + nRC8_Pos_Watch + nRC8_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
pswq_total ~~ nRC8_Pos_Watch + nRC8_Neg_Watch
'

fit_anxapp_rc8 <- cfa(rc8_anxapp_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_anxapp_rc8, fit.measures = TRUE, standardized = TRUE)

rc8_anxapp_params <- tidy(fit_anxapp_rc8)

rc8_anxapp_fit <- model_performance(fit_anxapp_rc8,
                                    metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anxapp_rc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_anx_meas <- lavTestLRT(mod_meas_rc8, rc8_anxapp_fit)
### comparison model does NOT fit the data significantly better

###################
## -- pos RC8 -- ##
###################
# define measurement model with neg RC8
prc8_meas_mod <- '
# neural responsivity factor
  NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
'

mod_meas_prc8 <- cfa(prc8_meas_mod,
                     data = dat, estimator = "MLR",
                     missing = "ML",
                     orthogonal = TRUE,
                     std.lv = TRUE)

summary(mod_meas_prc8, fit.measures = TRUE, standardized = TRUE)

prc8_meas_params <- tidy(mod_meas_prc8)

prc8_meas_fit <- model_performance(mod_meas_prc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(mod_meas_prc8,
         what = "diagram",
         style = "lisrel",
         whatLabels = "est",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

# define eci model with pos RC8
prc8_eci_mod <- '
# neural responsivity factor
  NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ pswq_total + masq_pa + masq_aa
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ pRC8_Pos_Watch + pRC8_Neg_Watch
'

fit_eci_prc8 <- cfa(prc8_eci_mod,
                    data = dat,
                    estimator = "MLR",
                    missing = "ML",
                    std.lv = TRUE)

summary(fit_eci_prc8, fit.measures = TRUE, standardized = TRUE)

prc8_eci_params <- tidy(fit_eci_prc8)

rc8_eci_fit <- model_performance(fit_eci_prc8,
                                 metrics = c("RMSEA", "SRMR", "CFI",
                                             "NNFI", "AIC", "BIC"))

semPaths(fit_eci_prc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_prc8_eci_meas <- lavTestLRT(mod_meas_prc8, fit_eci_prc8)

# define neural reactivity model with pos RC8
prc8_nr_mod <- '
# neural responsivity factor
  NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no association between nr and int factory but parameter estimate between nr and pa
  NR ~~ 0*INT + masq_pa
'

fit_nr_prc8 <- cfa(prc8_nr_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_nr_prc8, fit.measures = TRUE, standardized = TRUE)

prc8_nr_params <- tidy(fit_nr_prc8)

prc8_nr_fit <- model_performance(fit_nr_prc8,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_nr_prc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "PA", "AA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_prc8_nr_meas <- lavTestLRT(mod_meas_prc8, fit_nr_prc8)

# define internalizing model with RC8
prc8_int_mod <- '
# neural responsivity factor
  NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
  INT ~~ pRC8_Pos_Watch + pRC8_Neg_Watch
'

fit_int_prc8 <- cfa(prc8_int_mod,
                    data = dat,
                    estimator = "MLR",
                    missing = "ML",
                    std.lv = TRUE)

summary(fit_int_prc8, fit.measures = TRUE, standardized = TRUE)

prc8_int_params <- tidy(fit_int_prc8)

prc8_int_fit <- model_performance(fit_int_prc8,
                                  metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_int_prc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "PA", "AA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_prc8_int_meas <- lavTestLRT(mod_meas_prc8, fit_int_prc8)

# define anxiety model with pos RC8
prc8_anx_mod <- '
# neural responsivity factor
  NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
masq_aa ~~ pRC8_Pos_Watch + pRC8_Neg_Watch
'

fit_anx_prc8 <- cfa(prc8_anx_mod,
                    data = dat,
                    estimator = "MLR",
                    missing = "ML",
                    std.lv = TRUE)

summary(fit_anx_prc8, fit.measures = TRUE, standardized = TRUE)

prc8_anx_params <- tidy(fit_anx_prc8)

prc8_anx_fit <- model_performance(fit_anx_prc8,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anx_prc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_anx_meas <- lavTestLRT(mod_meas_prc8, fit_anx_prc8)
### comparison model does NOT fit the data significantly better

# define anxiety model with pos RC8
prc8_anxapp_mod <- '
# neural responsivity factor
  NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
pswq_total ~~ pRC8_Pos_Watch + pRC8_Neg_Watch
'

fit_anxapp_prc8 <- cfa(prc8_anxapp_mod,
                       data = dat,
                       estimator = "MLR",
                       missing = "ML",
                       std.lv = TRUE)

summary(fit_anxapp_prc8, fit.measures = TRUE, standardized = TRUE)

prc8_anxapp_params <- tidy(fit_anxapp_prc8)

prc8_anxapp_fit <- model_performance(fit_anxapp_prc8,
                                     metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anxapp_prc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_anx_meas <- lavTestLRT(mod_meas_prc8, fit_anxapp_prc8)
### comparison model does NOT fit the data significantly better

################
## -- RC2 -- ##
################
rc2_meas_mod <- '
# neural responsivity factor
  NR =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  NR ~~ 0*INT
'

fit_meas_rc2 <- cfa(rc2_meas_mod,
                    data = dat,
                    estimator = "MLR",
                    missing = "ML",
                    std.lv = TRUE)

summary(fit_meas_rc2, fit.measures = TRUE, standardized = TRUE)

rc2_meas_params <- tidy(fit_meas_rc2)

rc2_meas_fit <- model_performance(fit_meas_rc2,
                                  metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_meas_rc2,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

# define eci model with RC2
rc2_eci_mod <- '
# neural responsivity factor
  NR =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ RC2_Pos_Watch + RC2_Neg_Watch
'

fit_eci_rc2 <- cfa(rc2_eci_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_eci_rc2, fit.measures = TRUE, standardized = TRUE)

rc2_eci_params <- tidy(fit_eci_rc2)

rc2_eci_fit <- model_performance(fit_eci_rc2,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_eci_rc2,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_rc2_eci_meas <- lavTestLRT(fit_meas_rc2, rc2_eci_mod)
### comparison model does NOT fit the data significantly better

# define neural reactivity model with RC2
rc2_nr_mod <- '
# neural responsivity factor
  NR =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT + masq_pa
'

fit_nr_rc2 <- cfa(rc2_nr_mod,
                  data = dat,
                  estimator = "MLR",
                  missing = "ML",
                  std.lv = TRUE)

summary(fit_nr_rc2, fit.measures = TRUE, standardized = TRUE)

rc2_nr_params <- tidy(fit_nr_rc2)

rc2_nr_fit <- model_performance(fit_nr_rc2,
                                metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

## model comparison test with measurement model
chisq_rc2_nr_meas <- lavTestLRT(fit_meas_rc2, fit_nr_rc2)
### comparison model does NOT fit the data significantly better

# define internalizing model with RC2
rc2_int_mod <- '
# neural responsivity factor
  NR =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
INT ~~ RC2_Pos_Watch + RC2_Neg_Watch
'

fit_int_rc2 <- cfa(rc2_int_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_int_rc2, fit.measures = TRUE, standardized = TRUE)

rc2_int_params <- tidy(fit_int_rc2)

rc2_int_fit <- model_performance(fit_int_rc2,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_int_rc2,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chisq_rc2_int_meas <- lavTestLRT(fit_meas_rc2, fit_int_rc2)
### comparison model does NOT fit the data significantly better

# define anxiety model with RC2
rc2_anx_mod <- '
# neural responsivity factor
  NR =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
masq_aa ~~ RC2_Pos_Watch + RC2_Neg_Watch
'

fit_anx_rc2 <- cfa(rc2_anx_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_anx_rc2, fit.measures = TRUE, standardized = TRUE)

rc2_anx_params <- tidy(fit_anx_rc2)

rc2_anx_fit <- model_performance(fit_anx_rc2,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anx_rc2,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_rc2_anx_meas <- lavTestLRT(fit_meas_rc2, fit_anx_rc2)
### comparison model does NOT fit the data significantly better

# define anxiety model with pos RC8
rc2_anxapp_mod <- '
# neural responsivity factor
  NR =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
pswq_total ~~ RC2_Pos_Watch + RC2_Neg_Watch
'

fit_anxapp_rc2 <- cfa(rc2_anxapp_mod,
                       data = dat,
                       estimator = "MLR",
                       missing = "ML",
                       std.lv = TRUE)

summary(fit_anxapp_rc2, fit.measures = TRUE, standardized = TRUE)

rc2_anxapp_params <- tidy(fit_anxapp_rc2)

rc2_anxapp_fit <- model_performance(fit_anxapp_rc2,
                                    metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anxapp_rc2,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_anx_meas <- lavTestLRT(mod_meas_rc2, fit_anxapp_rc2)
### comparison model does NOT fit the data significantly better

################
## -- RC3 -- ##
################
# define measurement model with RC3
rc3_meas_mod <- '
# neural responsivity factor
  NR =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  NR ~~ 0*INT
'

fit_meas_rc3 <- cfa(rc3_meas_mod,
                    data = dat,
                    estimator = "MLR",
                    missing = "ML",
                    std.lv = TRUE)

summary(fit_meas_rc3, fit.measures = TRUE, standardized = TRUE)

rc3_meas_params <- tidy(fit_meas_rc3)

rc3_meas_fit <- model_performance(fit_meas_rc3,
                                  metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_meas_rc3)

# define eci model with RC3
rc3_eci_mod <- '
# neural responsivity factor
  NR =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ RC3_Pos_Watch + RC3_Neg_Watch
'

fit_eci_rc3 <- cfa(rc3_eci_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_eci_rc3, fit.measures = TRUE, standardized = TRUE)

rc3_eci_params <- tidy(fit_eci_rc3)

rc3_eci_fit <- model_performance(fit_eci_rc3,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

## model comparison test with measurement model
chisq_rc3_eci_meas <- lavTestLRT(fit_meas_rc3, fit_eci_rc3)
### comparison model does NOT fit the data significantly better

# define neural reactivity model with RC3
rc3_nr_mod <- '
# neural responsivity factor
  NR =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no association between nr and int factory but parameter estimate between nr and pa
  NR ~~ 0*INT + masq_pa
'

fit_nr_rc3 <- cfa(rc3_nr_mod,
                  data = dat,
                  estimator = "MLR",
                  missing = "ML",
                  std.lv = TRUE)

summary(fit_nr_rc3, fit.measures = TRUE, standardized = TRUE)

rc3_nr_params <- tidy(fit_nr_rc3)

rc3_nr_fit <- model_performance(fit_nr_rc3,
                                metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

## model comparison test with measurement model
chisq_rc3_nr_meas <- lavTestLRT(fit_meas_rc3, fit_nr_rc3)
### comparison model does NOT fit the data significantly better

# define internalizing model with RC3
rc3_int_mod <- '
# neural responsivity factor
  NR =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
  INT ~~ RC3_Pos_Watch + RC3_Neg_Watch
'

fit_int_rc3 <- cfa(rc3_int_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_int_rc3, fit.measures = TRUE, standardized = TRUE)

rc3_int_params <- tidy(fit_int_rc3)

rc3_int_fit <- model_performance(fit_int_rc3,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

## model comparison test with measurement model
chisq_int_meas <- lavTestLRT(fit_meas_rc3, fit_int_rc3)
### comparison model does NOT fit the data significantly better

# define anxiety model with RC3
rc3_anx_mod <- '
# neural responsivity factor
  NR =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
  masq_aa ~~ RC3_Pos_Watch + RC3_Neg_Watch
'

fit_anx_rc3 <- cfa(rc3_anx_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_anx_rc3, fit.measures = TRUE, standardized = TRUE)

rc3_anx_params <- tidy(fit_anx_rc2)

rc3_anx_fit <- model_performance(fit_anx_rc3,
                                 metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anx_rc3,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_rc3_anx_meas <- lavTestLRT(fit_meas_rc3, fit_anx_rc3)
### comparison model does NOT fit the data significantly better

# define anxiety model with pos RC8
rc3_anxapp_mod <- '
# neural responsivity factor
  NR =~ rc3_Neu_Watch + rc3_Pos_Watch + rc3_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
pswq_total ~~ rc3_Pos_Watch + rc3_Neg_Watch
'

fit_anxapp_prc8 <- cfa(prc8_anxapp_mod,
                       data = dat,
                       estimator = "MLR",
                       missing = "ML",
                       std.lv = TRUE)

summary(fit_anxapp_prc8, fit.measures = TRUE, standardized = TRUE)

prc8_anxapp_params <- tidy(fit_anxapp_prc8)

prc8_anxapp_fit <- model_performance(fit_anxapp_prc8,
                                     metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anxapp_prc8,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_anx_meas <- lavTestLRT(mod_meas_prc8, fit_anxapp_prc8)
### comparison model does NOT fit the data significantly better
# define anxiety model with rc3
rc3_anxapp_mod <- '
# neural responsivity factor
  NR =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
pswq_total ~~ RC3_Pos_Watch + RC3_Neg_Watch
'

fit_anxapp_rc3 <- cfa(rc3_anxapp_mod,
                       data = dat,
                       estimator = "MLR",
                       missing = "ML",
                       std.lv = TRUE)

summary(fit_anxapp_rc3, fit.measures = TRUE, standardized = TRUE)

rc3_anxapp_params <- tidy(fit_anxapp_rc3)

rc3_anxapp_fit <- model_performance(fit_anxapp_rc3,
                                     metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))

semPaths(fit_anxapp_rc3,
         what = "diagram",
         whatLabels = "est",
         style = "lisrel",
         nCharNodes = 0,
         intercepts = FALSE,
         sizeMan = 8,
         sizeLat = 10,
         nodeLabels = c("NEU", "POS", "NEG",
                        "AA", "PA", "WOR",
                        "NR", "INT"))

## model comparison test with measurement model
chi_nrc8_anx_meas <- lavTestLRT(mod_meas_rc3, fit_anxapp_rc3)
### comparison model does NOT fit the data significantly better

####################################
#############TABLES#################
####################################
library(kableExtra)
# p value formating function
num_format <- function(val) {sub("^(-?)0.", "\\1.", sprintf("%.3f", val))}

# parameter estimates for internalizing factor
int_param_table <-
  int_params %>%
    select(term, estimate, std.all, std.error, p.value) %>%
    filter(str_detect(term, "INT =~")) %>%
    mutate(across(.cols = c(estimate, std.all, std.error),
                  .fns = ~ sprintf("%.2f", .x))) %>%
    mutate(p.value = if_else(p.value < .001, "<.001",
                             num_format(p.value))) %>%
    unite("Est/Std", estimate, std.all, sep = "/") %>%
    rename("Path" = "term",
           "$SE$" = "std.error",
           "$p$" = "p.value") %>%
  mutate(Path = str_replace(Path, "=~", "→"),
         Path = str_replace(Path, "masq_aa", "Anx Aro"),
         Path = str_replace(Path, "masq_pa", "PA"),
         Path = str_replace(Path, "pswq_total", "Anx App"))

int_fac_table_wide <- bind_cols(int_param_table, # repeat these three times so they can be merged
                                int_param_table, # with the three ERP components
                                int_param_table)

# parameter estimates for erp measurement components
meas_params_list <- list(rc8_meas_params,
                         rc2_meas_params,
                         rc3_meas_params)

meas_params_list <-
  map(meas_params_list, ~{
    .x %>%
      select(term, estimate, std.all, std.error, p.value) %>%
      filter(str_detect(term, "NR =~")) %>%
      mutate(across(.cols = c(estimate, std.all, std.error),
                    .fns = ~ sprintf("%.2f", .x))) %>%
      mutate(p.value = if_else(p.value < .001, "<.001",
                               num_format(p.value))) %>%
      unite("Est/Std", estimate, std.all, sep = "/") %>%
      rename("Path" = "term",
             "$SE$" = "std.error",
             "$p$" = "p.value") %>%
      mutate(Path = c("NR $\\rightarrow$ Neutral",
                      "NR $\\rightarrow$ Positive",
                      "NR $\\rightarrow$ Negative"))
})

meas_params_wide <- bind_cols(meas_params_list)

# parameter estimates for erp eci components
eci_params_list <- list(rc8_eci_params,
                        rc2_eci_params,
                        rc3_eci_params)

eci_params_list <-
  map(eci_params_list, ~{
    .x %>%
      select(term, estimate, std.all, std.error, p.value) %>%
      filter(str_detect(term, "~~ masq_pa"),
             term != "masq_pa ~~ masq_pa") %>%
      mutate(across(.cols = c(estimate, std.all, std.error),
                    .fns = ~ sprintf("%.2f", .x))) %>%
      mutate(p.value = if_else(p.value < .001, "<.001",
                               num_format(p.value))) %>%
      unite("Est/Std", estimate, std.all, sep = "/") %>%
      rename("Path" = "term",
             "$SE$" = "std.error",
             "$p$" = "p.value") %>%
      mutate(Path = c("Positive $\\leftrightarrow$ PA",
                      "Negative $\\leftrightarrow$ PA"))
})

eci_params_wide <- bind_cols(eci_params_list)

# parameter estimates for NR components
nr_params_list <- list(rc8_nr_params,
                       rc2_nr_params,
                       rc3_nr_params)

nr_params_list <-
  map(nr_params_list, ~{
    .x %>%
      select(term, estimate, std.all, std.error, p.value) %>%
      filter(term == "NR ~~ masq_pa") %>%
      mutate(across(.cols = c(estimate, std.all, std.error),
                    .fns = ~ sprintf("%.2f", .x))) %>%
      mutate(p.value = if_else(p.value < .001, "<.001",
                               num_format(p.value))) %>%
      unite("Est/Std", estimate, std.all, sep = "/") %>%
      rename("Path" = "term",
             "$SE$" = "std.error",
             "$p$" = "p.value") %>%
      mutate(Path = "NR $\\leftrightarrow$ PA")
  })

nr_params_wide <- bind_cols(nr_params_list)

# parameter estimates for internalizing components
int_params_list <- list(nrc8_int_params,
                        rc2_int_params,
                        rc3_int_params)

int_params_list <-
  map(int_params_list, ~{
    .x %>%
      select(term, estimate, std.all, std.error, p.value) %>%
      filter(str_detect(term, "INT ~~"),
             term != "INT ~~ INT") %>%
      mutate(across(.cols = c(estimate, std.all, std.error),
                    .fns = ~ sprintf("%.2f", .x))) %>%
      mutate(p.value = if_else(p.value < .001, "<.001",
                               num_format(p.value))) %>%
      unite("Est/Std", estimate, std.all, sep = "/") %>%
      rename("Path" = "term",
             "$SE$" = "std.error",
             "$p$" = "p.value") %>%
      mutate(Path = c("INT $\\leftrightarrow$ Positive",
                      "INT $\\leftrightarrow$ Negative"))
  })

int_params_wide <- bind_cols(int_params_list)

# parameter estimates for anxiety components
anx_params_list <- list(rc8_anx_params,
                        rc2_anx_params,
                        rc3_anx_params)

anx_params_list <-
  map(anx_params_list, ~{
    .x %>%
      select(term, estimate, std.all, std.error, p.value) %>%
      filter(str_detect(term, "Watch ~~ masq_aa")) %>%
      mutate(across(.cols = c(estimate, std.all, std.error),
                    .fns = ~ sprintf("%.2f", .x))) %>%
      mutate(p.value = if_else(p.value < .001, "<.001",
                               num_format(p.value))) %>%
      unite("Est/Std", estimate, std.all, sep = "/") %>%
      rename("Path" = "term",
             "$SE$" = "std.error",
             "$p$" = "p.value") %>%
      mutate(Path = c("Positive $\\leftrightarrow$ Anx Aro",
                      "Negative $\\leftrightarrow$ Anx Aro"))
  })

anx_params_wide <- bind_cols(anx_params_list)

# parameter estimates for anxious apprehension components
anxapp_params_list <- list(rc8_anxapp_params,
                           rc2_anxapp_params,
                           rc3_anxapp_params)

anxapp_params_list <-
  map(anxapp_params_list, ~{
    .x %>%
      select(term, estimate, std.all, std.error, p.value) %>%
      filter(str_detect(term, "Watch ~~ pswq_total")) %>%
      mutate(across(.cols = c(estimate, std.all, std.error),
                    .fns = ~ sprintf("%.2f", .x))) %>%
      mutate(p.value = if_else(p.value < .001, "<.001",
                               num_format(p.value))) %>%
      unite("Est/Std", estimate, std.all, sep = "/") %>%
      rename("Path" = "term",
             "$SE$" = "std.error",
             "$p$" = "p.value") %>%
      mutate(Path = c("Positive $\\leftrightarrow$ Anx App",
                      "Negative $\\leftrightarrow$ Anx App"))
  })

anxapp_params_wide <- bind_cols(anxapp_params_list)

# bind all wide tables together
tab <- bind_rows(int_fac_table_wide,
                 meas_params_wide,
                 eci_params_wide,
                 nr_params_wide,
                 int_params_wide,
                 anx_params_wide,
                 anxapp_params_wide) %>%
       select(-c("Path...5",
                 "Path...9"))

names(tab) <- gsub("...[[:digit:]]+", "", names(tab))

# construct the table
kable(tab, "latex", escape = FALSE, booktabs = TRUE, align = c("l", rep("r", 9)), linesep = "") %>%
  kable_styling(font_size = 12, latex_options = c("scale_down")) %>%
  add_header_above(c(" ", "257 ms Comp" = 3, "371 ms Comp" = 3, "736 ms Comp" = 3),
                   bold = TRUE,
                   italic = TRUE) %>%
  pack_rows("Measurement Model", 1, 6) %>%
  pack_rows("ECI Model", 7, 8) %>%
  pack_rows("NR Model", 9, 9) %>%
  pack_rows("INT Model", 10, 11) %>%
  pack_rows("Anxious Arousal Model", 12, 13) %>%
  pack_rows("Anxious Apprehension Model", 14, 15) %>%
  row_spec(0, align = "c") %>%
  landscape() %>%
  footnote(general = "Path labels correspond to a parameter estimates in each model. Group headings in the 'path' denote the specific model that the parameter estimates belong to. Note that the factor loadings in the measurement model are present in each of the subsequent competing models. $\\SE$ = Standard Error, Est/Std = undstandardized and standardized parameter estimate, $\\\\rightarrow$ = latent factor loading estimate,
           \n$\\\\leftrightarrow$ = covariance estimate.",
           threeparttable = TRUE,
           escape = FALSE,
           general_title = "Note.",
           footnote_as_chunk = TRUE) %>%
  save_kable("sem_param_table.pdf")
