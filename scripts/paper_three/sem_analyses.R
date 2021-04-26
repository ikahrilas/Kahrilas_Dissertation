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

rc8_meas_fit <- fitMeasures(mod_meas_rc8, c("chisq.scaled",
                                            "df",
                                            "rmsea.scaled",
                                            "srmr",
                                            "cfi.scaled",
                                            "nnfi.scaled",
                                            "aic",
                                            "bic"))

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

rc8_eci_fit <-   fitMeasures(fit_eci_rc8, c("chisq.scaled",
                                             "df",
                                             "rmsea.scaled",
                                             "srmr",
                                             "cfi.scaled",
                                             "nnfi.scaled",
                                             "aic",
                                             "bic"))

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

rc8_nr_fit <- fitMeasures(fit_nr_rc8, c("chisq.scaled",
                                          "df",
                                          "rmsea.scaled",
                                          "srmr",
                                          "cfi.scaled",
                                          "nnfi.scaled",
                                          "aic",
                                          "bic"))

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

nrc8_int_fit <- fitMeasures(fit_int_rc8, c("chisq.scaled",
                                            "df",
                                            "rmsea.scaled",
                                            "srmr",
                                            "cfi.scaled",
                                            "nnfi.scaled",
                                            "aic",
                                            "bic"))

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

rc8_anx_fit <- fitMeasures(fit_anx_rc8, c("chisq.scaled",
                                           "df",
                                           "rmsea.scaled",
                                           "srmr",
                                           "cfi.scaled",
                                           "nnfi.scaled",
                                           "aic",
                                           "bic"))

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

rc8_anxapp_fit <- fitMeasures(fit_anxapp_rc8, c("chisq.scaled",
                                             "df",
                                             "rmsea.scaled",
                                             "srmr",
                                             "cfi.scaled",
                                             "nnfi.scaled",
                                             "aic",
                                             "bic"))

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
chi_nrc8_anxapp_meas <- lavTestLRT(mod_meas_rc8, rc8_anxapp_fit)
### comparison model does NOT fit the data significantly better

# RC8 brightening effect model
rc8_be_mod <- '
# neural responsivity factor
  NR =~ nRC8_Neu_Watch + nRC8_Pos_Watch + nRC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ pswq_total + masq_pa + masq_aa
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ nRC8_Pos_Watch
'

rc8_be_mod <- cfa(rc8_be_mod,
                  data = dat,
                  estimator = "MLR",
                  missing = "ML",
                  std.lv = TRUE)

summary(rc8_be_mod, fit.measures = TRUE, standardized = TRUE)

rc8_be_params <- tidy(rc8_be_mod)

rc8_be_fit <-   fitMeasures(rc8_be_mod, c("chisq.scaled",
                                            "df",
                                            "rmsea.scaled",
                                            "srmr",
                                            "cfi.scaled",
                                            "nnfi.scaled",
                                            "aic",
                                            "bic"))
modificationIndices(mod_meas_rc8, sort. = TRUE)
semPaths(rc8_be_mod,
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
chi_nrc8_be_meas <- lavTestLRT(mod_meas_rc8, rc8_be_mod)
chi_nrc8_be_eci <- lavTestLRT(fit_eci_rc8, rc8_be_mod)

# ###################
# ## -- pos RC8 -- ##
# ###################
# # define measurement model with neg RC8
# prc8_meas_mod <- '
# # neural responsivity factor
#   NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# # internalizing measurement model with MASQ subscales
#   INT =~ masq_aa + masq_pa + pswq_total
# '
#
# mod_meas_prc8 <- cfa(prc8_meas_mod,
#                      data = dat, estimator = "MLR",
#                      missing = "ML",
#                      orthogonal = TRUE,
#                      std.lv = TRUE)
#
# summary(mod_meas_prc8, fit.measures = TRUE, standardized = TRUE)
#
# prc8_meas_params <- tidy(mod_meas_prc8)
#
# prc8_meas_fit <- model_performance(mod_meas_prc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))
#
# semPaths(mod_meas_prc8,
#          what = "diagram",
#          style = "lisrel",
#          whatLabels = "est",
#          nCharNodes = 0,
#          intercepts = FALSE,
#          sizeMan = 8,
#          sizeLat = 10,
#          nodeLabels = c("NEU", "POS", "NEG",
#                         "AA", "PA", "WOR",
#                         "NR", "INT"))
#
# # define eci model with pos RC8
# prc8_eci_mod <- '
# # neural responsivity factor
#   NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# # internalizing measurement model with MASQ subscales
#   INT =~ pswq_total + masq_pa + masq_aa
# # orthogonal model
#   NR ~~ 0*INT
# # residual covariances (constrained to equality)
#   masq_pa ~~ pRC8_Pos_Watch + pRC8_Neg_Watch
# '
#
# fit_eci_prc8 <- cfa(prc8_eci_mod,
#                     data = dat,
#                     estimator = "MLR",
#                     missing = "ML",
#                     std.lv = TRUE)
#
# summary(fit_eci_prc8, fit.measures = TRUE, standardized = TRUE)
#
# prc8_eci_params <- tidy(fit_eci_prc8)
#
# rc8_eci_fit <- model_performance(fit_eci_prc8,
#                                  metrics = c("RMSEA", "SRMR", "CFI",
#                                              "NNFI", "AIC", "BIC"))
#
# semPaths(fit_eci_prc8,
#          what = "diagram",
#          whatLabels = "est",
#          style = "lisrel",
#          nCharNodes = 0,
#          intercepts = FALSE,
#          sizeMan = 8,
#          sizeLat = 10,
#          nodeLabels = c("NEU", "POS", "NEG",
#                         "AA", "PA", "WOR",
#                         "NR", "INT"))
#
# ## model comparison test with measurement model
# chi_prc8_eci_meas <- lavTestLRT(mod_meas_prc8, fit_eci_prc8)
#
# # define neural reactivity model with pos RC8
# prc8_nr_mod <- '
# # neural responsivity factor
#   NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# # internalizing measurement model with MASQ subscales
#   INT =~ masq_aa + masq_pa + pswq_total
# # no association between nr and int factory but parameter estimate between nr and pa
#   NR ~~ 0*INT + masq_pa
# '
#
# fit_nr_prc8 <- cfa(prc8_nr_mod,
#                    data = dat,
#                    estimator = "MLR",
#                    missing = "ML",
#                    std.lv = TRUE)
#
# summary(fit_nr_prc8, fit.measures = TRUE, standardized = TRUE)
#
# prc8_nr_params <- tidy(fit_nr_prc8)
#
# prc8_nr_fit <- model_performance(fit_nr_prc8,
#                                  metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))
#
# semPaths(fit_nr_prc8,
#          what = "diagram",
#          whatLabels = "est",
#          style = "lisrel",
#          nCharNodes = 0,
#          intercepts = FALSE,
#          sizeMan = 8,
#          sizeLat = 10,
#          nodeLabels = c("NEU", "POS", "NEG",
#                         "PA", "AA", "WOR",
#                         "NR", "INT"))
#
# ## model comparison test with measurement model
# chi_prc8_nr_meas <- lavTestLRT(mod_meas_prc8, fit_nr_prc8)
#
# # define internalizing model with RC8
# prc8_int_mod <- '
# # neural responsivity factor
#   NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# # internalNRizing measurement model with MASQ subscales
#   INT =~ masq_aa + masq_pa + pswq_total
# # no correlation with internalizing factor but correlated with masq_pa
#   NR ~~ 0*INT
# # estimate parameters among pos and neg watch and internalizing factor constrained to equality
#   INT ~~ pRC8_Pos_Watch + pRC8_Neg_Watch
# '
#
# fit_int_prc8 <- cfa(prc8_int_mod,
#                     data = dat,
#                     estimator = "MLR",
#                     missing = "ML",
#                     std.lv = TRUE)
#
# summary(fit_int_prc8, fit.measures = TRUE, standardized = TRUE)
#
# prc8_int_params <- tidy(fit_int_prc8)
#
# prc8_int_fit <- model_performance(fit_int_prc8,
#                                   metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))
#
# semPaths(fit_int_prc8,
#          what = "diagram",
#          whatLabels = "est",
#          style = "lisrel",
#          nCharNodes = 0,
#          intercepts = FALSE,
#          sizeMan = 8,
#          sizeLat = 10,
#          nodeLabels = c("NEU", "POS", "NEG",
#                         "PA", "AA", "WOR",
#                         "NR", "INT"))
#
# ## model comparison test with measurement model
# chi_prc8_int_meas <- lavTestLRT(mod_meas_prc8, fit_int_prc8)
#
# # define anxiety model with pos RC8
# prc8_anx_mod <- '
# # neural responsivity factor
#   NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# # internalNRizing measurement model with MASQ subscales
#   INT =~ masq_aa + masq_pa + pswq_total
# # no correlation with internalizing factor but correlated with masq_pa
#   NR ~~ 0*INT
# # estimate parameters among pos and neg watch and internalizing factor constrained to equality
# masq_aa ~~ pRC8_Pos_Watch + pRC8_Neg_Watch
# '
#
# fit_anx_prc8 <- cfa(prc8_anx_mod,
#                     data = dat,
#                     estimator = "MLR",
#                     missing = "ML",
#                     std.lv = TRUE)
#
# summary(fit_anx_prc8, fit.measures = TRUE, standardized = TRUE)
#
# prc8_anx_params <- tidy(fit_anx_prc8)
#
# prc8_anx_fit <- model_performance(fit_anx_prc8,
#                                  metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))
#
# semPaths(fit_anx_prc8,
#          what = "diagram",
#          whatLabels = "est",
#          style = "lisrel",
#          nCharNodes = 0,
#          intercepts = FALSE,
#          sizeMan = 8,
#          sizeLat = 10,
#          nodeLabels = c("NEU", "POS", "NEG",
#                         "AA", "PA", "WOR",
#                         "NR", "INT"))
#
# ## model comparison test with measurement model
# chi_nrc8_anx_meas <- lavTestLRT(mod_meas_prc8, fit_anx_prc8)
# ### comparison model does NOT fit the data significantly better
#
# # define anxiety model with pos RC8
# prc8_anxapp_mod <- '
# # neural responsivity factor
#   NR =~ pRC8_Neu_Watch + pRC8_Pos_Watch + pRC8_Neg_Watch
# # internalNRizing measurement model with MASQ subscales
#   INT =~ masq_aa + masq_pa + pswq_total
# # no correlation with internalizing factor but correlated with masq_pa
#   NR ~~ 0*INT
# # estimate parameters among pos and neg watch and internalizing factor constrained to equality
# pswq_total ~~ pRC8_Pos_Watch + pRC8_Neg_Watch
# '
#
# fit_anxapp_prc8 <- cfa(prc8_anxapp_mod,
#                        data = dat,
#                        estimator = "MLR",
#                        missing = "ML",
#                        std.lv = TRUE)
#
# summary(fit_anxapp_prc8, fit.measures = TRUE, standardized = TRUE)
#
# prc8_anxapp_params <- tidy(fit_anxapp_prc8)
#
# prc8_anxapp_fit <- model_performance(fit_anxapp_prc8,
#                                      metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC"))
#
# semPaths(fit_anxapp_prc8,
#          what = "diagram",
#          whatLabels = "est",
#          style = "lisrel",
#          nCharNodes = 0,
#          intercepts = FALSE,
#          sizeMan = 8,
#          sizeLat = 10,
#          nodeLabels = c("NEU", "POS", "NEG",
#                         "AA", "PA", "WOR",
#                         "NR", "INT"))
#
# ## model comparison test with measurement model
# chi_nrc8_anx_meas <- lavTestLRT(mod_meas_prc8, fit_anxapp_prc8)
# ### comparison model does NOT fit the data significantly better

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

rc2_meas_fit <- fitMeasures(fit_meas_rc2, c("chisq.scaled",
                                              "df",
                                              "rmsea.scaled",
                                              "srmr",
                                              "cfi.scaled",
                                              "nnfi.scaled",
                                              "aic",
                                              "bic"))

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

rc2_eci_fit <- fitMeasures(fit_eci_rc2, c("chisq.scaled",
                                           "df",
                                           "rmsea.scaled",
                                           "srmr",
                                           "cfi.scaled",
                                           "nnfi.scaled",
                                           "aic",
                                           "bic"))

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
chi_rc2_eci_meas <- lavTestLRT(fit_meas_rc2, fit_eci_rc2)
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

rc2_nr_fit <- fitMeasures(fit_nr_rc2, c("chisq.scaled",
                                         "df",
                                         "rmsea.scaled",
                                         "srmr",
                                         "cfi.scaled",
                                         "nnfi.scaled",
                                         "aic",
                                         "bic"))

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

rc2_int_fit <- fitMeasures(fit_int_rc2, c("chisq.scaled",
                                         "df",
                                         "rmsea.scaled",
                                         "srmr",
                                         "cfi.scaled",
                                         "nnfi.scaled",
                                         "aic",
                                         "bic"))

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

rc2_anx_fit <- fitMeasures(fit_anx_rc2, c("chisq.scaled",
                                          "df",
                                          "rmsea.scaled",
                                          "srmr",
                                          "cfi.scaled",
                                          "nnfi.scaled",
                                          "aic",
                                          "bic"))

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

# define anxiety model with rc2
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

rc2_anxapp_fit <- fitMeasures(fit_anxapp_rc2, c("chisq.scaled",
                                             "df",
                                             "rmsea.scaled",
                                             "srmr",
                                             "cfi.scaled",
                                             "nnfi.scaled",
                                             "aic",
                                             "bic"))

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
chi_nrc8_anxapp_meas <- lavTestLRT(fit_meas_rc2, fit_anxapp_rc2)
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

rc3_meas_fit <- fitMeasures(fit_meas_rc3, c("chisq.scaled",
                                              "df",
                                              "rmsea.scaled",
                                              "srmr",
                                              "cfi.scaled",
                                              "nnfi.scaled",
                                              "aic",
                                              "bic"))

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

rc3_eci_fit <- fitMeasures(fit_eci_rc3, c("chisq.scaled",
                                           "df",
                                           "rmsea.scaled",
                                           "srmr",
                                           "cfi.scaled",
                                           "nnfi.scaled",
                                           "aic",
                                           "bic"))

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

rc3_nr_fit <- fitMeasures(fit_nr_rc3, c("chisq.scaled",
                                         "df",
                                         "rmsea.scaled",
                                         "srmr",
                                         "cfi.scaled",
                                         "nnfi.scaled",
                                         "aic",
                                         "bic"))

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

rc3_int_fit <- fitMeasures(fit_int_rc3, c("chisq.scaled",
                                         "df",
                                         "rmsea.scaled",
                                         "srmr",
                                         "cfi.scaled",
                                         "nnfi.scaled",
                                         "aic",
                                         "bic"))

## model comparison test with measurement model
chisq_rc3_int_meas <- lavTestLRT(fit_meas_rc3, fit_int_rc3)
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

rc3_anx_fit <- fitMeasures(fit_anx_rc3, c("chisq.scaled",
                                          "df",
                                          "rmsea.scaled",
                                          "srmr",
                                          "cfi.scaled",
                                          "nnfi.scaled",
                                          "aic",
                                          "bic"))

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
chisq_rc3_anx_meas <- lavTestLRT(fit_meas_rc3, fit_anx_rc3)
### comparison model does NOT fit the data significantly better

# define anxious apprehension model with rc3
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

rc3_anxapp_fit <- fitMeasures(fit_anxapp_rc3, c("chisq.scaled",
                                             "df",
                                             "rmsea.scaled",
                                             "srmr",
                                             "cfi.scaled",
                                             "nnfi.scaled",
                                             "aic",
                                             "bic"))

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
chi_rc3_anxapp_meas <- lavTestLRT(fit_meas_rc3, fit_anxapp_rc3)
### comparison model does NOT fit the data significantly better

# save all analyses to work space
save.image("data/paper_three/sem_analyses_results.RData")





################
## -- RC5 -- ##
################
RC5_meas_mod <- '
# neural responsivity factor
  NR =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  NR ~~ 0*INT
'

fit_meas_RC5 <- cfa(RC5_meas_mod,
                    data = dat,
                    estimator = "MLR",
                    missing = "ML",
                    std.lv = TRUE)

summary(fit_meas_RC5, fit.measures = TRUE, standardized = TRUE)

RC5_meas_params <- tidy(fit_meas_RC5)

RC5_meas_fit <- fitMeasures(fit_meas_RC5, c("chisq.scaled",
                                            "df",
                                            "rmsea.scaled",
                                            "srmr",
                                            "cfi.scaled",
                                            "nnfi.scaled",
                                            "aic",
                                            "bic"))

semPaths(fit_meas_RC5,
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

# define eci model with RC5
RC5_eci_mod <- '
# neural responsivity factor
  NR =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ RC5_Pos_Watch + RC5_Neg_Watch
'

fit_eci_RC5 <- cfa(RC5_eci_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_eci_RC5, fit.measures = TRUE, standardized = TRUE)

RC5_eci_params <- tidy(fit_eci_RC5)

RC5_eci_fit <- fitMeasures(fit_eci_RC5, c("chisq.scaled",
                                          "df",
                                          "rmsea.scaled",
                                          "srmr",
                                          "cfi.scaled",
                                          "nnfi.scaled",
                                          "aic",
                                          "bic"))

semPaths(fit_eci_RC5,
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
chi_RC5_eci_meas <- lavTestLRT(fit_meas_RC5, fit_eci_RC5)
### comparison model does NOT fit the data significantly better

# define neural reactivity model with RC5
RC5_nr_mod <- '
# neural responsivity factor
  NR =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT + masq_pa
'

fit_nr_RC5 <- cfa(RC5_nr_mod,
                  data = dat,
                  estimator = "MLR",
                  missing = "ML",
                  std.lv = TRUE)

summary(fit_nr_RC5, fit.measures = TRUE, standardized = TRUE)

RC5_nr_params <- tidy(fit_nr_RC5)

RC5_nr_fit <- fitMeasures(fit_nr_RC5, c("chisq.scaled",
                                        "df",
                                        "rmsea.scaled",
                                        "srmr",
                                        "cfi.scaled",
                                        "nnfi.scaled",
                                        "aic",
                                        "bic"))

## model comparison test with measurement model
chisq_RC5_nr_meas <- lavTestLRT(fit_meas_RC5, fit_nr_RC5)
### comparison model does NOT fit the data significantly better

# define internalizing model with RC5
RC5_int_mod <- '
# neural responsivity factor
  NR =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
INT ~~ RC5_Pos_Watch + RC5_Neg_Watch
'

fit_int_RC5 <- cfa(RC5_int_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_int_RC5, fit.measures = TRUE, standardized = TRUE)

RC5_int_params <- tidy(fit_int_RC5)

RC5_int_fit <- fitMeasures(fit_int_RC5, c("chisq.scaled",
                                          "df",
                                          "rmsea.scaled",
                                          "srmr",
                                          "cfi.scaled",
                                          "nnfi.scaled",
                                          "aic",
                                          "bic"))

semPaths(fit_int_RC5,
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
chisq_RC5_int_meas <- lavTestLRT(fit_meas_RC5, fit_int_RC5)
### comparison model does NOT fit the data significantly better

# define anxiety model with RC5
RC5_anx_mod <- '
# neural responsivity factor
  NR =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
masq_aa ~~ RC5_Pos_Watch + RC5_Neg_Watch
'

fit_anx_RC5 <- cfa(RC5_anx_mod,
                   data = dat,
                   estimator = "MLR",
                   missing = "ML",
                   std.lv = TRUE)

summary(fit_anx_RC5, fit.measures = TRUE, standardized = TRUE)

RC5_anx_params <- tidy(fit_anx_RC5)

RC5_anx_fit <- fitMeasures(fit_anx_RC5, c("chisq.scaled",
                                          "df",
                                          "rmsea.scaled",
                                          "srmr",
                                          "cfi.scaled",
                                          "nnfi.scaled",
                                          "aic",
                                          "bic"))

semPaths(fit_anx_RC5,
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
chi_RC5_anx_meas <- lavTestLRT(fit_meas_RC5, fit_anx_RC5)
### comparison model does NOT fit the data significantly better

# define anxiety model with RC5
RC5_anxapp_mod <- '
# neural responsivity factor
  NR =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
pswq_total ~~ RC5_Pos_Watch + RC5_Neg_Watch
'

fit_anxapp_RC5 <- cfa(RC5_anxapp_mod,
                      data = dat,
                      estimator = "MLR",
                      missing = "ML",
                      std.lv = TRUE)

summary(fit_anxapp_RC5, fit.measures = TRUE, standardized = TRUE)

RC5_anxapp_params <- tidy(fit_anxapp_RC5)

RC5_anxapp_fit <- fitMeasures(fit_anxapp_RC5, c("chisq.scaled",
                                                "df",
                                                "rmsea.scaled",
                                                "srmr",
                                                "cfi.scaled",
                                                "nnfi.scaled",
                                                "aic",
                                                "bic"))

semPaths(fit_anxapp_RC5,
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
chi_nrc8_anxapp_meas <- lavTestLRT(fit_meas_RC5, fit_anxapp_RC5)
### comparison model does NOT fit the data significantly better
