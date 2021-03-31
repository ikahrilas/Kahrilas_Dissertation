# analyses

## load packages
library(tidyverse)
library(gt)
library(lavaan)
library(here)
library(semPlot)
library(tidySEM)
library(performance)

## read in data and pivot to wide format
dat <- read_csv(here("data", "paper_three", "dat_for_analyses_2021-03-26.csv")) %>%
  pivot_wider(names_from = block,
              values_from = RC2:RC17) %>%
  relocate(pid:race, RC2_Neg_Watch:RC17_NA) %>%
  select(-c(RC2_NA, RC3_NA, RC5_NA, RC7_NA, RC8_NA, RC17_NA)) %>%
  filter(!is.na(RC2_Neg_Watch))

# note that pswq_total + masq_pa + pss_total is giving awesome results
##########################
###INTERNALIZING MODEL###
##########################
internalizing_mod <- '
# internalizing measurement model
  internalizing =~ masq_aa + masq_pa + pswq_total
'

fit <- cfa(internalizing_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit, fit.measures = TRUE, standardized = TRUE)

#####################
## -- RC2
rc2_meas_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  neur_responsivity ~~ 0*internalizing
'

fit_meas_rc2 <- cfa(rc2_meas_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_meas_rc2, fit.measures = TRUE)

model_performance(fit_meas_rc2, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

semPaths(fit_meas_rc2, "std")

# define eci model with RC2
rc2_eci_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
    internalizing =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  neur_responsivity ~~ 0*internalizing
# residual covariances (constrained to equality)
  masq_pa ~~ e*RC2_Pos_Watch + RC2_Neg_Watch
'

fit_eci_rc2 <- cfa(rc2_eci_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_eci_rc2, fit.measures = TRUE, standardized = TRUE)

model_performance(fit_eci_rc2, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

semPaths(fit_eci_rc2)

## model comparison test with measurement model
lavTestLRT(fit_meas_rc2, rc2_eci_mod)
### comparison model does NOT fit the data significantly better

# define neural reactivity model with RC2
rc2_nr_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  neur_responsivity ~~ 0*internalizing + masq_pa
'

fit_nr_rc2 <- cfa(rc2_nr_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_nr_rc2, fit.measures = TRUE, standardized = TRUE)

## model comparison test with measurement model
lavTestLRT(fit_meas_rc2, fit_nr_rc2)
### comparison model does NOT fit the data significantly better

# define internalizing model with RC2
rc2_int_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  neur_responsivity ~~ 0*internalizing
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
internalizing ~~ e*RC2_Pos_Watch + e*RC2_Neg_Watch
'

fit_int_rc2 <- cfa(rc2_int_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_int_rc2, fit.measures = TRUE, standardized = TRUE)

## model comparison test with measurement model
lavTestLRT(fit_meas_rc2, fit_int_rc2)
### comparison model does NOT fit the data significantly better

################ MODIFICATION INDICES
modindices(fit_meas_rc2, sort = TRUE)

## correlate residual variance of masq_pa and positive RC2
rc2_pos_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC2_Neu_Watch + RC2_Pos_Watch + RC2_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  neur_responsivity ~~ 0*internalizing
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
RC2_Pos_Watch ~~ internalizing
'

fit_pos_rc2 <- sem(rc2_pos_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_pos_rc2, fit.measures = TRUE, standardized = TRUE)

## model comparison test with measurement model
lavTestLRT(fit_meas_rc2, fit_pos_rc2)
### comparison model does NOT fit the data significantly better

modindices(fit_pos_rc2, sort = TRUE)

#######################################
## -- RC3
# define measurement model with RC3
rc3_meas_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  neur_responsivity ~~ 0*internalizing
'

fit_meas_rc3 <- cfa(rc3_meas_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_meas_rc3, fit.measures = TRUE, standardized = TRUE)

semPaths(fit_meas_rc3)

# define eci model with RC3
rc3_eci_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  neur_responsivity ~~ 0*internalizing
# residual covariances (constrained to equality)
  masq_pa ~~ RC3_Pos_Watch + RC3_Neg_Watch
'

fit_eci_rc3 <- cfa(rc3_eci_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_eci_rc3, fit.measures = TRUE, standardized = TRUE)

## model comparison test with measurement model
lavTestLRT(fit_meas_rc3, fit_eci_rc3)
### comparison model does NOT fit the data significantly better

# define neural reactivity model with RC3
rc3_nr_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no association between nr and int factory but parameter estimate between nr and pa
  neur_responsivity ~~ 0*internalizing + masq_pa
'

fit_nr_rc3 <- cfa(rc3_nr_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_nr_rc3, fit.measures = TRUE, standardized = TRUE)

## model comparison test with measurement model
lavTestLRT(fit_meas_rc3, fit_nr_rc3)
### comparison model does NOT fit the data significantly better

# define internalizing model with RC3
rc3_int_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  neur_responsivity ~~ 0*internalizing
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
internalizing ~~ e*RC3_Pos_Watch + e*RC3_Neg_Watch
'

fit_int_rc3 <- cfa(rc3_int_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_int_rc3, fit.measures = TRUE, standardized = TRUE)

## model comparison test with measurement model
lavTestLRT(fit_meas_rc3, fit_int_rc3)
### comparison model does NOT fit the data significantly better

## MODIFICATION INDICES
modindices(fit_meas_rc3, sort = TRUE)


## correlate the factors
rc3_cor_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC3_Neu_Watch + RC3_Pos_Watch + RC3_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
'

fit_cor_rc3 <- sem(rc3_cor_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_cor_rc3, fit.measures = TRUE, standardized = TRUE)

## -- RC8
# define measurement model with RC8
rc8_meas_mod <- '
# neural responsivity factor
  NR =~ RC8_Neu_Watch + RC8_Pos_Watch + RC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
'

fit_meas_rc8 <- cfa(rc8_meas_mod, data = dat, estimator = "MLR", missing = "ML", orthogonal = TRUE)

model_performance(fit_meas_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

summary(fit_meas_rc8, fit.measures = TRUE, standardized = TRUE)

model_performance(fit_meas_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

semPaths(fit_meas_rc8,
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

# define eci model with RC8
rc8_eci_mod <- '
# neural responsivity factor
  NR =~ RC8_Neu_Watch + RC8_Pos_Watch + RC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ pswq_total + masq_pa + masq_aa
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ e*RC8_Pos_Watch + e*RC8_Neg_Watch
'

fit_eci_rc8 <- cfa(rc8_eci_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_eci_rc8, fit.measures = TRUE, standardized = TRUE)

model_performance(fit_eci_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

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
lavTestLRT(fit_meas_rc8, fit_eci_rc8)
### comparison model DOES fit the data significantly better

# define neural reactivity model with RC8
rc8_nr_mod <- '
# neural responsivity factor
  NR =~ RC8_Neu_Watch + RC8_Pos_Watch + RC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no association between nr and int factory but parameter estimate between nr and pa
  NR ~~ 0*INT + masq_pa
'

fit_nr_rc8 <- cfa(rc8_nr_mod, data = dat, estimator = "MLR", missing = "ML")

model_performance(fit_nr_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

summary(fit_nr_rc8, fit.measures = TRUE, standardized = TRUE)

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
lavTestLRT(fit_meas_rc8, fit_nr_rc8)
### comparison model DOES fit the data significantly better

# define internalizing model with RC8
rc8_int_mod <- '
# neural responsivity factor
  NR =~ RC8_Neu_Watch + RC8_Pos_Watch + RC8_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
INT ~~ e*RC8_Pos_Watch + e*RC8_Neg_Watch
'

fit_int_rc8 <- cfa(rc8_int_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_int_rc8, fit.measures = TRUE, standardized = TRUE)

model_performance(fit_int_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

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
lavTestLRT(fit_meas_rc8, fit_int_rc8)
### comparison model does NOT fit the data significantly better

# define anxiety model with RC8
rc8_anx_mod <- '
# neural responsivity factor
  NR =~ RC8_Neu_Watch + RC8_Pos_Watch + RC8_Neg_Watch
# internalNRizing measurement model with MASQ subscales
  INT =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  NR ~~ 0*INT
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
masq_aa ~~ e*RC8_Pos_Watch + e*RC8_Neg_Watch
'

fit_anx_rc8 <- cfa(rc8_anx_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_anx_rc8, fit.measures = TRUE, standardized = TRUE)

model_performance(fit_anx_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

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
lavTestLRT(fit_meas_rc8, fit_anx_rc8)
### comparison model does NOT fit the data significantly better



##### getting at the "brightening effect" a bit more
## no equality constraint between pa and neg/pos
rc8_be_mod <- '
# neural responsivity factor
  NR =~ RC8_Neu_Watch + RC8_Pos_Watch + RC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ pswq_total + masq_pa + masq_aa
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ RC8_Pos_Watch + RC8_Neg_Watch
'

fit_be_rc8 <- cfa(rc8_be_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_be_rc8, fit.measures = TRUE, standardized = TRUE)

model_performance(fit_be_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

semPaths(fit_be_rc8,
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

# compare with eci model
lavTestLRT(fit_be_rc8, fit_eci_rc8)

## just estimate PA and positive, no negative
rc8_bep_mod <- '
# neural responsivity factor
  NR =~ RC8_Neu_Watch + RC8_Pos_Watch + RC8_Neg_Watch
# internalizing measurement model with MASQ subscales
  INT =~ pswq_total + masq_pa + masq_aa
# orthogonal model
  NR ~~ 0*INT
# residual covariances (constrained to equality)
  masq_pa ~~ RC8_Pos_Watch
'

fit_bep_rc8 <- cfa(rc8_bep_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_bep_rc8, fit.measures = TRUE, standardized = TRUE)

model_performance(fit_bep_rc8, metrics = c("RMSEA", "SRMR", "CFI", "NNFI", "AIC", "BIC")) %>% gt()

semPaths(fit_bep_rc8,
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

# compare with be model
lavTestLRT(fit_bep_rc8, fit_be_rc8)

# compare with eci model
lavTestLRT(fit_bep_rc8, fit_eci_rc8)





################################
## RC7 measurement model
rc7_meas_mod <- '
# neural responsivity factor
  NR =~ RC7_Neu_Watch + RC7_Pos_Watch + RC7_Neg_Watch
# internalizing measurement model trio scores
  INT =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  NR ~~ 0*INT
'

fit_meas_rc7 <- cfa(rc7_meas_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_meas_rc7, fit.measures = TRUE, standardized = TRUE)

## RC7 eci model
rc7_eci_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC7_Neu_Watch + RC7_Pos_Watch + RC7_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  neur_responsivity ~~ 0*internalizing
# eci
  masq_pa ~~ e*RC7_Pos_Watch + e*RC7_Neg_Watch
'

fit_eci_rc7 <- cfa(rc7_eci_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_eci_rc7, fit.measures = TRUE, standardized = TRUE)

semPaths(fit_eci_rc7)

# define neural reactivity model with RC7
rc7_nr_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC7_Neu_Watch + RC7_Pos_Watch + RC7_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no association between nr and int factory but parameter estimate between nr and pa
  neur_responsivity ~~ 0*internalizing + masq_pa
'

fit_nr_rc7 <- cfa(rc7_nr_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_nr_rc7, fit.measures = TRUE, standardized = TRUE)

# define internalizing model with RC7
rc7_int_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC7_Neu_Watch + RC7_Pos_Watch + RC7_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  neur_responsivity ~~ 0*internalizing
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
internalizing ~~ e*RC7_Pos_Watch + e*RC7_Neg_Watch
'

fit_int_rc7 <- cfa(rc7_int_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_int_rc7, fit.measures = TRUE, standardized = TRUE)





### MODIFICATION INDICES
modindices(fit_meas_rc7, sort = TRUE)

## estimate covariance between negatie watch and masq aa
rc7_aa_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC7_Neu_Watch + RC7_Pos_Watch + RC7_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  neur_responsivity ~~ 0*internalizing
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
masq_aa ~~ RC7_Neg_Watch
'

fit_aa_rc7 <- sem(rc7_aa_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_aa_rc7, fit.measures = TRUE, standardized = TRUE)

semPaths(fit_aa_rc7)
## model comparison test with measurement model
lavTestLRT(fit_meas_rc7, fit_aa_rc7)
### comparison model does NOT fit the data significantly better


################################
## RC5 measurement model
rc5_meas_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  neur_responsivity ~~ 0*internalizing
'

fit_meas_rc5 <- cfa(rc5_meas_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_meas_rc5, fit.measures = TRUE, standardized = TRUE)

## RC5 eci model
rc5_eci_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# orthogonal model
  neur_responsivity ~~ 0*internalizing
# eci
  masq_pa ~~ e*RC5_Pos_Watch + e*RC5_Neg_Watch
'

fit_eci_rc5 <- cfa(rc5_eci_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_eci_rc5, fit.measures = TRUE, standardized = TRUE)

semPaths(fit_eci_rc5)

# define neural reactivity model with RC7
rc5_nr_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no association between nr and int factory but parameter estimate between nr and pa
  neur_responsivity ~~ 0*internalizing + masq_pa
'

fit_nr_rc5 <- cfa(rc5_nr_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_nr_rc5, fit.measures = TRUE, standardized = TRUE)

# define internalizing model with RC7
rc5_int_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC5_Neu_Watch + RC5_Pos_Watch + RC5_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  neur_responsivity ~~ 0*internalizing
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
internalizing ~~ RC5_Pos_Watch + RC5_Neg_Watch
'

fit_int_rc5 <- cfa(rc5_int_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_int_rc5, fit.measures = TRUE, standardized = TRUE)


### MODIFICATION INDICES
modindices(fit_meas_rc7, sort = TRUE)

## estimate covariance between negatie watch and masq aa
rc7_aa_mod <- '
# neural responsivity factor
  neur_responsivity =~ RC7_Neu_Watch + RC7_Pos_Watch + RC7_Neg_Watch
# internalizing measurement model with MASQ subscales
  internalizing =~ masq_aa + masq_pa + pswq_total
# no correlation with internalizing factor but correlated with masq_pa
  neur_responsivity ~~ 0*internalizing
# estimate parameters among pos and neg watch and internalizing factor constrained to equality
masq_aa ~~ RC7_Neg_Watch
'

fit_aa_rc7 <- sem(rc7_aa_mod, data = dat, estimator = "MLR", missing = "ML")

summary(fit_aa_rc7, fit.measures = TRUE, standardized = TRUE)

semPaths(fit_aa_rc7)
## model comparison test with measurement model
lavTestLRT(fit_meas_rc7, fit_aa_rc7)
### comparison model does NOT fit the data significantly better


