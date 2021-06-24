# results tables for paper three

# load packages
library(kableExtra)
library(tidyverse)
library(lavaan)
library(broom)

# read in sem analyses results
load("data/paper_three/sem_analyses_results.RData")

## -- parameter table -- ##
# p value formatting function
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
  mutate(Path = str_replace(Path, "INT", "PsySx"),
         Path = str_replace(Path, "=~", "→"),
         Path = str_replace(Path, "masq_aa", "Anx Aro"),
         Path = str_replace(Path, "masq_pa", "PA"),
         Path = str_replace(Path, "pswq_total", "Anx App"))

int_fac_table_wide <- bind_cols(int_param_table, # repeat these three times so they can be merged
                                int_param_table, # with the three ERP components
                                int_param_table)

# parameter estimates for erp measurement components
meas_params_list <- list(rc8_meas_params,
                         rc2_meas_params %>% filter(term != "NR ~~ INT"),
                         rc3_meas_params %>% filter(term != "NR ~~ INT"))

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
      mutate(Path = c("Reactivity $\\rightarrow$ Neutral",
                      "Reactivity $\\rightarrow$ Positive",
                      "Reactivity $\\rightarrow$ Negative"))
  })

meas_params_wide <- bind_cols(meas_params_list)

# parameter estimates for erp eci components
eci_params_list <- list(rc8_eci_params %>% filter(term != "NR ~~ INT"),
                        rc2_eci_params %>% filter(term != "NR ~~ INT"),
                        rc3_eci_params %>% filter(term != "NR ~~ INT"))

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

# # parameter estimates for NR components
# nr_params_list <- list(rc8_nr_params,
#                        rc2_nr_params,
#                        rc3_nr_params)
#
# nr_params_list <-
#   map(nr_params_list, ~{
#     .x %>%
#       select(term, estimate, std.all, std.error, p.value) %>%
#       filter(term == "NR ~~ masq_pa") %>%
#       mutate(across(.cols = c(estimate, std.all, std.error),
#                     .fns = ~ sprintf("%.2f", .x))) %>%
#       mutate(p.value = if_else(p.value < .001, "<.001",
#                                num_format(p.value))) %>%
#       unite("Est/Std", estimate, std.all, sep = "/") %>%
#       rename("Path" = "term",
#              "$SE$" = "std.error",
#              "$p$" = "p.value") %>%
#       mutate(Path = "NR $\\leftrightarrow$ PA")
#   })
#
# nr_params_wide <- bind_cols(nr_params_list)

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
      mutate(Path = c("PsySx $\\leftrightarrow$ Positive",
                      "PsySx $\\leftrightarrow$ Negative"))
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
  pack_rows("Psychopathology Symptoms Model", 9, 10) %>%
  pack_rows("Anxious Arousal Model", 11, 12) %>%
  pack_rows("Anxious Apprehension Model", 13, 14) %>%
  row_spec(0, align = "c") %>%
  landscape() %>%
  footnote(general = "Path labels correspond to parameter estimates in each model. Group headings
  in the 'path' column denote the specific model that the following parameter estimates are
  unique to. PsySx = Psychological Symptoms latent factor, ECI = emotion context insensitivity, $\\SE$ = standard error,
  Est/Std =
    \nundstandardized and standardized parameter estimate, $\\\\rightarrow$ = latent factor loading,
  $\\\\leftrightarrow$ = covariance.",
           threeparttable = TRUE,
           escape = FALSE,
           general_title = "Note.",
           footnote_as_chunk = TRUE) %>%
  save_kable("images/paper_3/sem_param_table.pdf")

## -- fit index table -- ##
# define function that inputs chi square difference test output is input and outputs
# latex formatted test statistic and p value
chisq_format_fn <- function(chi){
  p <- chi$`Pr(>Chisq)`[2]
  if (p < .001) {
    p <- "< .001"
  } else {
    p <- paste("=", sub("^(-?)0.", "\\1.", sprintf("%.3f", p)))
  }
  df_diff <- chi$`Df diff`[2]
  chi_diff <- sprintf("%.2f", chi$`Chisq diff`[2])
  paste0("$\\Delta$", "${\\chi}^2$","(", df_diff, ") ", "= ", chi_diff, ", ", "$p$ ", p)
}

# table for rc8 fit indices
rc8_fit_table <-
  bind_rows(rc8_meas_fit,
            rc8_eci_fit,
            nrc8_int_fit,
            rc8_anx_fit,
            rc8_anxapp_fit) %>%
    mutate("Model" = c("Measurement",
                       "Emotion Context Insensitivity",
                       "Psychopathology Symptoms",
                       "Anxious Arousal",
                       "Anxious Apprehension")) %>%
    relocate(Model, everything()) %>%
    mutate("$\\Delta{\\chi}^2$" =
             c("Nested Measurement Model",
               chisq_format_fn(chi_nrc8_eci_meas),
               chisq_format_fn(chi_nrc8_int_meas),
               chisq_format_fn(chi_nrc8_anx_meas),
               chisq_format_fn(chi_nrc8_anxapp_meas))
           )

# table for rc2 fit indices
rc2_fit_table <-
  bind_rows(rc2_meas_fit,
            rc2_eci_fit,
            rc2_int_fit,
            rc2_anx_fit,
            rc2_anxapp_fit) %>%
  mutate("Model" = c("Measurement",
                     "Emotion Context Insensitivity",
                     "Psychopathology Symptoms",
                     "Anxious Arousal",
                     "Anxious Apprehension")) %>%
  relocate(Model, everything()) %>%
  mutate("$\\Delta{\\chi}^2$" =
           c("Nested Measurement Model",
             chisq_format_fn(chi_rc2_eci_meas),
             chisq_format_fn(chisq_rc2_int_meas),
             chisq_format_fn(chi_rc2_anx_meas),
             chisq_format_fn(chi_nrc8_anxapp_meas))
        )

# table for rc3 fit indices
rc3_fit_table <-
  bind_rows(rc3_meas_fit,
            rc3_eci_fit,
            rc3_int_fit,
            rc3_anx_fit,
            rc3_anxapp_fit) %>%
  mutate("Model" = c("Measurement",
                     "Emotion Context Insensitivity",
                     "Psychopathology Symptoms",
                     "Anxious Arousal",
                     "Anxious Apprehension")) %>%
  relocate(Model, everything()) %>%
  mutate("$\\Delta{\\chi}^2$" =
           c("Nested Measurement Model",
             chisq_format_fn(chisq_rc3_eci_meas),
             chisq_format_fn(chisq_rc3_int_meas),
             chisq_format_fn(chisq_rc3_anx_meas),
             chisq_format_fn(chi_rc3_anxapp_meas))
  )

# put all the tables together
fit_tab <-
  bind_rows(
    rc8_fit_table,
    rc2_fit_table,
    rc3_fit_table
  ) %>%
    mutate(nnfi.scaled = if_else(nnfi.scaled > 1.0, "1.00", sprintf("%.2f", nnfi.scaled)),
           across(.cols = c(chisq.scaled,  # round to two decimal places
                            rmsea.scaled,
                            srmr,
                            cfi.scaled,
                            aic,
                            bic),
                  .fns = ~ sprintf("%.2f", .x))) %>%
           rename("${\\chi}^2$" = "chisq.scaled",
                  "$df$" = "df",
                  "RMSEA" = "rmsea.scaled",
                  "SRMR" = "srmr",
                  "CFI" = "cfi.scaled",
                  "NNFI" = "nnfi.scaled",
                  "AIC" = "aic",
                  "BIC" = "bic")

# construct the table
kable(fit_tab, "latex", escape = FALSE, booktabs = TRUE,
      align = c("l", "r", "c", "c", "c", "c", "c", "c", "c", "l"), linesep = "") %>%
  kable_styling(font_size = 12, latex_options = c("scale_down")) %>%
  pack_rows("257 ms Component", 1, 5) %>%
  pack_rows("371 ms Component", 6, 10) %>%
  pack_rows("736 ms Component", 11, 15) %>%
  row_spec(0, align = c("l", rep("c", times = 9))) %>%
  landscape() %>%
  footnote(general = "${\\\\chi}^2$ = chi-square, $df$ = degrees of freedom, RMSEA = root mean square error of
           approximation, SRMR = standardized root mean square residual, CFI = comparative fit index,
           NNFI = non-normed fit index, AIC = Akaike Information Criteria, BIC = Bayesian Information Criteria,
           $\\\\Delta{\\\\chi}^2$ = Satorra-Bentler scaled difference chi-square test comparing each model with its \\\\newline respective
           nested measurement model.",
           threeparttable = TRUE,
           escape = FALSE,
           general_title = "Note.",
           footnote_as_chunk = TRUE) %>%
  save_kable("images/paper_3/sem_fit_table.pdf")

