# Reading in data and analyses for Hassan's Mulcahy project

## load packages
library(tidyverse)
library(MBESS)
library(haven)
library(lmerTest)
library(performance)
library(here)
library(effectsize)
library(boot)
library(parameters)

## read in data
dat_hs <-
  read_sav("data/paper_three/SMILE full data 9 2 2020.sav") %>%
  select(StudyID_T1,
         Group,
         Group_R,
         contains("MASQ"),
         contains("PHQ"),
         contains("SBI"))

names(dat_hs) <- tolower(names(dat_hs))

dat_hs$group[dat_hs$group == "PSA"] <- "HS"
dat_hs <- dat_hs[!dat_hs$group == "WL+", ]
dat_hs$group[dat_hs$group == "HSN"] <- "HS"

dat_hs <- dat_hs %>%
  rename("pid" = "studyid_t1")

dat_hs <- dat_hs %>%
  mutate(masq_pa_t1 = masq2_t1 + masq4_t1 + masq5_t1 + masq7_t1 + masq11_t1 + masq14_t1 + masq19_t1 + masq23_t1 + masq26_t1 + masq28_t1 + masq32_t1 + masq34_t1 + masq36_t1 + masq37_t1,
         masq_na_t1 = masq9_t1 + masq13_t1 + masq17_t1 + masq21_t1 + masq29_t1 + masq30_t1 + masq35_t1 + masq38_t1,
         masq_aa_t1 = masq1_t1 + masq3_t1+ masq6_t1 + masq8_t1 + masq10_t1 + masq12_t1 + masq15_t1 + masq16_t1 + masq18_t1 + masq20_t1 + masq22_t1 + masq24_t1 + masq25_t1 + masq27_t1 + masq31_t1 + masq33_t1 + masq39_t1,
         masq_pa_t3 = masq2_t3 + masq4_t3 + masq5_t3 + masq7_t3 + masq11_t3 + masq14_t3 + masq19_t3 + masq23_t3 + masq26_t3 + masq28_t3 + masq32_t3 + masq34_t3 + masq36_t3 + masq37_t3,
         masq_na_t3 = masq9_t3 + masq13_t3 + masq17_t3 + masq21_t3 + masq29_t3 + masq30_t3 + masq35_t3 + masq38_t3,
         masq_aa_t3 = masq1_t3 + masq3_t3 + masq6_t3 + masq8_t3 + masq10_t3 + masq12_t3 + masq15_t3 + masq16_t3 + masq18_t3 + masq20_t3 + masq22_t3 + masq24_t3 + masq25_t3 + masq27_t3 + masq31_t3 + masq33_t3 + masq39_t3,
         phq_t1 = phq1_t1 + phq2_t1 + phq3_t1 + phq4_t1 + phq5_t1 + phq6_t1 + phq7_t1 + phq8_t1 + phq9_t1,
         phq_t3 = phq1_t3 + phq2_t3 + phq3_t3 + phq4_t3 + phq5_t3 + phq6_t3 + phq7_t3 + phq8_t3 + phq9_t3,
         phq_diff = phq_t3 - phq_t1,
         pa_diff = masq_pa_t3 - masq_pa_t1,
         na_diff = masq_na_t3 - masq_na_t1,
         aa_diff = masq_aa_t3 - masq_aa_t1,
         )

# models

## -- DEPRESSION
## Bootstrap 95% CI for beta coefficient
## function to obtain R-Squared from the data
beta <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data = d)
  return(summary(fit)$coef[2])
}

## bootstrapping with 10000 replications
results <- boot(data = dat_hs,
                statistic = beta,
                R = 10000,
                formula = phq_diff ~ group)

lm(phq_diff ~ group, data = dat_hs) %>%
  bootstrap_parameters(iterations = 10000)

## view results
results
plot(results)

## get 95% confidence interval
boot.ci(results, type = "bca")

##-- PA
## bootstrapping with 10000 replications
results <- boot(data = dat_hs,
                statistic = beta,
                R = 10000,
                formula = pa_diff ~ group)

## view results
results
plot(results)

## get 95% confidence interval
boot.ci(results, type = "bca")


## -- NA
## bootstrapping with 10000 replications
results <- boot(data = dat_hs,
                statistic = beta,
                R = 10000,
                formula = na_diff ~ group)

## view results
results
plot(results)

## get 95% confidence interval
boot.ci(results, type = "bca")



# ## read in data
# dat_hs <- read_csv("data/paper_three/smile_questionnaire_data.csv")
# nrow(dat_hs)
# dat_hs <-
#   dat_hs %>%
#   select(StudyID_T1,
#          Group,
#          Group_R,
#          contains("MASQ"),
#          contains("SBI"),
#          contains("PHQ"))
#
# names(dat_hs) <- tolower(names(dat_hs))
#
# dat_hs <- dat_hs %>%
#   mutate(across(masq1_t1:phq9_t7, .fns = ~ as.numeric(.x))) %>%
#   na_if(999) %>%
#   na_if(666)
#
# dat_hs$group[dat_hs$group == "PSA"] <- "HS"
# dat_hs$group[dat_hs$group == "HS No Orientation"] <- "HS"
# dat_hs <- dat_hs[!dat_hs$group == "WL+",]
#
# dat_hs <-
#   dat_hs %>%
#   rename("pid" = studyid_t1) %>%
#   mutate(masq_pa_t1 = masq2_t1 + masq4_t1 + masq5_t1 + masq7_t1 + masq11_t1 + masq14_t1 + masq19_t1 + masq23_t1 + masq26_t1 + masq28_t1 + masq32_t1 + masq34_t1 + masq36_t1 + masq37_t1,
#          masq_na_t1 = masq9_t1 + masq13_t1 + masq17_t1 + masq21_t1 + masq29_t1 + masq30_t1 + masq35_t1 + masq38_t1,
#          masq_aa_t1 = masq1_t1 + masq3_t1+ masq6_t1 + masq8_t1 + masq10_t1 + masq12_t1 + masq15_t1 + masq16_t1 + masq18_t1 + masq20_t1 + masq22_t1 + masq24_t1 + masq25_t1 + masq27_t1 + masq31_t1 + masq33_t1 + masq39_t1,
#          masq_pa_t3 = masq2_t3 + masq4_t3 + masq5_t3 + masq7_t3 + masq11_t3 + masq14_t3 + masq19_t3 + masq23_t3 + masq26_t3 + masq28_t3 + masq32_t3 + masq34_t3 + masq36_t3 + masq37_t3,
#          masq_na_t3 = masq9_t3 + masq13_t3 + masq17_t3 + masq21_t3 + masq29_t3 + masq30_t3 + masq35_t3 + masq38_t3,
#          masq_aa_t3 = masq1_t3 + masq3_t3 + masq6_t3 + masq8_t3 + masq10_t3 + masq12_t3 + masq15_t3 + masq16_t3 + masq18_t3 + masq20_t3 + masq22_t3 + masq24_t3 + masq25_t3 + masq27_t3 + masq31_t3 + masq33_t3 + masq39_t3,
#          phq_t1 = phq1_t1 + phq2_t1 + phq3_t1 + phq4_t1 + phq5_t1 + phq6_t1 + phq7_t1 + phq8_t1 + phq9_t1,
#          phq_t3 = phq1_t3 + phq2_t3 + phq3_t3 + phq4_t3 + phq5_t3 + phq6_t3 + phq7_t3 + phq8_t3 + phq9_t3,
#          pa_diff = masq_pa_t3 - masq_pa_t1,
#          na_diff = masq_na_t3 - masq_na_t1,
#          aa_diff = masq_aa_t3 - masq_aa_t1,
#          phq_diff = phq_t3 - phq_t1
#         )
#
# dat_hs <- dat_hs %>%
#   arrange(pid) %>%
#   filter(!is.na(pid))
#
# dat_ar <- dat_ar %>%
#   arrange(pid) %>%
#   filter(!is.na(pid))
#
# view(dat_hs)
# view(dat_ar)
#
# dat_hs %>%
#   full_join(dat_ar, by = "pid") %>%
#   relocate(phq1_t1.x, phq1_t1.y) %>%
#   view()
#
#
# # analyses
# dep_mod <- lm(phq_diff ~ group, data = dat_hs)
# summary(dep_mod)
#
#
# full_join(dat_ar, dat_hs, by = "pid") %>%
#   select(group.x, group.y, pid) %>%
#   arrange(pid) %>%
#   view()
#
# my_dat <- dat_ar %>%
#   select(pid, group) %>%
#   arrange(pid) %>%
#   pull()
#
# ar_dat <- dat_hs %>%
#   select(pid, group) %>%
#   arrange(pid) %>%
#   pull()
#
#
# ## race code
# dat_hs <-
#   dat_hs %>%
#     rowwise() %>%
#     mutate(race = if_else(sum(`Q9African-American_T1`,
#                             Q9AmericanIndianorAlaskaNat_T1,
#                             `Q9Asian-American_T1`,
#                             Q9NativeHawaiianorOtherPaci_T1,
#                             Q9HispanicorLatino_T1,
#                             `Q9Non-HispanicCaucasian/Whi_T1`,
#                             Q9Other_T1
#     ) >= 2, "Multi-Racial", if_else(
#       `Q9African-American_T1` == 1, "African American", if_else(
#         Q9AmericanIndianorAlaskaNat_T1 == 1, "American Indian or Alaska Native", if_else(
#           `Q9Asian-American_T1` == 1, "Asian American", if_else(
#             Q9NativeHawaiianorOtherPaci_T1 == 1, "Native Hawaiian or Other Pacific Islander", if_else(
#               Q9HispanicorLatino_T1 == 1, "Hispanic or Latinx", if_else(
#                 Q9Other_T1 == 1, Q9OtherOther_T1, "Caucasian"
#               )
#             )
#           )
#         )
#       )
#     )
#     ),
#     race = if_else(race == "666", "Other", race)
#     ) %>%
#     rename(pid = StudyID_T1)
#
# dat_hs <- dat_hs %>%
#   select("pid",
#          "Cohort",
#          "Group",
#          "race",
#          contains(c("Gender", "age", "PHQ", "MASQ", "SBI")),
#          -contains("Language")) %>%
#   rename("age" = "Age_T1",
#          "gender" = "Gender_T1")
#
# ## replace 666 and 999 with NAs
# dat_hs <- map_df(dat_hs, ~ {
#   na_if(.x, 999) %>%
#     na_if(666)
# })
#
# # convert select variables to numeric - they're character for some reason
# dat_hs$MASQ8_T1 <- as.numeric(dat_hs$MASQ8_T1)
# dat_hs$SBI17_T1 <- as.numeric(dat_hs$SBI17_T1)
# dat_hs$SBI23_T3 <- as.numeric(dat_hs$SBI23_T3)
#
# ## clean up variable names
# # names(dat_hs) <- str_remove_all(names(dat_hs), "...160")
# # names(dat_hs) <- str_remove_all(names(dat_hs), "...163")
# # names(dat_hs) <- str_remove_all(names(dat_hs), "...170")
# # names(dat_hs) <- str_remove_all(names(dat_hs), "...173")
#
# # names to lower case
# names(dat_hs) <- tolower(names(dat_hs))
#
# # create subscales and contrast scores
# dat_hs <- dat_hs %>%
#   mutate(masq_pa_t1 = masq2_t1 + masq4_t1 + masq5_t1 + masq7_t1 + masq11_t1 + masq14_t1 + masq19_t1 + masq23_t1 + masq26_t1 + masq28_t1 + masq32_t1 + masq34_t1 + masq36_t1 + masq37_t1,
#          masq_na_t1 = masq9_t1 + masq13_t1 + masq17_t1 + masq21_t1 + masq29_t1 + masq30_t1 + masq35_t1 + masq38_t1,
#          masq_aa_t1 = masq1_t1 + masq3_t1+ masq6_t1 + masq8_t1 + masq10_t1 + masq12_t1 + masq15_t1 + masq16_t1 + masq18_t1 + masq20_t1 + masq22_t1 + masq24_t1 + masq25_t1 + masq27_t1 + masq31_t1 + masq33_t1 + masq39_t1,
#          masq_pa_t3 = masq2_t3 + masq4_t3 + masq5_t3 + masq7_t3 + masq11_t3 + masq14_t3 + masq19_t3 + masq23_t3 + masq26_t3 + masq28_t3 + masq32_t3 + masq34_t3 + masq36_t3 + masq37_t3,
#          masq_na_t3 = masq9_t3 + masq13_t3 + masq17_t3 + masq21_t3 + masq29_t3 + masq30_t3 + masq35_t3 + masq38_t3,
#          masq_aa_t3 = masq1_t3 + masq3_t3 + masq6_t3 + masq8_t3 + masq10_t3 + masq12_t3 + masq15_t3 + masq16_t3 + masq18_t3 + masq20_t3 + masq22_t3 + masq24_t3 + masq25_t3 + masq27_t3 + masq31_t3 + masq33_t3 + masq39_t3,
#          phq_t1 = phq1_t1 + phq2_t1 + phq3_t1 + phq4_t1 + phq5_t1 + phq6_t1 + phq7_t1 + phq8_t1 + phq9_t1,
#          phq_t3 = phq1_t3 + phq2_t3 + phq3_t3 + phq4_t3 + phq5_t3 + phq6_t3 + phq7_t3 + phq8_t3 + phq9_t3,
#          pa_diff = masq_pa_t3 - masq_pa_t1,
#          na_diff = masq_na_t3 - masq_na_t1,
#          aa_diff = masq_aa_t3 - masq_aa_t1,
#          phq_diff = phq_t3 - phq_t1
#   )
#
# # combine headspace groups
# dat_hs$group[dat_hs$group == "PSA"] <- "HS"
# dat_hs$group[dat_hs$group == "HS No Orientation"] <- NA
#
# # analyses
# ## depression
# dep_mod <- lm(phq_diff ~ group, data = dat_hs)
# summary(dep_mod)
#
# summary(dat_ar %>%
#           select(phqdiff,
#                  pa_diff,
#                  na_diff,
#                  aa_diff))
#
# summary(dat_hs %>%
#           select(phq_diff,
#                  pa_diff,
#                  na_diff,
#                  aa_diff))
# table(dat_hs$group)
#
# table(dat_ar$group_rr)
#
#
# dep_mod_ar <- lm(phqdiff ~ group_rr, data = dat_ar)
# summary(dep_mod_ar)
#
#
#
# ### MLM analyses if you wanna
# ## convert dataframe to long form
# dat_hs_long <- dat_hs %>%
#   pivot_longer(cols = PHQ1_T1:last_col(),
#                names_to = c(".value", "time"),
#                names_sep = "_"
#   )
#
# # names to lower case
# names(dat_hs_long) <- tolower(names(dat_hs_long))
#
# # masq subscales
# dat_hs_long <- dat_hs_long %>%
#   mutate(masq_pa = masq2 + masq4 + masq5 + masq7 + masq11 + masq14 + masq19 + masq23 + masq26 + masq28 + masq32 + masq34 + masq36 + masq37,
#          masq_na = masq9 + masq13 + masq17 + masq21 + masq29 + masq30 + masq35 + masq38,
#          masq_aa = masq1 + masq3 + masq6 + masq8 + masq10 + masq12 + masq15 + masq16 + masq18 + masq20 + masq22 + masq24 + masq25 + masq27 + masq31 + masq33 + masq39
#   )
# # phq total scale
# dat_hs_long <- dat_hs_long %>%
#   mutate(phq_total = phq1 + phq2 + phq3 + phq4 + phq5 + phq6 + phq7 + phq8 + phq9)
#
# # change group variable to just HS and WL
# dat_hs_long <- dat_hs_long %>%
#   mutate(
#     group = case_when(
#       group == "HS" ~ "HS",
#       group == "HS No Orientation" ~ "HS",
#       group == "PSA" ~ "HS",
#       group == "WL" ~ "WL"
#     ))
#
# dat_hs_long$group <- relevel(as.factor(dat_hs_long$group), ref = "WL")
#
# # # histograms
# # ## depression histogram
# # ggplot(dat_hs_long, aes(phq_total)) +
# #   geom_histogram(fill = "white", color = "black") +
# #   facet_grid(group ~ time)
# # ## masq pa histogram
# # ggplot(dat_hs_long, aes(masq_pa)) +
# #   geom_histogram(fill = "white", color = "black") +
# #   facet_grid(group ~ time)
# # ## masq na histogram
# # ggplot(dat_hs_long, aes(masq_na)) +
# #   geom_histogram(fill = "white", color = "black") +
# #   facet_grid(group ~ time)
# # ## masq aa histogram
# # ggplot(dat_hs_long, aes(masq_aa)) +
# #   geom_histogram(fill = "white", color = "black") +
# #   facet_grid(group ~ time)
#
# dat <- dat_hs_long %>%
#   select(pid:time, masq_pa:phq_total) %>%
#   pivot_wider(names_from = "time",
#               values_from = masq_pa:phq_total) %>%
#   mutate(pa_diff = masq_pa_T3 - masq_pa_T1,
#          na_diff = masq_na_T3 - masq_na_T1,
#          aa_diff = masq_aa_T3 - masq_aa_T1,
#          phq_diff = phq_total_T3 - phq_total_T1) %>%
#   relocate(pid:age, pa_diff:phq_diff)
#
# # models
# ## phq as outcome
# phq_mod <- lm(phq_diff ~ group, data = dat)
# summary(phq_mod)
#
#





## depression as outcome
# mod_dep <- lmer(phq_total ~ time*group + (1|pid), data = dat_hs_long)
# anova(mod_dep)
# summary(mod_dep)
# effectsize(mod_dep)
# check_model(mod_dep, panel = FALSE)
# confint.merMod(mod_dep, parm = "beta_", method = "boot", nsim = 5000)
# ## PA as outcome
# mod_pa <- lmer(masq_pa ~ time*group + (1|pid), data = dat_hs_long)
# anova(mod_pa)
# summary(mod_pa)
# effectsize(mod_pa)
# check_model(mod_pa, panel = FALSE)
# confint.merMod(mod_pa, parm = "beta_", method = "boot", nsim = 5000)
# ## NA as outcome
# mod_na <- lmer(masq_na ~ time*group + (1|pid), data = dat_hs_long)
# anova(mod_na)
# summary(mod_na)
# effectsize(mod_na)
# check_model(mod_na, panel = FALSE)
# confint.merMod(mod_na, parm = "beta_", method = "boot", nsim = 5000)
# ## AA as outcome
# mod_aa <- lmer(masq_aa ~ time*group + (1|pid), data = dat_hs_long)
# anova(mod_aa)
# summary(mod_aa)
# effectsize(mod_aa)
# check_model(mod_aa, panel = FALSE)
# confint.merMod(mod_aa, parm = "beta_", method = "boot", nsim = 5000)
dat_hs$group <- droplevels(as.factor(dat_hs$group))
dat_hs$group <- relevel(dat_hs$group, "WL")
dat_hs <-
  dat_hs %>%
    filter(group %in% c("WL", "HS")) %>%
  group_by(group) %>%
  mutate(mean_phq = mean(phq_diff))
# Visualizations for results
# define ggplot function
ggplot(dat_hs, aes(group, phq_diff)) +
  geom_violin(fill = "gold", width = 0.5) +
  geom_jitter(color = "maroon", width = 0.1) +
  geom_smooth(aes(as.numeric(group)), method = "lm", se = FALSE, color = "black") +
  theme_classic()

## plot for depression
dep <- plot_fun(dat_hs, aes(group, phq_diff), "Depression")
ggsave(filename = here("images", "paper_3", "mulcahy_poster_images", "dep.png"), plot = dep, device = "png", width = 14, height = 5)
## plot for pa
pa <- plot_fun(dat_hs_long, aes(time, masq_pa), "Positive Affectivity")
ggsave(filename = here("images", "paper_3", "mulcahy_poster_images", "pa.png"), plot = pa, device = "png", width = 14, height = 5)
## plot for na
na <- plot_fun(dat_hs_long, aes(time, masq_na), "Negative Affectivity")
ggsave(filename = here("images", "paper_3", "mulcahy_poster_images", "na.png"), plot = na, device = "png", width = 14, height = 5)


