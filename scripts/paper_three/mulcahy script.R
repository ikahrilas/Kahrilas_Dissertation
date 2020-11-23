# Reading in data and analyses for Hassan's Mulcahy project

## load packages
library(tidyverse)
library(readxl)
library(MBESS)
library(naniar)
library(janitor)
library(DataExplorer)
library(lmerTest)
library(performance)
library(effectsize)

## read in data
dat <- read_excel("data/paper_three/headspace_questionnaire_data.xlsx") %>%
  select("Initials",
         "Cohort",
         "Code Activation Date",
         "Group",
         contains(c("T1", "T2", "T3"))) %>%
  filter(Group != "WL+")

## race code
dat_hs <- dat %>%
  rowwise() %>%
  mutate(race = if_else(sum(`Q9African-American_T1`,
                            Q9AmericanIndianorAlaskaNat_T1,
                            `Q9Asian-American_T1`,
                            Q9NativeHawaiianorOtherPaci_T1,
                            Q9HispanicorLatino_T1,
                            `Q9Non-HispanicCaucasian/Whi_T1`,
                            Q9Other_T1
  ) >=2, "Multi-Racial", if_else(
    `Q9African-American_T1` == 1, "African American", if_else(
      Q9AmericanIndianorAlaskaNat_T1 == 1, "American Indian or Alaska Native", if_else(
        `Q9Asian-American_T1` == 1, "Asian American", if_else(
          Q9NativeHawaiianorOtherPaci_T1 == 1, "Native Hawaiian or Other Pacific Islander", if_else(
            Q9HispanicorLatino_T1 == 1, "Hispanic or Latinx", if_else(
              Q9Other_T1 == 1, Q9OtherOther_T1, "Caucasian"
            )
          )
        )
      )
    )
  )
  ),
  race = if_else(race == "666", "Other", race)
  ) %>%
  rename(pid = StudyID_T1)

dat_hs <- dat_hs %>%
  select("Initials",
         "Cohort",
         "Group",
         "race",
         contains(c("Gender", "age", "PHQ", "MASQ", "SBI", "PSWQ", "PANAS", "BSI")),
         -contains("Language")) %>%
  rename("age" = "Age_T1",
         "gender" = "Gender_T1")

## replace 666 and 999 with NAs
dat_hs <- map_df(dat_hs, ~ {
  na_if(.x, 999) %>%
    na_if(666)
})

## clean up variable names
names(dat_hs) <- str_remove_all(names(dat_hs), "...160")
names(dat_hs) <- str_remove_all(names(dat_hs), "...163")
names(dat_hs) <- str_remove_all(names(dat_hs), "...170")
names(dat_hs) <- str_remove_all(names(dat_hs), "...173")

## convert MASQ8 to numeric
dat_hs$MASQ8_T1 <- as.numeric(dat_hs$MASQ8_T1)
dat_hs$SBI17_T1 <- as.numeric(dat_hs$SBI17_T1)
dat_hs$SBI23_T3 <- as.numeric(dat_hs$SBI23_T3)

## replace initials column with pid
dat_hs <- dat_hs %>%
  select(-Initials) %>%
  mutate(pid = as.character(1:nrow(.))) %>%
  relocate(pid, everything())

## convert dataframe to long form
dat_hs_long <- dat_hs %>%
  pivot_longer(cols = PHQ1_T1:last_col(),
               names_to = c(".value", "time"),
               names_sep = "_"
  )

# names to lower case
names(dat_hs_long) <- tolower(names(dat_hs_long))

# masq subscales
dat_hs_long <- dat_hs_long %>%
  mutate(masq_pa = masq2 + masq4 + masq5 + masq7 + masq11 + masq14 + masq19 + masq23 + masq26 + masq28 + masq32 + masq34 + masq36 + masq37,
         masq_na = masq9 + masq13 + masq17 + masq21 + masq29 + masq30 + masq35 + masq38,
         masq_aa = masq1 + masq3 + masq6 + masq8 + masq10 + masq12 + masq15 + masq16 + masq18 + masq20 + masq22 + masq24 + masq25 + masq27 + masq31 + masq33 + masq39
  )
# phq total scale
dat_hs_long <- dat_hs_long %>%
  mutate(phq_total = phq1 + phq2 + phq3 + phq4 + phq5 + phq6 + phq7 + phq8 + phq9)

# change group variable to just HS and WL
dat_hs_long <- dat_hs_long %>%
  mutate(
    group = case_when(
      group == "PSA" ~ "HS",
      group == "HS No Orientation" ~ "HS",
      TRUE ~ group
    ))
dat_hs_long$group <- relevel(as.factor(dat_hs_long$group), ref = "WL")

# histograms
## depression histogram
ggplot(dat_hs_long, aes(phq_total)) +
  geom_histogram(fill = "white", color = "black") +
  facet_grid(group ~ time)
## masq pa histogram
ggplot(dat_hs_long, aes(masq_pa)) +
  geom_histogram(fill = "white", color = "black") +
  facet_grid(group ~ time)
## masq na histogram
ggplot(dat_hs_long, aes(masq_na)) +
  geom_histogram(fill = "white", color = "black") +
  facet_grid(group ~ time)
## masq aa histogram
ggplot(dat_hs_long, aes(masq_aa)) +
  geom_histogram(fill = "white", color = "black") +
  facet_grid(group ~ time)

# models
## depression as outcome
mod_dep <- lmer(phq_total ~ time*group + (1|pid), data = dat_hs_long)
anova(mod_dep)
summary(mod_dep)
effectsize(mod_dep)
check_model(mod_dep, panel = FALSE)
confint.merMod(mod_dep, parm = "beta_", method = "boot", nsim = 5000)
## PA as outcome
mod_pa <- lmer(masq_pa ~ time*group + (1|pid), data = dat_hs_long)
anova(mod_pa)
summary(mod_pa)
effectsize(mod_pa)
check_model(mod_pa, panel = FALSE)
confint.merMod(mod_pa, parm = "beta_", method = "boot", nsim = 5000)
## NA as outcome
mod_na <- lmer(masq_na ~ time*group + (1|pid), data = dat_hs_long)
anova(mod_na)
summary(mod_na)
effectsize(mod_na)
check_model(mod_na, panel = FALSE)
confint.merMod(mod_na, parm = "beta_", method = "boot", nsim = 5000)
## AA as outcome
mod_aa <- lmer(masq_aa ~ time*group + (1|pid), data = dat_hs_long)
anova(mod_aa)
summary(mod_aa)
effectsize(mod_aa)
check_model(mod_aa, panel = FALSE)
confint.merMod(mod_aa, parm = "beta_", method = "boot", nsim = 5000)
