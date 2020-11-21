# Reading in data and analyses for Hassan's Mulcahy project

## load packages
library(tidyverse)
library(readxl)
library(MBESS)
library(naniar)
library(janitor)

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
per_dataset <- per_dataset %>%
  mutate(masq_pa = masq_2 + masq_4 + masq_5 + masq_7 + masq_11 + masq_14 + masq_19 + masq_23 + masq_26 + masq_28 + masq_32 + masq_34 + masq_36 + masq_37,
         masq_na = masq_9 + masq_13 + masq_17 + masq_21 + masq_29 + masq_30 + masq_35 + masq_38,
         masq_aa = masq_1 + masq_3 + masq_6 + masq_8 + masq_10 + masq_12 + masq_15 + masq_16 + masq_18 + masq_20 + masq_22 + masq_24 + masq_25 + masq_27 + masq_31 + masq_33 + masq_39
  )
