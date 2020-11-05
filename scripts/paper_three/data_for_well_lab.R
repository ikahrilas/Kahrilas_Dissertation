# output Headspace data with variables of interest for Well Lab

## load packages
library(tidyverse)
library(readxl)
library(MBESS)
library(naniar)

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
         contains(c("PHQ", "MASQ", "SBI", "PSWQ", "PANAS", "BSI")))

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
dat$MASQ8_T1 <- as.numeric(dat$MASQ8_T1)
