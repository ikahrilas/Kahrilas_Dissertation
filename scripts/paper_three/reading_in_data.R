# reading in headspace questionnaire data

## load packages
library(tidyverse)
library(readxl)
library(MBESS)
library(naniar)
library(haven)
library(data.table)

## read in data
dat <-
  read_sav("data/paper_three/questionnaire/SMILE full data 11 10 2020 FOR DATA CLEANING.sav") %>%
  select("Initials",
         "Cohort",
         "Group",
         contains(c("T1", "T2", "T3"))) %>%
  filter(Group != "WL+") %>%
  select("StudyID_T1",
         "Group",
         Gender_T1,
         Age_T1,
         SexualOrientation_T1,
         Q9AfricanAmerican_T1:Q9OtherOther_T1,
         Religion_T1,
         FreeText4_T1,
         ParentalHouseholdIncome_T1,
         contains(c("PHQ", "SBI", "MASQ", "PSWQ", "PANAS", "PSS")))

## race code
dat_hs <-
  dat %>%
  rowwise() %>%
  mutate(race = if_else(sum(`Q9AfricanAmerican_T1`,
                            Q9AmericanIndianorAlaskaNat_T1,
                            `Q9AsianAmerican_T1`,
                            Q9NativeHawaiianorOtherPaci_T1,
                            Q9HispanicorLatino_T1,
                            `Q9NonHispanicCaucasianWhi_T1`,
                            Q9Other_T1
                            ) >=2, "Multi-Racial", if_else(
                              `Q9AfricanAmerican_T1` == 1, "African American", if_else(
                                Q9AmericanIndianorAlaskaNat_T1 == 1, "American Indian or Alaska Native", if_else(
                                  `Q9AsianAmerican_T1` == 1, "Asian American", if_else(
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

dat_hs <- dat_hs %>% replace_with_na_all(condition = ~.x == 666)
dat_hs <- dat_hs %>% replace_with_na_all(condition = ~.x == 999)

## for paper 3, just retain time one data
dat_hs <-
  dat_hs %>%
    select(pid,
           contains("T1"),
           race,
           -c(FreeText4_T1, Q9AfricanAmerican_T1:Q9OtherOther_T1))

names(dat_hs) <- str_remove(tolower(names(dat_hs)), "_t1")

## create group variable to denote which study participants are from
dat_hs <-
  dat_hs %>%
  mutate(group = "hs")

####################
# derive subscales #
####################

## -- SBI
## get rid of variables brynn made
dat_hs <-
  dat_hs %>%
  select(-c(sbi_mean:sbi_remin_total,
            masq_anx_mean:masq_dep_total))

### reverse score SBI items. Items are a 7-item likert scale, so subtract each score from 8.
sbi_rev_items <- paste0("sbi", seq(2, 24, by = 2)) #variable containing SBI items to be reversed, which is every other item from 2 to 24

dat_hs <-
  dat_hs %>%
  mutate(across(.cols = c(sbi_rev_items),
                .fns = ~ 8 - .x,
                .names = "{.col}_r"))

dat_hs <-
  dat_hs %>%
  mutate(anticipating = (sbi1 + sbi7 + sbi13 + sbi19 + sbi4_r + sbi10_r + sbi16_r + sbi22_r) / 8,
         savoring_moment = (sbi5 + sbi11 + sbi17 + sbi23 + sbi2_r + sbi8_r + sbi14_r + sbi20_r) / 8,
         reminiscing  = (sbi3 + sbi9 + sbi15 + sbi21 + sbi6_r + sbi12_r + sbi18_r + sbi24_r) / 8,
         sbi_tot = (sbi1 + sbi7 + sbi13 + sbi19 + sbi4_r + sbi10_r + sbi16_r + sbi22_r + sbi5 + sbi11 + sbi17 +
                      sbi23 + sbi2_r + sbi8_r + sbi14_r + sbi20_r + sbi3 + sbi9 + sbi15 + sbi21 + sbi6_r +
                      sbi12_r + sbi18_r + sbi24_r) / 24)

## -- MASQ
dat_hs <-
  dat_hs %>%
  mutate(masq_pa = masq2 + masq4 + masq5 + masq7 + masq11 + masq14 + masq19 + masq23 + masq26 + masq28 + masq32 + masq34 + masq36 + masq37,
         masq_na = masq9 + masq13 + masq17 + masq21 + masq29 + masq30 + masq35 + masq38,
         masq_aa = masq1 + masq3 + masq6 + masq8 + masq10 + masq12 + masq15 + masq16 + masq18 + masq20 + masq22 + masq24 + masq25 + masq27 + masq31 + masq33 + masq39,
         # note that positive items are subtracted from six
         masq_ad = (6 - masq2) + (6 - masq4) + (6 - masq5) + (6 - masq7) + (6 - masq11) + (6 - masq14) + (6 - masq19) + (6 - masq23) + (6 - masq26) + (6 - masq28) + (6 - masq32) + (6 - masq34) + (6 - masq36) + (6 - masq37) + masq9 + masq13 + masq17 + masq21 + masq29 + masq30 + masq35 + masq38,
         )

## -- PHQ-9
dat_hs <-
  dat_hs %>%
  mutate(phq_total =
           (phq1 - 1) +
           (phq2 - 1) +
           (phq3 - 1) +
           (phq4 - 1) +
           (phq5 - 1) +
           (phq6 - 1) +
           (phq7 - 1) +
           (phq8 - 1) +
           (phq9 - 1))
glimpse(dat_hs)

## -- PSWQ
dat_hs <-
  dat_hs %>%
  mutate(pswq_total = (6 - pswq1) + pswq2 + (6 - pswq3) +
           pswq4 + pswq5 + pswq6 + pswq7 + (6 - pswq8) +
           pswq9 + (6 - pswq10) + (6 - pswq11) + pswq12 +
           pswq13 + pswq14 + pswq15 + pswq16)

## -- PSS
dat_hs <-
  dat_hs %>%
  mutate(pss_total = (pss1 - 1) + (4 - (pss2 - 1)) + (4 - (pss3 - 1)) + (pss4 - 1)) %>% glimpse()

## -- PANAS
dat_hs <-
  dat_hs %>%
    mutate(panas_pos_total = panas1 + panas3 + panas5 + panas9 + panas10 + panas12 + panas14 +
             panas16 + panas17 + panas19,
           panas_neg_total = panas2 + panas4 + panas6 + panas7 + panas8 + panas11 + panas13 + panas15 +
             panas18 + panas20)

##-- read in PER data for harmonization
per_dat <- read_csv(here::here("data", "paper_two", "created_data", "per_measures.csv"))

## -- RESCORE 27-ITEM PANAS TO 20-ITEM
per_dat <-
  per_dat %>%
  mutate(panas_pos_total = t_panas_16 + t_panas_5 + t_panas_23 + t_panas_11 + t_panas_15 +
                           t_panas_24 + t_panas_9 + t_panas_14 + t_panas_20 + t_panas_19,
         panas_neg_total = t_panas_13 + t_panas_8 + t_panas_2 + t_panas_27 + t_panas_6 +
                         t_panas_25 + t_panas_26 + t_panas_4 + t_panas_7 + t_panas_3)

# rename variables to match per dataset
setnames(per_dat,
         old = c("t_panas_16", "t_panas_5", "t_panas_23", "t_panas_11", "t_panas_15",
                  "t_panas_24", "t_panas_9", "t_panas_14", "t_panas_20", "t_panas_19"),
         new = c("panas1", "panas3", "panas5", "panas9", "panas10", "panas12", "panas14",
                 "panas16", "panas17", "panas19"))

setnames(per_dat,
         old = c("t_panas_13",  "t_panas_8",  "t_panas_2",  "t_panas_27",  "t_panas_6",
                   "t_panas_25",  "t_panas_26",  "t_panas_4",  "t_panas_7",  "t_panas_3"),
         new = c("panas2", "panas4", "panas6", "panas7", "panas8", "panas11", "panas13", "panas15",
                   "panas18", "panas20"))

## -- SCORE 14-ITEM PSS TO 4-ITEM
per_dat <-
  per_dat %>%
  mutate(pss_total = (pss_2 - 1) + (4 - (pss_6 - 1)) + (4 - (pss_7 - 1)) + (pss_14 - 1))

## rename pss items to match hs
setnames(per_dat,
         old = c("pss_2",  "pss_6",  "pss_7",  "pss_14"),
         new = c("pss1", "pss2", "pss3", "pss4"))

## -- PHQ9
# Subtract 9 from total phq score to make consistent with hs data and remove duplicate cases
per_dat <-
  per_dat %>%
    mutate(depression = depression - 9) %>%
    rename(phq_total = depression) %>%
    select(-c(block, valence, difficulty, arousal)) %>%
    distinct() %>%
  mutate(group = "per")

## -- PSWQ
per_dat <-
  per_dat %>%
  rename_with(~ str_remove(.x, "_"), .cols = c(pswq_1:pswq_16))

per_dat <-
  per_dat %>%
  mutate(pswq_total = (6 - pswq1) + pswq2 + (6 - pswq3) +
           pswq4 + pswq5 + pswq6 + pswq7 + (6 - pswq8) +
           pswq9 + (6 - pswq10) + (6 - pswq11) + pswq12 +
           pswq13 + pswq14 + pswq15 + pswq16)

## -- MASQ AD SCALE
per_dat <-
  per_dat %>%
  mutate(masq_ad = (6 - masq_2) + (6 - masq_4) + (6 - masq_5) + (6 - masq_7) +
           (6 - masq_11) + (6 - masq_14) + (6 - masq_19) + (6 - masq_23) +
           (6 - masq_26) + (6 - masq_28) + (6 - masq_32) + (6 - masq_34) +
           (6 - masq_36) + (6 - masq_37) + masq_9 + masq_13 + masq_17 + masq_21 +
           masq_29 + masq_30 + masq_35 + masq_38)

## merge the datasets
dat_hs <-
  dat_hs %>%
  select(pid, gender, age, race,
         contains("phq"),
         contains("masq"),
         contains("pswq"),
         contains("panas"),
         contains("pss"),
         anticipating,
         reminiscing,
         savoring_moment,
         group) %>%
    select(-c(paste0("phq", 1:9, "r")),
           -c(paste0("masq", c(2, 4, 5, 7, 11, 14, 19, 23, 26, 28, 32, 34, 36, 37), "r")),
           -c(paste0("pswq", c(1, 3, 8, 10, 11), "r")),
           -c(paste0("pss", c(1, 2, 3, 4), "r")),
           -c(paste0("pss", c(2, 3), "rr")),
           -c(phq_mean, pss_mean, phq9t1, masq_mean, masq_total, pswq_mean, pswqt1, panas_mean, panas_total,
              panasnt1, panas_pos_mean, panas_neg_mean, panaspt1)) %>% glimpse()

per_dat <-
  per_dat %>%
  mutate(sex = if_else(sex == 1, "Male", if_else(sex == 2, "Female", "Other")),
         group = "per") %>%
   select(pid, age, sex, Race,
          contains("phq"),
          -c(phq_10, phq_11, phq_12),
          contains("pswq"),
          contains("panas"),
          contains("pss"),
          -contains("t_panas"),
          -contains("s_panas"),
          -c("pss_1", "pss_3", "pss_4", "pss_5", "pss_8", "pss_9", "pss_10", "pss_11", "pss_12", "pss_13"),
          contains("masq"),
          anticipating,
          savoring_moment,
          reminiscing,
          pss_total,
          group) %>%
  rename(gender = sex,
         race = Race) %>%
  rename_with(~ str_remove(.x, "_"), .cols = c(phq_1:phq_9, masq_1:masq_39)) %>% glimpse()

sum(!(names(per_dat) %in% names(dat_hs))) # all column names match

# merge together
total_dat <- bind_rows(dat_hs, per_dat)

# write to workspace
write_csv(total_dat, file = here::here("data", "paper_three", "total_questionnaire_data.csv"))
