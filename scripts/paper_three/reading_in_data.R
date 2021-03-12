# reading in headspace questionnaire data

## load packages
library(tidyverse)
library(readxl)
library(MBESS)
library(naniar)
library(haven)

## read in data
dat <-
  read_sav("data/paper_three/questionnaire/SMILE full data 9 2 2020.sav") %>%
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
         contains(c("PHQ", "SBI", "MASQ", "PSWQ", "PANAS")))

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
  mutate(group = "headspace")

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
         masq_aa = masq1 + masq3 + masq6 + masq8 + masq10 + masq12 + masq15 + masq16 + masq18 + masq20 + masq22 + masq24 + masq25 + masq27 + masq31 + masq33 + masq39
  )

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

##-- read in PER data for harmonization
per_dat <- read_csv(here::here("data", "paper_two", "created_data", "per_measures.csv"))

## subtract 9 from total phq score to make consistent with hs data and remove duplicate cases
per_dat <-
  per_dat %>%
    mutate(depression = depression - 9) %>%
    rename(phq_total = depression) %>%
    select(-c(block, valence, difficulty, arousal)) %>%
    distinct() %>%
  mutate(group = "per")

## trim down both datasets and merge them
dat_hs <-
  dat_hs %>%
  select(pid, gender, age, race, phq_total, pswq_total, panas_pos_total, panas_neg_total,
         anticipating:masq_aa, group)

per_dat <-
  per_dat %>%
  mutate(sex = if_else(sex == 1, "Male", if_else(sex == 2, "Female", "Other")),
         pswq_total = pswq_1 + pswq_2 + pswq_3 + pswq_4 + pswq_5 +
           pswq_6 + pswq_7 + pswq_8 + pswq_9 + pswq_10 + pswq_11 +
           pswq_12 + pswq_13 + pswq_14 + pswq_15 + pswq_16) %>%
   select(pid, age, sex, Race, phq_total, pswq_total, pos_affectivity, neg_affectivity, masq_pa, masq_na, masq_aa,
          anticipating, savoring_moment, reminiscing, sbi_tot, group) %>%
  rename(gender =  sex,
         panas_pos_total = pos_affectivity,
         panas_neg_total = neg_affectivity,
         race = Race)

total_dat <- bind_rows(dat_hs, per_dat)
# write to workspace
write_csv(total_dat, file = here::here("data", "paper_three", "total_questionnaire_data.csv"))

## -- Code for deriving omega reliability coefficients
#### anticipating
ci.reliability(sbi_aim_one %>%
                 select(sbi_1, sbi_7, sbi_13, sbi_19, sbi_4, sbi_10, sbi_16, sbi_22),
               type = "omega")

#### savoring the moment
ci.reliability(sbi_aim_one %>%
                 select(sbi_5 ,sbi_11, sbi_17, sbi_23, sbi_2, sbi_8, sbi_14, sbi_20),
               type = "omega")
#### reminiscing
ci.reliability(sbi_aim_one %>%
                 select(sbi_3, sbi_9, sbi_15, sbi_21, sbi_6, sbi_12, sbi_18, sbi_24),
               type = "omega")

### depression
dep_aim_one <- dat_hs %>%
  select(contains("phq") & contains("t1"))

names(dep_aim_one) <- paste0("phq_", 1:9)

tmp <- per_dataset %>%
  select(contains("phq")) %>%
  select(-c("phq_10", "phq_11", "phq_12"))

dep_aim_one <- bind_rows(dep_aim_one, tmp)

ci.reliability(dep_aim_one, type = "omega")

### PA
masq_pa_items <- c("masq_2", "masq_4", "masq_5", "masq_7", "masq_11", "masq_14", "masq_19", "masq_23", "masq_26", "masq_28", "masq_32", "masq_34", "masq_36", "masq_37")
masq_na_items <- c("masq_9", "masq_13", "masq_17", "masq_21", "masq_29", "masq_30", "masq_35", "masq_38")
masq_aa_items <- c("masq_1", "masq_3", "masq_6", "masq_8", "masq_10", "masq_12", "masq_15", "masq_16", "masq_18", "masq_20", "masq_22", "masq_24", "masq_25", "masq_27", "masq_31", "masq_33", "masq_39")

pa_aim_one <- dat_hs %>%
  select(contains("masq") & contains("t1"))

names(pa_aim_one) <- paste0("masq_", 1:39)

pa_aim_one <- pa_aim_one %>% select(masq_pa_items)

tmp <- per_dataset %>% select(masq_pa_items)

pa_aim_one <- bind_rows(pa_aim_one, tmp)

ci.reliability(pa_aim_one, type = "omega")

### NA
na_aim_one <- dat_hs %>%
  select(contains("masq") & contains("t1"))

names(na_aim_one) <- paste0("masq_", 1:39)

na_aim_one <- na_aim_one %>% select(masq_na_items)

tmp <- per_dataset %>% select(masq_na_items)

na_aim_one <- bind_rows(na_aim_one, tmp)

ci.reliability(na_aim_one, type = "omega")

### AA
aa_aim_one <- dat_hs %>%
  select(contains("masq") & contains("t1"))

names(aa_aim_one) <- paste0("masq_", 1:39)

aa_aim_one <- aa_aim_one %>%
  select(masq_aa_items) %>%
  mutate(masq_8 = as.numeric(masq_8))

tmp <- per_dataset %>% select(masq_aa_items)

aa_aim_one <- bind_rows(aa_aim_one, tmp)

ci.reliability(aa_aim_one, type = "omega")







## create smaller dataset of smile data that only contains demographic information
hs_dem <- dat_hs %>%
  select(pid, Group, SexualOrientation_T1, Gender_T1, race, ParentalHouseholdIncome_T1, Age_T1) %>%
  rename(household_income = ParentalHouseholdIncome_T1,
         age = Age_T1,
         sex = Gender_T1,
         orientation = SexualOrientation_T1) %>%
  mutate(race = str_to_title(race))

## demographics for manuscript
### age
range(hs_dem$age)
mean(hs_dem$age)
sd(hs_dem$age)

### gender
hs_dem %>%
  group_by(sex) %>%
  summarize(n = n(), perc = (n() / nrow(.)) * 100)

### sexual orientation
hs_dem %>%
  group_by(orientation) %>%
  summarize(n = n(), perc = (n() / nrow(.)) * 100)

### race
hs_dem %>%
  group_by(race) %>%
  summarize(n = n(), perc = (n() / nrow(.)) * 100)

### ses
hs_dem %>%
  group_by(household_income) %>%
  summarize(n = n(), perc = (n() / nrow(.)) * 100)

# read in PER data for harmonization
per_dataset <- read_xlsx("data/paper_two/questionnaire/PER Questionnaires.xlsx")

## without running this code, R is unable to recognize some column names
names(per_dataset) <- enc2native(names(per_dataset))

## rename some of the race categories to match the smile data
per_dataset <- per_dataset %>%
  mutate(Race = if_else(Race == "caucasian" & ethnicity == 1, "Hispanic or Latinx", Race),
         Race = if_else(Race == "multiracial", "Multi-Racial", Race),
         Race = if_else(Race == "african_american", "African American", Race),
         Race = str_to_title(Race)) %>%
  rename(race = Race)

## create smaller dataset of per data that only contains deompgraphic information
per_dem <- per_dataset %>%
  select(pid, age, race, sex, household_income) %>%
  mutate(sex = if_else(sex == 1, "Male", if_else(
    sex == 2, "Female", "Other"
  )),
  household_income = if_else(household_income == 1, "Under $10,000", if_else(
    household_income == 2, "$10,000-$19,999", if_else(
      household_income == 3, "$20,000-$29,999", if_else(
        household_income == 4, "$30,000-$49,999", if_else(
          household_income == 5, "$50,000-$69,999", if_else(
            household_income == 6, "$70,000-$89,999", if_else(
              household_income == 7, "$90,000-$109,999", if_else(
                household_income == 8, "$110,000-$129,999", "Over $130,000"
              )
            )
            )
          )
        )
      )
    )
  ))

## demographics for SES
per_dem %>%
  group_by(household_income) %>%
  summarize(total = n(),
            perc = (total / 52) * 100)

# vector containing EEG pid at time one
eeg_pid <- c(22585442,
             22585444,
             22585436,
             22585437,
             22585441,
             22585446,
             22585439,
             22585448,
             22585452,
             22585443,
             22585440,
             22585450,
             22585449,
             22585438,
             22585464,
             22585462,
             22585466,
             22585469,
             22585470,
             22585458,
             22586479,
             22585460,
             22585471,
             22585474,
             22585473,
             22585465,
             22585472,
             22585480,
             22585475,
             22585467,
             22585478,
             22585491,
             22585482,
             22585496,
             22585490,
             22585484,
             22585488,
             22585500,
             22585499,
             22585498,
             22585483,
             22585492,
             22585501,
             22585502,
             22585503,
             22585487,
             22585511,
             22585513,
             22585512,
             22585510,
             22585515,
             22585520,
             22585519,
             22585508,
             22585524,
             22585517,
             22585522,
             22585545,
             22585553,
             22585541,
             22585548,
             22585549,
             22585551,
             22585560,
             22585543,
             22585565,
             22585569,
             22585564,
             22585568,
             22585570,
             22585563,
             22585573,
             22585572,
             22585577)

# vector of eeg pids at time 2
t2 <- c(22585437,
        22585440,
        22585444,
        22585441,
        22585448,
        22585439,
        22585449,
        22585446,
        22585452,
        22585450,
        22585438,
        22585443,
        22585442,
        22585436,
        22585470,
        22585471,
        22585473,
        22586479,
        22585458,
        22585474,
        22585472,
        22585460,
        22585462,
        22585469,
        22585465,
        22585475,
        22585467,
        22585480,
        22585478,
        22585501,
        22585503,
        22585499,
        22585484,
        22585490,
        22585492,
        22585500,
        22585483,
        22585491,
        22585488,
        22585520,
        22585526,
        22585524,
        22585515,
        22585511,
        22585519,
        22585508,
        22585517,
        22585522,
        22585513,
        22585549,
        22585541,
        22585545,
        22585548,
        22585558,
        22585551,
        22585552,
        22585544,
        22585560,
        22585543,
        22585561)

# vector with pids that have t1 and t1 EEG
eeg_pid_comp <- eeg_pid[eeg_pid %in% t2]

# data frame of headspace data with only EEG participants
hs_comp_eeg <- hs_dem %>%
  filter(pid %in% eeg_pid_comp)

## demographic data for manuscript
### age
range(hs_comp_eeg$age)
mean(hs_comp_eeg$age)
sd(hs_comp_eeg$age)

### gender
hs_comp_eeg %>%
  group_by(sex) %>%
  summarize(n = n(), perc = (n() / nrow(.)) * 100)

### sexual orientation
hs_comp_eeg %>%
  group_by(orientation) %>%
  summarize(n = n(), perc = (n() / nrow(.)) * 100)

### race
hs_comp_eeg %>%
  group_by(race) %>%
  summarize(n = n(), perc = (n() / nrow(.)) * 100)

### ses
hs_comp_eeg %>%
  group_by(household_income) %>%
  summarize(n = n(), perc = (n() / nrow(.)) * 100)


# reliability coefficients for manuscript
## aim one
### sbi
sbi_aim_one <- dat_hs %>%
  select(contains("sb") & contains("t1"))

names(sbi_aim_one) <- paste0("sbi_", 1:24)

sbi_aim_one <- sbi_aim_one %>%
  mutate(sbi_17 = as.numeric(sbi_17))

##### reverse score SBI items. Items are a 7-item likert scale, so subtract each score from 8.
sbi_rev_items <- paste("sbi_", seq(2, 24, by = 2), sep = "") #variable containing SBI items to be reversed, which is every other item from 2 to 24

tmp <- per_dataset %>%
  select(contains("sbi"))

sbi_aim_one <- bind_rows(sbi_aim_one, tmp)

sbi_aim_one <- sbi_aim_one %>%
  mutate(across(.cols = sbi_rev_items, .fns = ~ 8 - .x))

#### anticipating
ci.reliability(sbi_aim_one %>%
                 select(sbi_1, sbi_7, sbi_13, sbi_19, sbi_4, sbi_10, sbi_16, sbi_22),
               type = "omega")

#### savoring the moment
ci.reliability(sbi_aim_one %>%
                 select(sbi_5 ,sbi_11, sbi_17, sbi_23, sbi_2, sbi_8, sbi_14, sbi_20),
               type = "omega")
#### reminiscing
ci.reliability(sbi_aim_one %>%
                 select(sbi_3, sbi_9, sbi_15, sbi_21, sbi_6, sbi_12, sbi_18, sbi_24),
               type = "omega")

### depression
dep_aim_one <- dat_hs %>%
  select(contains("phq") & contains("t1"))

names(dep_aim_one) <- paste0("phq_", 1:9)

tmp <- per_dataset %>%
  select(contains("phq")) %>%
  select(-c("phq_10", "phq_11", "phq_12"))

dep_aim_one <- bind_rows(dep_aim_one, tmp)

ci.reliability(dep_aim_one, type = "omega")

### PA
masq_pa_items <- c("masq_2", "masq_4", "masq_5", "masq_7", "masq_11", "masq_14", "masq_19", "masq_23", "masq_26", "masq_28", "masq_32", "masq_34", "masq_36", "masq_37")
masq_na_items <- c("masq_9", "masq_13", "masq_17", "masq_21", "masq_29", "masq_30", "masq_35", "masq_38")
masq_aa_items <- c("masq_1", "masq_3", "masq_6", "masq_8", "masq_10", "masq_12", "masq_15", "masq_16", "masq_18", "masq_20", "masq_22", "masq_24", "masq_25", "masq_27", "masq_31", "masq_33", "masq_39")

pa_aim_one <- dat_hs %>%
  select(contains("masq") & contains("t1"))

names(pa_aim_one) <- paste0("masq_", 1:39)

pa_aim_one <- pa_aim_one %>% select(masq_pa_items)

tmp <- per_dataset %>% select(masq_pa_items)

pa_aim_one <- bind_rows(pa_aim_one, tmp)

ci.reliability(pa_aim_one, type = "omega")

### NA
na_aim_one <- dat_hs %>%
  select(contains("masq") & contains("t1"))

names(na_aim_one) <- paste0("masq_", 1:39)

na_aim_one <- na_aim_one %>% select(masq_na_items)

tmp <- per_dataset %>% select(masq_na_items)

na_aim_one <- bind_rows(na_aim_one, tmp)

ci.reliability(na_aim_one, type = "omega")

### AA
aa_aim_one <- dat_hs %>%
  select(contains("masq") & contains("t1"))

names(aa_aim_one) <- paste0("masq_", 1:39)

aa_aim_one <- aa_aim_one %>%
  select(masq_aa_items) %>%
  mutate(masq_8 = as.numeric(masq_8))

tmp <- per_dataset %>% select(masq_aa_items)

aa_aim_one <- bind_rows(aa_aim_one, tmp)

ci.reliability(aa_aim_one, type = "omega")

## aim two
### sbi
sbi_aim_two_t1 <- dat_hs %>%
  select(contains("sb") & contains("t1"))

sbi_aim_two_t2 <- dat_hs %>%
  select(contains("sb") & contains("t2"))

names(sbi_aim_two_t1) <- paste0("sbi_", 1:24)

names(sbi_aim_two_t2) <- paste0("sbi_", 1:24)

sbi_aim_two_t1 <- sbi_aim_two_t1 %>%
  mutate(sbi_17 = as.numeric(sbi_17))

sbi_aim_two_t2 <- sbi_aim_two_t2 %>%
  mutate(sbi_17 = as.numeric(sbi_17))

##### reverse score SBI items. Items are a 7-item likert scale, so subtract each score from 8.
sbi_rev_items <- paste("sbi_", seq(2, 24, by = 2), sep = "") #variable containing SBI items to be reversed, which is every other item from 2 to 24

sbi_aim_two_t1 <- sbi_aim_two_t1 %>%
  mutate(across(.cols = sbi_rev_items, .fns = ~ 8 - .x))

sbi_aim_two_t2 <- sbi_aim_two_t2 %>%
  mutate(across(.cols = sbi_rev_items, .fns = ~ 8 - .x))

#### anticipating time 1
ci.reliability(sbi_aim_two_t1 %>%
                 select(sbi_1, sbi_7, sbi_13, sbi_19, sbi_4, sbi_10, sbi_16, sbi_22),
               type = "omega")

#### anticipating time 2
ci.reliability(sbi_aim_two_t2 %>%
                 select(sbi_1, sbi_7, sbi_13, sbi_19, sbi_4, sbi_10, sbi_16, sbi_22),
               type = "omega")

#### savoring the moment time 1
ci.reliability(sbi_aim_two_t1 %>%
                 select(sbi_5 ,sbi_11, sbi_17, sbi_23, sbi_2, sbi_8, sbi_14, sbi_20),
               type = "omega")

#### savoring the moment time 2
ci.reliability(sbi_aim_two_t2 %>%
                 select(sbi_5 ,sbi_11, sbi_17, sbi_23, sbi_2, sbi_8, sbi_14, sbi_20),
               type = "omega")


#### reminiscing time 1
ci.reliability(sbi_aim_two_t1 %>%
                 select(sbi_3, sbi_9, sbi_15, sbi_21, sbi_6, sbi_12, sbi_18, sbi_24),
               type = "omega")

#### reminiscing time 2
ci.reliability(sbi_aim_two_t2 %>%
                 select(sbi_3, sbi_9, sbi_15, sbi_21, sbi_6, sbi_12, sbi_18, sbi_24),
               type = "omega")

### depression
dep_aim_two_t1 <- dat_hs %>%
  select(contains("phq") & contains("t1"))

dep_aim_two_t2 <- dat_hs %>%
  select(contains("phq") & contains("t2"))

names(dep_aim_two_t1) <- paste0("phq_", 1:9)

names(dep_aim_two_t2) <- paste0("phq_", 1:9)

#### time one
ci.reliability(dep_aim_two_t1, type = "omega")

#### time two
ci.reliability(dep_aim_two_t2, type = "omega")


### PA
masq_pa_items <- c("masq_2", "masq_4", "masq_5", "masq_7", "masq_11", "masq_14", "masq_19", "masq_23", "masq_26", "masq_28", "masq_32", "masq_34", "masq_36", "masq_37")
masq_na_items <- c("masq_9", "masq_13", "masq_17", "masq_21", "masq_29", "masq_30", "masq_35", "masq_38")
masq_aa_items <- c("masq_1", "masq_3", "masq_6", "masq_8", "masq_10", "masq_12", "masq_15", "masq_16", "masq_18", "masq_20", "masq_22", "masq_24", "masq_25", "masq_27", "masq_31", "masq_33", "masq_39")

pa_aim_two_t1 <- dat_hs %>%
  select(contains("masq") & contains("t1"))

pa_aim_two_t2 <- dat_hs %>%
  select(contains("masq") & contains("t2"))

names(pa_aim_two_t1) <- paste0("masq_", 1:39)

names(pa_aim_two_t2) <- paste0("masq_", 1:39)

pa_aim_two_t1 <- pa_aim_two_t1 %>% select(masq_pa_items)

pa_aim_two_t2 <- pa_aim_two_t2 %>% select(masq_pa_items)

#### time one
ci.reliability(pa_aim_two_t1, type = "omega")

#### time two
ci.reliability(pa_aim_two_t2, type = "omega")

### NA
na_aim_two_t1 <- dat_hs %>%
  select(contains("masq") & contains("t1"))

na_aim_two_t2 <- dat_hs %>%
  select(contains("masq") & contains("t2"))

names(na_aim_two_t1) <- paste0("masq_", 1:39)

names(na_aim_two_t2) <- paste0("masq_", 1:39)

na_aim_two_t1 <- na_aim_two_t1 %>% select(masq_na_items)

na_aim_two_t2 <- na_aim_two_t2 %>% select(masq_na_items)

#### time one
ci.reliability(na_aim_two_t1, type = "omega")

#### time two
ci.reliability(na_aim_two_t2, type = "omega")

### AA
aa_aim_two_t1 <- dat_hs %>%
  select(contains("masq") & contains("t1"))

aa_aim_two_t2 <- dat_hs %>%
  select(contains("masq") & contains("t2"))

names(aa_aim_two_t1) <- paste0("masq_", 1:39)

names(aa_aim_two_t2) <- paste0("masq_", 1:39)

aa_aim_two_t1 <- aa_aim_two_t1 %>% select(masq_aa_items)

aa_aim_two_t2 <- aa_aim_two_t2 %>% select(masq_aa_items)

#### time one
ci.reliability(aa_aim_two_t1, type = "omega")

#### time two
ci.reliability(aa_aim_two_t2, type = "omega")

