# reading in headspace questionnaire data

## load packages
library(tidyverse)
library(readxl)

## read in data
dat <- read_excel("data/paper_three/headspace_questionnaire_data.xlsx") %>%
  select("Initials",
         "Cohort",
         "Code Activation Date",
         "Group",
         contains(c("T1", "T2", "T3")))

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

# read in PER data for harmonization
per_dataset <- read_xlsx("data/paper_two/questionnaire/PER Questionnaires.xlsx")

# without running this code, R is unable to recognize some column names
names(per_dataset) <- enc2native(names(per_dataset))

# rename some of the race categories to match the smile data
per_dataset <- per_dataset %>%
  mutate(Race = if_else(Race == "caucasian" & ethnicity == 1, "Hispanic or Latinx", Race),
         Race = if_else(Race == "multiracial", "Multi-Racial", Race),
         Race = if_else(Race == "african_american", "African American", Race),
         Race = str_to_title(Race)) %>%
  rename(race = Race)

# create smaller dataset of smile data that only contains demographic information
hs_dem <- dat_hs %>%
  select(pid, Gender_T1, race, ParentalHouseholdIncome_T1, Age_T1) %>%
  rename(household_income = ParentalHouseholdIncome_T1,
         age = Age_T1,
         sex = Gender_T1)

# create smaller dataset of per data that only contains deompgraphic information
per_dem <- per_dataset %>%
  select(pid, age, race, sex, household_income) %>%
  mutate(sex = if_else(sex == 1, "Male", if_else(
    sex == 2, "Female", "Other"
  )),
  household_income = if_else(household_income == 1, "Under $25,000", if_else(
    household_income == 2, "$25,000-$50,000", if_else(
      household_income == 3, "$50,000-$75,000", if_else(
        household_income == 4, "$75,000-$100,000", if_else(
          household_income == 5, "$100,000-$150,000", if_else(
            household_income == 6, "$150,000-$200,000", "Over $200,000"
            )
          )
        )
      )
    )
  )) %>% View()

