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
dat <- dat %>%
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
         race = if_else(race == "666", "Other", race))

