## Loading and cleaning data
## Author: Ian J. Kahrirlas
# This script loads in and cleans participant pool data for the manuscript Savoring the moment: A link between
# affectivity and depression. 
# Run this script first to reproduce results.

# load in packages
library(tidyverse)
library(readr)


# load in and clean up the Fall 2014 and Spring 2015 data
fall_2014 <- c(list.files(pattern = "Individual Differences"), 
               list.files(pattern = "FINAL SPRING 2015"))

fall_2014 <- map_df(fall_2014, read_csv)

# eliminate extraneous columns
fall_2014 <- fall_2014[-(1:13)]
fall_2014 <- fall_2014[-c(3, 4, 11:15)]

# make new names for demographic information
fall_2014 <- fall_2014 %>% rename("age" = names(fall_2014[1]),
                                  "sex" = names(fall_2014[2]),
                                  "ethnicity" = names(fall_2014[3]),
)

# columns 9 to 30 pertain to the BRIEF questionnaire, so we can label them as such
names(fall_2014)[9:30] <- paste("BRIEF", 1:22, sep = "_")

# eliminate questionnaires not relelvant to the present study hypotheses
fall_2014 <- fall_2014 %>% 
  select(-c(names(fall_2014)[31:102]))

# label columns 31-46 as PSWQ items
names(fall_2014)[31:46] <- paste("PSWQ", 1:16, sep = "_")


# label MASQ items
names(fall_2014)[47:85] <- paste("MASQ", 1:39, sep = "_")

# label PHQ items
names(fall_2014)[86:97] <- paste("PHQ", 1:12, sep = "_")

# eliminate pain measures
fall_2014 <- fall_2014 %>%
  select(-c(98:117))

# label SBI items
names(fall_2014)[98:121] <- paste("SBI", 1:24, sep = "_")

# eliminate the last item
fall_2014 <- fall_2014 %>%
  select(-122)

# Eliminate first row of data, which doesn't contain anything
fall_2014 <- fall_2014[-1,]

# all of the data are of the character class. This can be remedied by using the `map_dbl` and `as.numeric` functions. 
# this combination of functions iterates over each column of a data frame and turns them into numeric class data

# this code illustrates the problem at hand, that all variables are of character class
map_chr(fall_2014, class)
# this changes all the variables to nuermical class
fall_2014 <- map_df(fall_2014, as.numeric)

# fix the race varialbe
names(fall_2014)[4] <- "X19" # to reduce verbose code, rename the variable
fall_2014 <- fall_2014 %>%
  select(X19:X23) %>% # limit scope of operations
  rowwise() %>% # perform vector-wise operations on rows rather than columns
  mutate(race = if_else(sum(X19, X20, X21, X22, X23, na.rm = TRUE) >= 6, "biracial", # if the sum of all the race variables is 6 or greater, participant is biracial
                        if_else(!is.na(X19), "american indian", 
                                if_else(!is.na(X20), "asian",
                                        if_else(!is.na(X21), "pacific islander",
                                                if_else(!is.na(X22), "african american",
                                                        if_else(!is.na(X23), "caucasian", NA_character_)
                                                        )
                                                )
                                        )
                                )
                        )
         ) %>%
  select(race) %>% # only select the newly created race variable
  bind_cols(., fall_2014) %>% # bind everything back to the original dataset
  select(-c(paste0("X", 19:23))) # get rid of the old race variables

# here it is now
select(fall_2014, race)

# turn it into a function
tidy_race <- function(df) {
  names(df)[4] <- "X19"
  df <- df %>%
    select(X19:X23) %>%
    rowwise() %>% #
    mutate(race = if_else(sum(X19, X20, X21, X22, X23, na.rm = TRUE) >= 6, "biracial",
                          if_else(!is.na(X19), "american indian", 
                                  if_else(!is.na(X20), "asian",
                                          if_else(!is.na(X21), "pacific islander",
                                                  if_else(!is.na(X22), "african american",
                                                          if_else(!is.na(X23), "caucasian", NA_character_)
                                                  )
                                          )
                                  )
                          )
    )
    ) %>%
    select(race) %>%
    bind_cols(., df) %>%
    select(-c(paste0("X", 19:23)))
}

# load in and clean Fall 2015 data
fall_2015 <- read_csv(list.files(pattern = "FINAL FALL 2015"))

# Omit unnecessary variables
fall_2015 <- fall_2015[-c(1:13, 16, 17, 24:28, 51:117, 185:204, 229:237)]

# name variables
names(fall_2015)[c(1, 2, 3)] <- c("age", "sex", "ethnicity")
names(fall_2015)[9:30] <- paste("BRIEF", 1:22, sep = "_")
names(fall_2015)[31:46] <- paste("PSWQ", 1:16, sep = "_")
names(fall_2015)[47:85] <- paste("MASQ", 1:39, sep = "_")
names(fall_2015)[86:97] <- paste("PHQ", 1:12, sep = "_")
names(fall_2015)[98:121] <- paste("SBI", 1:24, sep = "_")

# clean data using similar procedures as followed with Fall 2014 and Spring 2015.
fall_2015 <- fall_2015[-1,]
fall_2015 <- map_df(fall_2015, as.numeric)
fall_2015 <- tidy_race(fall_2015)

# define function for cleaning
sad_rename_clean <- function(df) {
  names(df)[c(1, 2, 3)] <- c("age", "sex", "ethnicity")
  names(df)[31:46] <- paste("PSWQ", 1:16, sep = "_")
  names(df)[47:85] <- paste("MASQ", 1:39, sep = "_")
  names(df)[9:30] <- paste("BRIEF", 1:22, sep = "_")
  names(df)[86:97] <- paste("PHQ", 1:12, sep = "_")
  names(df)[98:121] <- paste("SBI", 1:24, sep = "_")
  df <- df[-1,]
  df <- map_df(df, as.numeric)
  df <- tidy_race(df)
}

# load in and clean Spring 2016 data
spring_2016 <- read_csv(list.files(pattern = "FINAL SPRING 2016"))
spring_2016 <- spring_2016[-c(1:13, 16, 17, 24:28, 51:117, 185:204, 229:247)]

# use the sad_rename_clean function to name the columns and clean it up a bit
spring_2016 <- sad_rename_clean(spring_2016)
glimpse(spring_2016)

# load in and clean Fall 2016 data
fall_2016 <- read_csv(list.files(pattern = "FINAL FALL 2016"))

# delete unnecessary columns
fall_2016 <- fall_2016[-c(1:13, 16, 17, 24:28, 51:117, 185:204, 229:246, 301)]

# use functions to clean data
fall_2016 <- sad_rename_clean(fall_2016)

# Rename PANAS questionnaire
# "S_" is for the state PANAS
names(fall_2016)[118:144] <- paste("S_PANAS", 1:27, sep = "_")
# "T_" is for the trait PANAS
names(fall_2016)[145:171] <- paste("T_PANAS", 1:27, sep = "_")

# load in and clean Spring 2017 to Spring 2019 data
data_to_merge <- c("sad_fall_2017.csv",
                   "sad_fall_2018.csv",
                   "sad_spring_2017.csv",
                   "sad_spring_2018.csv",
                   "sad_spring_2019.csv")

spring_2017_spring_2019 <- map_df(data_to_merge, read_csv)

# get rid of extra variables in new data
spring_2017_spring_2019 <- spring_2017_spring_2019 %>% 
  select(-c(color_blind, handedness, native_english:household_member))

# fix race variable
spring_2017_spring_2019 <- spring_2017_spring_2019 %>% 
  select(`American Indian or Alaskan Native`:`White / Caucasian`) %>% 
  rowwise() %>% 
  mutate(race = if_else(sum(`American Indian or Alaskan Native`, 
                            Asian,
                            `Native Hawaiian or Other Pacific Islander`,
                            `Black or African American`,
                            `White / Caucasian`, 
                            na.rm = TRUE) >= 6, "biracial",
                        if_else(!is.na(`American Indian or Alaskan Native`), "american indian",
                                if_else(!is.na(Asian), "asian",
                                        if_else(!is.na(`Native Hawaiian or Other Pacific Islander`), "pacific islander",
                                                if_else(!is.na(`Black or African American`), "african american",
                                                               if_else(!is.na(`White / Caucasian`), "caucasian", NA_character_)
                                                )
                                        )
                                )
                        )
  )
         ) %>% 
  select(race) %>%
  bind_cols(., spring_2017_spring_2019) %>% 
  select(-c(`American Indian or Alaskan Native`:`White / Caucasian`))

# change age to numerical class
spring_2017_spring_2019$age <- as.numeric(spring_2017_spring_2019$age)

# merge all of the data together
data_list <- list(fall_2014,
                  fall_2015,
                  fall_2016,
                  spring_2016,
                  spring_2017_spring_2019)
sad_data <- map_df(data_list, bind_rows)

# check out the data}
glimpse(sad_data)

# omit brief and PANAS
sad_data <- sad_data %>% 
  select(-starts_with("BRIEF")) %>%
  select(-contains("PANAS"))
names(sad_data)

# Let's inspect the age variable for suspect entries
sad_data %>% 
  distinct(age) %>%
  pull(.)

# there are missing values, some who erroneously input their ID numbers as their age, 
# a 15 year old, and someone who entered zero. We will replace these with NAs and omit the 15 year old
sad_data <- sad_data %>% filter(!age == 15)

sad_data <- sad_data %>% mutate(age = na_if(age, 7087106163),
                                age = na_if(age, 6306244561),
                                age = na_if(age, 0))

# r omit all cases with missing values
sad_data <- na.omit(sad_data)

# proof that 0% of all variables are missing
map_dbl(sad_data, ~ {
  (sum(is.na(.x))/sum(!is.na(.x))) * 100
})

# save to workspace for other scripts to work
write_csv(sad_data, "total_data.csv")