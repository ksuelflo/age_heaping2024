library(survival)
library(dplyr)
library(flexsurv)
library(knitr)
library(pssst)
library(logquad5q0)
library(SUMMER)
library(tidyr)
library(stringr)
library(readstata13)
library(rdhs)
library(foreach)
library(doMC)
library(ggplot2)
library(gghighlight)
library(pssst)
registerDoMC(cores = 6)

#Getting country codes and country names for DHS countries
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
str(ids)

#Getting all surveys beginning in year 2005 (any survey that starts from 2005 or later), only DHS surveys (not AIDs or malaria, etc.) and 
#Only sub-Saharan African countries. 
survs <- dhs_surveys(surveyYearStart = 2005, surveyType = "DHS", countryIds = c("AO", "BJ", "BT", "BF", "BU", "CM", "CF", "TD", "CG", "CD", "CI", 
                                                                       "ER", "KM", "SZ", "ET", "GA", "GM", "GH", "GN", "KE", "LS", "LB", 
                                                                       "MD", "MW", "ML", "MR", "MZ", "NM", "NI", "NG", "OS", "RW", "SN",
                                                                       "SL", "ZA", "SD", "ST", "TZ", "TG", "UG", "ZM", "ZW"))
#Getting dataset info: Filtering so that it only includes Birth recode files as that is the only thing we want.
datasets <- dhs_datasets(surveyIds = survs$SurveyId, 
                         fileFormat = "flat")%>%
  filter(FileType == "Births Recode")

#Logging in (effectively) to DHS
set_rdhs_config(email = "ksuelflo@macalester.edu",
                project = "Evaluation of Existing Statistical Methods Regarding Age Heaping")

# download datasets. This puts them in a cache somewhere (don't quite understand this part). `downloads` is a list of filepaths corresponding
# to each RDS file (df).

downloads <- get_datasets(datasets$FileName)

#Setting up containers
years_df <- rep(0, length(downloads))

#skipping 47 for now: It has an error with both v023 and v024 strata. And skipping 53, 83

downloads[47]
downloads[53]
downloads[83]

for (i in c(47,53,83)){

  #Reading in data
  print(i)
  births <- readRDS(downloads[[i]])
  
  #Survey year
  year <- min(births$v007)
  years_df[i] <- year
  
  tryCatch(
    expr = {
      df_birth <- format_dhs(df = births, survey_year = year, period_boundaries = c(year-5, year-4, year-3, year-2, year-1, year))
      print(str_c("finished number ", i, " download.", " Using strata(v023)"))
    },
    error = function(e){
      df_birth <- format_dhs(df = births, survey_year = year, period_boundaries = c(year-5, year-4, year-3, year-2, year-1, year), strata = c("v024", "v025"))
      print(str_c("finished number ", i, " download.", " Using strata(v024, v025)"))
    }
  )
  df_container[[i]] <- df_birth
  
}


#Using `ids` and `downloads` to scrape out the country name for each survey.

countries <- rep("", length(downloads))
#All 90 codes for each df
ssa_codes <- str_sub(names(downloads), 1, 2)

#Looping through each survey, and finding corresponding country for the ith code.
for (i in seq_along(countries)){
  code <- ssa_codes[i]
  countries[i] <- ids$CountryName[ids$DHS_CountryCode == code]
}

#Container for summary deaths by month for each survey
summary_deaths <- vector(mode = "list", length = length(df_container))

#This for loop wrangles data into monthly deaths by period. Included are all intervals of death.
for (i in c(47,53,83)){
  #Checking dfs for NULL: There are many that had error with `format_dhs`, so they are null. We will skip those.
  if (is.null(df_container[[i]])){
    next
  }
  
  print(i)
  summary_heap_df <- df_container[[i]]%>%
    filter(I_i ==1)%>%
    mutate(avg_death = (t_0i + t_1i)/ 2,
           death_date = avg_death*30 + b_i,
           died_in_period = if_else(y_p <= death_date & death_date < (y_p + 365), 1, 0))%>%
    group_by(p)%>%
    summarize(int_0_1 = sum(0 <= avg_death & avg_death < 1 & died_in_period == 1),
              int_1_2 = sum(1 <= avg_death & avg_death < 2 & died_in_period == 1),
              int_2_3 = sum(2 <= avg_death & avg_death < 3 & died_in_period == 1),
              int_3_4 = sum(3 <= avg_death & avg_death < 4 & died_in_period == 1),
              int_4_5 = sum(4 <= avg_death & avg_death < 5 & died_in_period == 1),
              int_5_6 = sum(5 <= avg_death & avg_death < 6 & died_in_period == 1),
              int_6_7 = sum(6 <= avg_death & avg_death < 7 & died_in_period == 1),
              int_7_8 = sum(7 <= avg_death & avg_death < 8 & died_in_period == 1),
              int_8_9 = sum(8 <= avg_death & avg_death < 9 & died_in_period == 1),
              int_9_10 = sum(9 <= avg_death & avg_death < 10 & died_in_period == 1),
              int_10_11 = sum(10 <= avg_death & avg_death < 11 & died_in_period == 1),
              int_11_12 = sum(11 <= avg_death & avg_death < 12 & died_in_period == 1),
              int_12_13 = sum(12 <= avg_death & avg_death < 13 & died_in_period == 1),
              int_13_14 = sum(13 <= avg_death & avg_death < 14 & died_in_period == 1),
              int_14_15 = sum(14 <= avg_death & avg_death < 15 & died_in_period == 1),
              int_15_16 = sum(15 <= avg_death & avg_death < 16 & died_in_period == 1),
              int_16_17 = sum(16 <= avg_death & avg_death < 17 & died_in_period == 1),
              int_17_18 = sum(17 <= avg_death & avg_death < 18 & died_in_period == 1),
              int_18_19 = sum(18 <= avg_death & avg_death < 19 & died_in_period == 1),
              int_19_20 = sum(19 <= avg_death & avg_death < 20 & died_in_period == 1),
              int_20_21 = sum(20 <= avg_death & avg_death < 21 & died_in_period == 1),
              int_21_22 = sum(21 <= avg_death & avg_death < 22 & died_in_period == 1),
              int_22_23 = sum(22 <= avg_death & avg_death < 23 & died_in_period == 1),
              int_23_24 = sum(23 <= avg_death & avg_death < 24 & died_in_period == 1),
              int_24_36 = sum(24 <= avg_death & avg_death < 36 & died_in_period == 1),
              int_36_48 = sum(36 <= avg_death & avg_death < 48 & died_in_period == 1),
              int_48_60 = sum(48 <= avg_death & avg_death < 60 & died_in_period == 1),
              int_60_inf = sum(60 <= avg_death & avg_death < Inf & died_in_period == 1))%>%
    mutate(country = countries[i], #Use `i` instead when in the for loop.
           survey_year = years_df[i], 
           .before = p)%>% #Use `i` instead when in the for loop.
    mutate(p_begin = survey_year - 6 + p,
           p_end = p_begin + 1,
           .after = p)
  
  summary_deaths[[i]] <- summary_heap_df
}

saveRDS(summary_deaths, file = "../output/summary.rds")

summary_deaths <- readRDS(file = "../output/summary.rds")
#Now to calculate heaping indexes. Some info about heaping index

#Hill and Choi 2006 index:
# Note that they used 5-9 days, 7 days as the heaped day. This was for early neonatal deaths. However, we are
# extrapolating this method to use at 12 months. 
# Divide the deaths at 12 by 1/5 the total deaths from ages 10-14 months. Can experiment with different ranges here. 

#Pullum and Becker 2014 index:
# Calculate the average deaths from 10-14 months. Then the relative difference between the actual number of deaths
# at 12 months and the average, times 100, is the index. For example, if the average deaths was 50, and the actual
# number was 150, then then index would be 200 (200% higher).

#loop that calculated all the heap indexes for each survey.

#container for new dfs

#START HERE 
# summary_deaths <- readRDS(file = "../output/summary.rds")

indexed_dfs <- vector(mode = "list", length = length(summary_deaths))

#Calculating heap indexes
for (i in c(47,53,83)){
  #Checking if the df didn't read in, we skip it.
  if (is.null(summary_deaths[[i]])){
    next
  }
  
  df <- summary_deaths[[i]]%>%
    mutate(hillChoi1014 = int_12_13*5/(int_10_11 + int_11_12 + int_12_13 + int_13_14 + int_14_15),
           pullum1014 = (int_12_13 - ((int_10_11 + int_11_12 + int_12_13 + int_13_14 + int_14_15)/5))/
             ((int_10_11 + int_11_12 + int_12_13 + int_13_14 + int_14_15)/5)*100)
  indexed_dfs[[i]] <- df
}

#Setting up folders for plots and tables
#Also making a container for folders to make it easier

plot_folders <- rep("", length(indexed_dfs))
for (i in seq_along(indexed_dfs)){
  #dir.create(path = str_c("../plots/", str_replace_all(countries[i], pattern = "\\s", replacement = "-"), "_", years_df[i]))
  plot_folders[i] <- str_c("../plots/", str_replace_all(countries[i], pattern = "\\s", replacement = "-"), "_", unique(indexed_dfs[[i]]$survey_year))
}

#Making a vector of ages that is numeric, to help make plots:
numeric_ages <- rep((1:24), 5)

i<-1
j<-1
#Making plots for every period within every survey

