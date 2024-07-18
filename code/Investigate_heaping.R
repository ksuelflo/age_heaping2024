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

folders <- list.files(path = "../../data")
years <- as.numeric(str_extract(folders, "\\d+"))
dfs <- vector(mode = "list", length = length(folders))
countries <- str_extract(folders, ".*(?=_(\\d)*)")

#This for loop processes data and adds it to big container `dfs`
for (i in seq_along(folders)){
  print(i)
  #Getting file path
  path <- str_subset(string = list.files(path = str_c("../../data/", folders[i])), pattern = ".DTA") 
  
  #Reading in data
  births <- read.dta13(str_c("../../data/", folders[i], "/", path), generate.factors = TRUE)
  
  #Survey year
  year <- years[i]
  
  #formats raw birth recode file.
  dat <- format_dhs(df = births, survey_year = year, period_boundaries = c(year-5, year-4, year-3, year-2, year-1, year))
  dfs[[i]] <- dat
}

#Setting up for loop
for (i in seq_along(dfs)){

}

#Going to practice outside for loop

df <- dfs[[1]]

summary_heap_df <- df%>%
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
  mutate(country = countries[1], #Use `i` instead when in the for loop.
         survey_year = years[1], 
         .before = p)%>% #Use `i` instead when in the for loop.
  mutate(p_begin = survey_year - 6 + p,
         p_end = p_begin + 1,
         .after = p)
  

View(summary_heap_df)
summary_heap_df

x <- 1:10
cut(x, breaks = 2)