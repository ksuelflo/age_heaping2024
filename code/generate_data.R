library(tidyverse)

#Ages at heap, ranges, and proportions are all split up into two columns. This makes things easier. 

lnorm_mean_vec <- c(15,20)
lnorm_sd_vec <- c(2,2.2)
sample_size_vec <- c(1000,5000)
age_1_vec <- 12
age_2_vec <- c(NA, 6,60)
range_1_vec <- c("9,21", "6,18", "8,24")
range_2_vec <- c("3,9", "21,27", "57, 63")
period_length_vec <- c(12)
periods_vec <- 5
proportion_1_vec <- c(.01,.05,.1)
proportion_2_vec <- c(.01,.05,.1)



raw_expand <- expand.grid(lnorm_mean = lnorm_mean_vec, lnorm_sd = lnorm_sd_vec, sample_size = sample_size_vec, age_1 = age_1_vec, 
                          age_2 = age_2_vec, range_1 = range_1_vec,range_2 = range_2_vec, period_length = period_length_vec, 
                          periods = periods_vec, proportion_1 = proportion_1_vec, proportion_2 = proportion_2_vec, stringsAsFactors = FALSE)

clean_params <- raw_expand%>%
  filter(age_2 == (as.numeric(str_extract(range_2, "\\d*(?=,)")) + 3) | is.na(age_2))%>%
  mutate(range_2 = if_else(is.na(age_2), NA, range_2),
         proportion_2 = if_else(is.na(age_2), NA, proportion_2))



nrow(clean_params)

# sim_param_param_.....











