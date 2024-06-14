library(tidyverse)
source("sim_functions.R")

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



view(clean_params)

# sim_param_param_.....

#Looping through each simulation setting
for (i in 1:nrow(clean_params)){
  
  #Simulating 1000 times (getting 1000 data frames) for each setting
  for (j in 1: 1000){
  
    #Formatting clean_params for general_sim()
    mus_lnorm <- rep(clean_params$lnorm_mean[i], 5)
    sds_lnorm <- rep(clean_params$lnorm_sd[i], 5)
    period_lengths <- rep(clean_params$period_length, 5)
    matrix_params <- cbind(mus_lnorm, sds_lnorm, period_lengths)
    
    #Simulated data BEFORE age heaping
    before_heaping <- general_sim(num_child = clean_params$sample_size[i], param_matrix = matrix_params, distribution = "lognormal")
  
    #Formatting data for sim_age_heap()
    ages <- c(clean_params$age_1, clean_params$age_2)
    proportions <- c(clean_params$proportion_1, clean_params$proportion_2)
    bottom_range <- c(as.numeric(str_extract(clean_params$range_1[i], "\\d*(?=,)")), as.numeric(str_extract(clean_params$range_2[i], "\\d*(?=,)")))
    top_range <- c(as.numeric(str_extract(clean_params$range_1[i], "(?<=,)\\d*")), as.numeric(str_extract(clean_params$range_2[i], "(?<=,)\\d*")))
    ranges <- cbind(bottom_range, top_range)
    
    #Simulated data AFTER age heaping
    sim_with_heap <- sim_age_heap(ages_at_heaping = ages, proportion_heap = proportions, range_heap = ranges, sim_data = before_heaping)
    
    write_rds(sim_with_heap, str_c("sim_", i, "_", j, "_", clean_params$lnorm_mean, "_", clean_params$lnorm_sd, "_", clean_params$sample_size,
                                   "_", clean_params$age_1, "_", clean_params$age_2, "_", clean_params$range_1))
    
  }
}

#NAMING CONVENTION

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_age2_range1_range2_periodLength_periods_proportion1_proportion2

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
#           "sim_1_1_15_2_1000_12_NA_9,21_NA_12_5_.01_NA.rds"


strdgs <- str_c("hi", as.character(NA))




