source("sim_functions.R")
library(tidyr)

#Ages at heap, ranges, and proportions are all split up into two columns. This makes things easier. 

#Simulation settings as found in overleaf.

#For lnorm: We only want (12, 1.9), (20,2.1), (15, 2.1), 20(2.3)

lnorm_mean_vec <- c(12,15,20)
lnorm_sd_vec <- c(1.9,2.1,2.3)
sample_size_vec <- c(100,500)
age_1_vec <- 12
# age_2_vec <- c(NA, 6,60)
range_1_vec <- c("9,21", "6,18", "8,24")
# range_2_vec <- c(NA, "3,9", "57, 60")
period_length_vec <- c(60)
periods_vec <- 5
proportion_1_vec <- c(0,.10,.20,.50)
# proportion_2_vec <- c(.05,.10)


#Every possible combination of sim settings found using expand.grid. Even ones that aren't plausible.

clean_params <- expand.grid(lnorm_mean = lnorm_mean_vec, lnorm_sd = lnorm_sd_vec, sample_size = sample_size_vec, age_1 = age_1_vec, 
                          range_1 = range_1_vec, period_length = period_length_vec, 
                          periods = periods_vec, proportion_1 = proportion_1_vec, stringsAsFactors = FALSE)

#Filtering out combos to only include plausible combinations.
clean_params_2 <- clean_params%>%
  mutate(double = lnorm_sd + lnorm_mean)

#Filtering out combos to only include plausible combinations.
clean_params_3 <- clean_params_2%>%
  filter(double %in% c(13.9, 17.3, 17.1, 22.3))%>%
  dplyr::select(-double)
view(clean_params_3)

#Code for second age heaping.

# clean_params <- raw_expand%>%
#   filter(age_2 == (as.numeric(str_extract(range_2, "\\d*(?=,)")) + 3) | is.na(age_2))%>%
#   mutate(range_2 = if_else(is.na(age_2), NA, range_2),
#          proportion_2 = if_else(is.na(age_2), NA, proportion_2)
#          )



#Looping through each simulation setting
for (i in 1:nrow(clean_params_3)){
  
  if (i != 54){
    next
  }
  
  start <- Sys.time()
  
  #Simulating 1000 times (getting 1000 data frames) for each setting
  for (j in 1: 1){
    
    #Setting seed for reproducibility
    seed <- 1000*i + j
    set.seed(seed)
    
    #Formatting clean_params for general_sim()
    mus_lnorm <- rep(clean_params_3$lnorm_mean[i], 5)
    sds_lnorm <- exp(rep(clean_params_3$lnorm_sd[i], 5))
    period_lengths <- rep(clean_params_3$period_length[i], 5)
    matrix_params <- cbind(mus_lnorm, sds_lnorm, period_lengths)
    
    #Simulated data BEFORE age heaping
    before_heaping <- general_sim(num_child = clean_params_3$sample_size[i], param_matrix = matrix_params, distribution = "lognormal")
    
    print("Before heaping went well.")
    #Formatting data for sim_age_heap()
    ages <- c(clean_params_3$age_1[i], clean_params_3$age_2[i])
    proportions <- c(clean_params_3$proportion_1[i], clean_params_3$proportion_2[i])
    bottom_range <- c(as.numeric(str_extract(clean_params_3$range_1[i], "\\d*(?=,)")), as.numeric(str_extract(clean_params_3$range_2[i], "\\d*(?=,)")))
    top_range <- c(as.numeric(str_extract(clean_params_3$range_1[i], "(?<=,)\\d*")), as.numeric(str_extract(clean_params_3$range_2[i], "(?<=,)\\d*")))
    ranges <- cbind(bottom_range, top_range)
  
    #Simulated data AFTER age heaping
    sim_with_heap <- sim_age_heap(ages_at_heaping = ages, proportion_heap = proportions, range_heap = ranges, sim_data = before_heaping)

    
    print("Heaping went well.")
    #Clean data (interval censor is main bit).
    cleaned_data <- clean_data(sim_with_heap)
    
    # Checking for NA values, to convert them to character "NA". THIS IS FOR FILE NAMING, THAT IS ALL. THIS IS FOR MULTIPLE AGES
    # if (is.na(clean_params$age_2[i])){
    #   age_2 <- "NA"
    #   proportion_2 <- "NA"
    #   range_2 <- "NA"
    # }
    # #If not NA, just use normal.
    # else{
    #   age_2 <- clean_params$age_2[i]
    #   proportion_2 <- clean_params$proportion_2[i]
    #   range_2 <- clean_params$range_2[i]
    # }
    
    #Build out file name.
    filename <- str_c("sim_ROW=", i, "_numSIM=", j, "_lnormmean=", clean_params_3$lnorm_mean[i], "_lnormsd=", clean_params_3$lnorm_sd[i],
                      "_samplesize=", clean_params_3$sample_size[i], "_age=", clean_params_3$age_1[i], "_range=", clean_params_3$range_1[i], 
                      "_periodlength=", clean_params_3$period_length[i], "_periods=", clean_params_3$periods[i], 
                      "_proportion=", clean_params_3$proportion_1[i], "_seed=", seed)
    #Removing both periods and commas. Based on context of numbers, we don't need them, as they mess up files.
    filename <- str_replace_all(filename, ",", "")
    filename <- str_replace_all(filename, "\\.", "")
    print(filename)
    saveRDS(cleaned_data, file = str_c("../data/", filename, ".rds"))
  }
  
  end <- Sys.time()
  print(start-end)
}


# res_surv_weibull <- surv_synthetic(df = cleaned_data, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
#                                    I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
#                                    t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
# 
# cleaned_data_2 <- cleaned_data%>%
#   mutate(indicator = if_else(period == 5 & binom_heap == 1 & max_age < left_interval, 1, 0))
# 
# cleaned_data_2%>%
#   filter(indicator == 1)



#NAMING CONVENTION

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_range1_periodLength_periods_proportion1_seed

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
# "sim_ROW=1_numSIM=1_lnormmean=15_lnormsd=2_samplesize=1000_age1=12_range1=921_periodLength=12_periods=5_proportion1=01_seed=1"


