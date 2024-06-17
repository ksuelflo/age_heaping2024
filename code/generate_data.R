source("sim_functions.R")

#Ages at heap, ranges, and proportions are all split up into two columns. This makes things easier. 

#Simulation settings as found in overleaf.

lnorm_mean_vec <- c(15,20)
lnorm_sd_vec <- c(2,2.2)
sample_size_vec <- c(1000,5000)
age_1_vec <- 12
age_2_vec <- c(NA, 6,60)
range_1_vec <- c("9,21", "6,18", "8,24")
range_2_vec <- c(NA, "3,9", "57, 60")
period_length_vec <- c(12)
periods_vec <- 5
proportion_1_vec <- c(.01,.05,.1)
proportion_2_vec <- c(.01,.05,.1)


#Every possible combination of sim settings found using expand.grid. Even ones that aren't plausible.

raw_expand <- expand.grid(lnorm_mean = lnorm_mean_vec, lnorm_sd = lnorm_sd_vec, sample_size = sample_size_vec, age_1 = age_1_vec, 
                          age_2 = age_2_vec, range_1 = range_1_vec,range_2 = range_2_vec, period_length = period_length_vec, 
                          periods = periods_vec, proportion_1 = proportion_1_vec, proportion_2 = proportion_2_vec, stringsAsFactors = FALSE)

#Filtering out combos to only include plausible combinations.

clean_params <- raw_expand%>%
  filter(age_2 == (as.numeric(str_extract(range_2, "\\d*(?=,)")) + 3) | is.na(age_2))%>%
  mutate(range_2 = if_else(is.na(age_2), NA, range_2),
         proportion_2 = if_else(is.na(age_2), NA, proportion_2)
         )

#Looping through each simulation setting
for (i in 1:nrow(clean_params)){
  #Simulating 1000 times (getting 1000 data frames) for each setting
  for (j in 1: 1000){
    
    #Setting seed for reproducibility
    seed <- 1000*i + j
    set.seed(seed)
    
    #Formatting clean_params for general_sim()
    mus_lnorm <- rep(clean_params$lnorm_mean[i], 5)
    sds_lnorm <- exp(rep(clean_params$lnorm_sd[i], 5))
    period_lengths <- rep(clean_params$period_length[i], 5)
    matrix_params <- cbind(mus_lnorm, sds_lnorm, period_lengths)
    
    #Simulated data BEFORE age heaping
    before_heaping <- general_sim(num_child = clean_params$sample_size[i], param_matrix = matrix_params, distribution = "lognormal")
    
    print("Before heaping went well.")
    #Formatting data for sim_age_heap()
    ages <- c(clean_params$age_1[i], clean_params$age_2[i])
    proportions <- c(clean_params$proportion_1[i], clean_params$proportion_2[i])
    bottom_range <- c(as.numeric(str_extract(clean_params$range_1[i], "\\d*(?=,)")), as.numeric(str_extract(clean_params$range_2[i], "\\d*(?=,)")))
    top_range <- c(as.numeric(str_extract(clean_params$range_1[i], "(?<=,)\\d*")), as.numeric(str_extract(clean_params$range_2[i], "(?<=,)\\d*")))
    ranges <- cbind(bottom_range, top_range)
  
    #Simulated data AFTER age heaping
    sim_with_heap <- sim_age_heap(ages_at_heaping = ages, proportion_heap = proportions, range_heap = ranges, sim_data = before_heaping)

    
    print("Heaping went well.")
    #Clean data (interval censor is main bit.
    cleaned_data <- clean_data(sim_with_heap)
    
    #Checking for NA values, to convert them to character "NA". THIS IS FOR FILE NAMING, THAT IS ALL.
    if (is.na(clean_params$age_2[i])){
      age_2 <- "NA"
      proportion_2 <- "NA"
      range_2 <- "NA"
    }
    #If not NA, just use normal.
    else{
      age_2 <- clean_params$age_2[i]
      proportion_2 <- clean_params$proportion_2[i]
      range_2 <- clean_params$range_2[i]
    }
    
    #Build out file name.
    filename <- str_c("sim_", i, "_", j, "_", clean_params$lnorm_mean[i], "_", clean_params$lnorm_sd[i], "_", clean_params$sample_size[i],
                      "_", clean_params$age_1[i], "_", age_2, "_", clean_params$range_1[i], "_", range_2, "_", clean_params$period_length[i],
                      "_", clean_params$periods[i], "_", clean_params$proportion_1[i], "_", proportion_2, "_", seed)
    #Removing both periods and commas. Based on context of numbers, we don't need them, as they mess up files.
    filename <- str_replace_all(filename, ",", "")
    filename <- str_replace_all(filename, "\\.", "")
    print(filename)
    saveRDS(cleaned_data, file = str_c("../data/", filename))
    
    # break
  }
  # break
}

#NAMING CONVENTION

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_age2_range1_range2_periodLength_periods_proportion1_proportion2_seed

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
#           "sim_1_1_15_2_1000_12_NA_9,21_NA_12_5_.01_NA_1.rds"


