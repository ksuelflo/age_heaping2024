source("sim_functions.R")

# Catch the “task ID”
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(str_c("The task ID is: ", task_id))

#Ages at heap, ranges, and proportions are all split up into two columns. This makes things easier. 

#Simulation settings as found in overleaf.

#For lnorm: We only want (12, 1.9), (20, 2.1), (15, 2.1), (20, 2.3)

lnorm_mean_vec <- c(12,15,20)
lnorm_sd_vec <- c(1.9,2.1,2.3)
sample_size_vec <- c(300,1000)
age_1_vec <- 12
range_1_vec <- c("9,21", "6,18", "8,24")
period_length_vec <- c(60)
periods_vec <- 5
proportion_1_vec <- c(0,.10,.20,.50)

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

#Setting up folder name for the given task_id
foldername <- str_c("sim_ROW=", task_id, "_lnormmean=", clean_params_3$lnorm_mean[task_id], "_lnormsd=", clean_params_3$lnorm_sd[task_id],
              "_samplesize=", clean_params_3$sample_size[task_id], "_age=", clean_params_3$age_1[task_id], "_range=", clean_params_3$range_1[task_id],
              "_periodlength=", clean_params_3$period_length[task_id], "_periods=", clean_params_3$periods[task_id],
            "_proportion=", clean_params_3$proportion_1[task_id])

#Create the directory if it is not already: If it is, the code won't break, it just won't recreate the directory.
dir.create(path = str_c("../data/", foldername))


#Simulating 1000 times (getting 1000 data frames) for the given setting (based on task_id)
for (j in 1:1000){

  #Setting seed for reproducibility
  seed <- 1000*task_id + j
  set.seed(seed)
  
  #Formatting clean_params for general_sim()
  mus_lnorm <- rep(clean_params_3$lnorm_mean[task_id], 5)
  sds_lnorm <- exp(rep(clean_params_3$lnorm_sd[task_id], 5))
  period_lengths <- rep(clean_params_3$period_length[task_id], 5)
  matrix_params <- cbind(mus_lnorm, sds_lnorm, period_lengths)
  
  #Simulated data BEFORE age heaping
  before_heaping <- general_sim(num_child = clean_params_3$sample_size[task_id], param_matrix = matrix_params, distribution = "lognormal")
  
  #Formatting data for sim_age_heap()
  ages <- c(clean_params_3$age_1[task_id], clean_params_3$age_2[task_id])
  proportions <- c(clean_params_3$proportion_1[task_id], clean_params_3$proportion_2[task_id])
  bottom_range <- c(as.numeric(str_extract(clean_params_3$range_1[task_id], "\\d*(?=,)")), as.numeric(str_extract(clean_params_3$range_2[task_id], "\\d*(?=,)")))
  top_range <- c(as.numeric(str_extract(clean_params_3$range_1[task_id], "(?<=,)\\d*")), as.numeric(str_extract(clean_params_3$range_2[task_id], "(?<=,)\\d*")))
  ranges <- cbind(bottom_range, top_range)

  #Simulated data AFTER age heaping
  sim_with_heap <- sim_age_heap(ages_at_heaping = ages, proportion_heap = proportions, range_heap = ranges, sim_data = before_heaping)

  #Clean data (interval censor is main bit).
  cleaned_data <- clean_data(sim_with_heap)
  
  #Build out file name.
  filename <- str_c("sim_ROW=", task_id+96, "_numSIM=", j, "_lnormmean=", clean_params_3$lnorm_mean[task_id], "_lnormsd=", clean_params_3$lnorm_sd[task_id],"_samplesize=", clean_params_3$sample_size[task_id], "_age=", clean_params_3$age_1[task_id], "_range=", clean_params_3$range_1[task_id],"_periodlength=", clean_params_3$period_length[task_id], "_periods=", clean_params_3$periods[task_id], "_proportion=", clean_params_3$proportion_1[task_id], "_seed=", seed)   
  #Removing both periods and commas. Based on context of numbers, we don't need them, as they mess up files.
  filename <- str_replace_all(filename, ",", "")
  filename <- str_replace_all(filename, "\\.", "")

  #Save the generated data frame to the correct directory.
  saveRDS(cleaned_data, file = str_c("../data/", folder, "/", filename, ".rds"))
}

#NAMING CONVENTION

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_range1_periodLength_periods_proportion1_seed

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
# "sim_ROW=1_numSIM=1_lnormmean=15_lnormsd=2_samplesize=1000_age1=12_range1=921_periodLength=12_periods=5_proportion1=01_seed=1"
