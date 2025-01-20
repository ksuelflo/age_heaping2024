source("sim_functions.R")

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_range1_periodLength_periods_proportion1_seed

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
# "sim_ROW=1_numSIM=1_lnormmean=15_lnormsd=2_samplesize=1000_age1=12_range1=921_periodLength=12_periods=5_proportion1=01_seed=1"

# Catch the "task ID"
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#Getting all generated data frames' file paths.
all_folders <- list.files(path = "../data")
data_folder <- str_subset(all_folders, pattern = str_c("ROW=", task_id, "_"))
model_folder <- str_c("../models2/", str_replace(data_folder, "sim", "res"))
summary_folder <- str_c("../summaries2/", str_replace(data_folder, "sim", "summary"))

model_folder <- str_remove_all(model_folder, "\\.(?=\\d)")
summary_folder <- str_remove_all(summary_folder, "\\.(?=\\d)")

#NOTE: data folder, model folder, and summary folder will have 2 folders in them: this is due to an error I made in the file paths when
#generating the data. They are both the exact same sim setting, just with different sample sizes. So keep this in mind!

#CONT. If you only wanted to run the sample sizes 300, and 1000, add this line of code before `model_folder` is created:
#data_folder <- str_subset(data_folder, pattern = "samplesize=1000|samplesize=300")


#You'll have to create directories for the models and summaries if you haven't already.

#IMPORTANT: this data line is used to help fit logquad model.
data("fin1933")
	
#Getting all the different dfs for this task_id's sim setting
files <- list.files(path = str_c("../data/", data_folder))

#The below 2 lines of code are for if not all the models have been run during a job run on the cluster. The for loop will begin wherever
#The next data frame to be modeled is. 
summaries <- list.files(path = summary_folder)
start_int <- length(summaries) + 1 

#For this task_id, we loop through each generated df (looping through all 1000 different seeds). 
for (j in start_int:1000){
  sim <- readRDS(str_c("../data/", data_folder, "/", files[j]))
  
  #fit surv_synthetic() with lognormal distribution WITHOUT adjusting for age heaping.
  res_surv_lnorm <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
  #Table of summary results. We will rbind to this data frame throughout the rest of the models.
  summary_results <- get_uncertainty_surv_synthetic(res_surv_lnorm)%>%
    mutate(model = "surv_synth_lnorm", .before = period)
  #write(str_c("../data/", all_folders[i], "/", files[j]), file = "../output/do_sim_output.txt")
  #------------------------------------------------
  
  #fit surv_synthetic() with lognormal distribution ADJUSTING for age heaping (interval censored at 9-21 months)
  
  add_intervals <- sim%>%
    mutate(left_interval = if_else(right_interval > 9 & right_interval <= 21, 9, left_interval),
           right_interval = if_else(right_interval > 9 & right_interval <= 21, 21, right_interval),
           across_boundary = if_else(age_at_begin > left_interval & age_at_begin < right_interval, 1, 0))%>%
    group_by(id)%>%
    mutate(across_boundary = if_else(sum(across_boundary) > 0, 1, 0))
  
  res_surv_lnorm_adjusted <- surv_synthetic(df = add_intervals, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                            I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                            t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
  
  summary_results <- rbind(summary_results, get_uncertainty_surv_synthetic(res_surv_lnorm_adjusted)%>%
                            mutate(model = "surv_synth_lnorm_adj", .before = period))
  #------------------------------------------------
  
  #fit surv_synthetic() with Weibull distribution WITHOUT age heaping adjustment
   res_surv_weibull <-  tryCatch(expr = {
  surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                     I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                     t_1i = "right_interval", numerical_grad = TRUE, dist = "weibull")
    summary_results <- rbind(summary_results, get_uncertainty_surv_synthetic(res_surv_weibull, distribution = "weibull")%>%
                           mutate(model = "surv_synth_weibull", .before = period))

},
error = function(e){
  message("got an error Weibull")
  return(NULL)
})

print("Made it past Weibull")
  #------------------------------------------------
  
  # reformat data to work with discrete hazards model - 6 age periods NO ADJUSTMENT for age heaping
  
  dh_data <- reformat_sims_disc_haz(sims = sim,
                                    agegroup_splits = c(0,1,12,24,36,48,60))
  
  res_discrete_hazard <- fit_disc_haz(dh_data)

  summary_results <- rbind(summary_results, summarize_disc_haz(res_discrete_hazard)%>%
                             mutate(model = "discrete_hazards", .before = period))
  print("made it past discrete hazard")
  #------------------------------------------------
  
  #SETUP for Logquad
  
  # get 22 age inputs for logquad model from monthly discrete hazards fit 
  dh_data_monthly <- reformat_sims_disc_haz(sims = sim,
                                            agegroup_splits = 0:60)
  
  res_discrete_hazard_monthly <- fit_disc_haz(dh_data_monthly)
  
  log_quad_pars <- get_log_quad_params(res_discrete_hazard_monthly)
  
  #log quad fitting
  
  #Where lower age always is 0, upper age is xq0, where x is upper age. n = upperage[i] - upperage[i-1]. rate is nqx. type is always qx
  #sex is always total. Fit is always min, besides last observation (5q0), which is match.
  
  #Calculating weights
  weight <- c(fin1933$n[1:22]/(365.25*5), NA)
  
  log_quad_pars <- log_quad_pars%>%
    mutate(lower_age = 0,
           type = "qx", 
           sex = "total",
           fit = rep(fin1933$fit, 5),
           weight = rep(weight, 5))
  
  #Need to have all 7 columns of `fin1933`
  
  #--------------------------------------------------------
  
  #Logquad WITHOUT adjusting for age heaping
  
  res_log_quad <- vector("list", length = 5)
  
  #Setting up summary table for logquad
  NMR <- rep(0,5)
  IMR <- rep(0,5)
  U5MR <- rep(0,5)
  NMR_lower <- rep(0,5)
  NMR_upper <- rep(0,5)
  IMR_lower <- rep(0,5)
  IMR_upper <- rep(0,5)
  U5MR_lower <- rep(0,5)
  U5MR_upper <- rep(0,5)
  period <- 1:5
  model <- rep("logquad", 5)
 
  #Looping through each period
  for (k in 1:5){
    input <- format_data(
      lower_age = log_quad_pars$lower_age[log_quad_pars$period == k],
      upper_age = log_quad_pars$upper_age[log_quad_pars$period == k],
      rate      = log_quad_pars$rate[log_quad_pars$period == k],
      type      = log_quad_pars$type[log_quad_pars$period == k],
      sex       = log_quad_pars$sex[log_quad_pars$period == k],
      fit       = log_quad_pars$fit[log_quad_pars$period == k],
      weight    = log_quad_pars$weight[log_quad_pars$period == k])
    
    tryCatch(expr= {
      res_log_quad[[k]] <- lagrange5q0(data = input)
      NMR[k] <- res_log_quad[[k]]$predictions[5,4]
      IMR[k] <- res_log_quad[[k]]$predictions[16,4]
      U5MR[k] <- res_log_quad[[k]]$predictions[22,4]
      NMR_lower[k] <- res_log_quad[[k]]$predictions[5,5]
      NMR_upper[k] <- res_log_quad[[k]]$predictions[5,6]
      IMR_lower[k] <- res_log_quad[[k]]$predictions[16,5]
      IMR_upper[k] <- res_log_quad[[k]]$predictions[16,6]
      U5MR_lower[k] <- res_log_quad[[k]]$predictions[22,5]
      U5MR_upper[k] <- res_log_quad[[k]]$predictions[22,6]
    },
    error = function(e){
      NMR[k] <- NA
      IMR[k] <- NA
      U5MR[k] <- NA
      NMR_lower[k] <- NA
      NMR_upper[k] <- NA
      IMR_lower[k] <- NA
      IMR_upper[k] <- NA
      U5MR_lower[k] <- NA
      U5MR_upper[k] <- NA
    }
    )
    
  }
  print("made it past first log quad try catch")
  
  summary_lquad <- data.frame(model = model, period = period, NMR = NMR, IMR = IMR, U5MR = U5MR,
                              NMR_lower = NMR_lower, NMR_upper = NMR_upper, IMR_lower = IMR_lower, 
                              IMR_upper = IMR_upper, U5MR_lower = U5MR_lower, U5MR_upper = U5MR_upper)
  summary_results <- rbind(summary_results, summary_lquad)
  names(res_log_quad) <- c("period_1", "period_2", "period_3", "period_4", "period_5")
  
  #--------------------------------------------------------
  
  #Logquad ADJUSTING for age heaping (removing x/qx which is in the interval 9,21 months)
  
  #Filtering out x/qxs
  
sum_adj <- sum(log_quad_pars$weight[12:18])


  log_quad_pars_adj <- log_quad_pars%>%
    filter(upper_age < (9*30.4375) | upper_age > (21*30.4375))
log_quad_pars_adj$weight[11] <- log_quad_pars_adj$weight[11] + sum_adj   
  res_log_quad_adj <- vector("list", length = 5)
  #Setting up summary table for logquad
  NMR_adj <- rep(0,5)
  IMR_adj <- rep(0,5)
  U5MR_adj <- rep(0,5)
  # NMR_lower_adj <- rep(0,5)
  # NMR_upper_adj <- rep(0,5)
  # IMR_lower_adj <- rep(0,5)
  # IMR_upper_adj <- rep(0,5)
  # U5MR_lower_adj <- rep(0,5)
  # U5MR_upper_adj <- rep(0,5)
  NMR_lower_adj <- rep(NA,5)
  NMR_upper_adj <- rep(NA,5)
  IMR_lower_adj <- rep(NA,5)
  IMR_upper_adj <- rep(NA,5)
  U5MR_lower_adj <- rep(NA,5)
  U5MR_upper_adj <- rep(NA,5)
  period_adj <- 1:5
  model_adj <- rep("logquad", 5)


  
  #Looping through each period
  for (w in 1:5){
    input_adj <- format_data_mine(
      lower_age = log_quad_pars_adj$lower_age[log_quad_pars_adj$period == w],
      upper_age = log_quad_pars_adj$upper_age[log_quad_pars_adj$period == w],
      rate      = log_quad_pars_adj$rate[log_quad_pars_adj$period == w],
      type      = log_quad_pars_adj$type[log_quad_pars_adj$period == w],
      sex       = log_quad_pars_adj$sex[log_quad_pars_adj$period == w],
      fit       = log_quad_pars_adj$fit[log_quad_pars_adj$period == w],
      weight    = log_quad_pars_adj$weight[log_quad_pars_adj$period == w])
    
    tryCatch({
      res_log_quad_adj[[w]] <- lagrange5q0(data = input_adj)
      NMR_adj[w] <- res_log_quad_adj[[w]]$predictions[5,4]
      IMR_adj[w] <- res_log_quad_adj[[w]]$predictions[16,4]
      U5MR_adj[w] <- res_log_quad_adj[[w]]$predictions[22,4]
      # NMR_lower_adj[j] <- res_log_quad_adj[[j]]$predictions[5,5]
      # NMR_upper_adj[j] <- res_log_quad_adj[[j]]$predictions[5,6]
      # IMR_lower_adj[j] <- res_log_quad_adj[[j]]$predictions[16,5]
      # IMR_upper_adj[j] <- res_log_quad_adj[[j]]$predictions[16,6]
      # U5MR_lower_adj[j] <- res_log_quad_adj[[j]]$predictions[22,5]
      # U5MR_upper_adj[j] <- res_log_quad_adj[[j]]$predictions[22,6]
    },
    error = function(e){
      NMR_adj[w] <- NA
      IMR_adj[w] <- NA
      U5MR_adj[w] <- NA
      # NMR_lower_adj[j] <- NA
      # NMR_upper_adj[j] <- NA
      # IMR_lower_adj[j] <- NA
      # IMR_upper_adj[j] <- NA
      # U5MR_lower_adj[j] <- NA
      # U5MR_upper_adj[j] <- NA
    }
    )
    
  } 
  summary_lquad_adj <- data.frame(model = model_adj, period = period_adj, NMR = NMR_adj, IMR = IMR_adj, U5MR = U5MR_adj,
                                  NMR_lower = NMR_lower_adj, NMR_upper = NMR_upper_adj, IMR_lower = IMR_lower_adj, 
                                  IMR_upper = IMR_upper_adj, U5MR_lower = U5MR_lower_adj, U5MR_upper = U5MR_upper_adj)
  
  summary_results <- rbind(summary_results, summary_lquad_adj)
  
  names(res_log_quad_adj) <- c("period_1", "period_2", "period_3", "period_4", "period_5")
  
  #------------------------------------------------
  
  #Saving all model objects/summaries into a list
  list_res <- list(surv_lnorm = res_surv_lnorm, surv_lnorm_adjusted = res_surv_lnorm_adjusted, 
                   surv_weibull = res_surv_weibull, discrete_hazard = res_discrete_hazard, log_quad = res_log_quad,
                   log_quad_adjusted = res_log_quad_adj)
  #Save rds object to correct directory
  saveRDS(list_res, str_c(model_folder, "/", str_replace(files[j], "sim", "res")))

  saveRDS(summary_results, str_c(summary_folder, "/", str_replace(files[j], "sim", "summary")))
}

#Run the for loop above.

#Then use warnings() to look through the warnings output from logquad. 

