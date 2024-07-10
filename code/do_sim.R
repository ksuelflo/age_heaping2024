source("sim_functions.R")
library(tibble)
library(tidyr)

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_range1_periodLength_periods_proportion1_seed

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
# "sim_ROW=1_numSIM=1_lnormmean=15_lnormsd=2_samplesize=1000_age1=12_range1=921_periodLength=12_periods=5_proportion1=01_seed=1"


#Getting all generated data frames' file paths.
all_files <- list.files(path = "/Users/kylesuelflow/Macalester-Stuff/Research-Taylor/age_heaping2024/data")

#IMPORTANT: this data line is used to help fit logquad model.
data("fin1933")
log_quad_fail_count <- 0
log_quad_fail_count_adj <- 0

inputs_row_51 <- vector("list", length = 5)
res_row_51 <- vector("list", length = 5)

for (i in seq_along(all_files)){
  
  # print(all_files[i])
  
  if (i != 47){
    next
  }
  # else{
  #   print(all_files[i])
  #   break
  # }
  # 
  sim <- readRDS(str_c("../data/", all_files[i]))
  
  #fit surv_synthetic() with lognormal distribution WITHOUT adjusting for age heaping.
  
  res_surv_lnorm <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")

  #Table of summary results. We will rbind to this data frame throughout the rest of the models.
  summary_results <- get_uncertainty_surv_synthetic(res_surv_lnorm)%>%
    mutate(model = "surv_synth_lnorm", .before = period)

  print("fit lnorm model")
  #------------------------------------------------
  
  #fit surv_synthetic() with lognormal distribution ADJUSTING for age heaping (interval censored at 9-21 months)
  
  
  
  #A_i needs to be updated here.
  
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
  
  print("fit lnorm adjusted model")
  #------------------------------------------------
  
  #fit surv_synthetic() with Weibull distribution WITHOUT age heaping adjustment
  
  res_surv_weibull <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "weibull")
  
  summary_results <- rbind(summary_results, get_uncertainty_surv_synthetic(res_surv_weibull, distribution = "weibull")%>%
                             mutate(model = "surv_synth_weibull", .before = period))
  
  print("fit weibull model")
  #------------------------------------------------
  
  # reformat data to work with discrete hazards model - 6 age periods NO ADJUSTMENT for age heaping
  
  dh_data <- reformat_sims_disc_haz(sims = sim,
                         agegroup_splits = c(0,1,12,24,36,48,60))
  
  res_discrete_hazard <- fit_disc_haz(dh_data)
  
  summary_results <- rbind(summary_results, summarize_disc_haz(res_discrete_hazard)%>%
                             mutate(model = "discrete_hazards", .before = period))

  print("fit discrete hazards")
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
  for (j in 1:5){
    input <- format_data(
      lower_age = log_quad_pars$lower_age[log_quad_pars$period == j],
      upper_age = log_quad_pars$upper_age[log_quad_pars$period == j],
      rate      = log_quad_pars$rate[log_quad_pars$period == j],
      type      = log_quad_pars$type[log_quad_pars$period == j],
      sex       = log_quad_pars$sex[log_quad_pars$period == j],
      fit       = log_quad_pars$fit[log_quad_pars$period == j],
      weight    = log_quad_pars$weight[log_quad_pars$period == j])
    
    inputs_row_51[[j]] <- input
    
    tryCatch({
      res_log_quad[[j]] <- lagrange5q0(data = input)
      res_row_51[[j]] <- lagrange5q0(data = input)
      
      NMR[j] <- res_log_quad[[j]]$predictions[5,4]
      IMR[j] <- res_log_quad[[j]]$predictions[16,4]
      U5MR[j] <- res_log_quad[[j]]$predictions[22,4]
      NMR_lower[j] <- res_log_quad[[j]]$predictions[5,5]
      NMR_upper[j] <- res_log_quad[[j]]$predictions[5,6]
      IMR_lower[j] <- res_log_quad[[j]]$predictions[16,5]
      IMR_upper[j] <- res_log_quad[[j]]$predictions[16,6]
      U5MR_lower[j] <- res_log_quad[[j]]$predictions[22,5]
      U5MR_upper[j] <- res_log_quad[[j]]$predictions[22,6]
    },
      error = function(e){
        message("log quad failed. Model equals NA.")
        print(e)
        NMR[j] <- NA
        IMR[j] <- NA
        U5MR[j] <- NA
        NMR_lower[j] <- NA
        NMR_upper[j] <- NA
        IMR_lower[j] <- NA
        IMR_upper[j] <- NA
        U5MR_lower[j] <- NA
        U5MR_upper[j] <- NA
        log_quad_fail_count <- log_quad_fail_count + 1
      }
    )
    
  }
  
  summary_lquad <- data.frame(model = model, period = period, NMR = NMR, IMR = IMR, U5MR = U5MR,
                              NMR_lower = NMR_lower, NMR_upper = NMR_upper, IMR_lower = IMR_lower, 
                              IMR_upper = IMR_upper, U5MR_lower = U5MR_lower, U5MR_upper = U5MR_upper)
  
  summary_results <- rbind(summary_results, summary_lquad)
  
  names(res_log_quad) <- c("period_1", "period_2", "period_3", "period_4", "period_5")
  
  print("fit logquad no adjust")
  
  #--------------------------------------------------------
  
  #Logquad ADJUSTING for age heaping (removing x/qx which is in the interval 9,21 months)
  
  #Filtering out x/qxs
  
  log_quad_pars_adj <- log_quad_pars%>%
    filter(upper_age < (9*30.4375) | upper_age > (21*30.4375))%>%
    mutate(weight = rep(c(fin1933$n[setdiff(1:22, 12:18)]/(sum(fin1933$n[setdiff(1:22, 12:18)])), NA), times = 5))
  
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
  NMR_lower <- rep(NA,5)
  NMR_upper <- rep(NA,5)
  IMR_lower <- rep(NA,5)
  IMR_upper <- rep(NA,5)
  U5MR_lower <- rep(NA,5)
  U5MR_upper <- rep(NA,5)
  period_adj <- 1:5
  model_adj <- rep("logquad", 5)
  
  #Looping through each period
  for (j in 1:5){
    input_adj <- format_data(
      lower_age = log_quad_pars_adj$lower_age[log_quad_pars_adj$period == j],
      upper_age = log_quad_pars_adj$upper_age[log_quad_pars_adj$period == j],
      rate      = log_quad_pars_adj$rate[log_quad_pars_adj$period == j],
      type      = log_quad_pars_adj$type[log_quad_pars_adj$period == j],
      sex       = log_quad_pars_adj$sex[log_quad_pars_adj$period == j],
      fit       = log_quad_pars_adj$fit[log_quad_pars_adj$period == j],
      weight    = log_quad_pars_adj$weight[log_quad_pars_adj$period == j])
    
    tryCatch({
      res_log_quad_adj[[j]] <- lagrange5q0(data = input_adj)
      NMR_adj[j] <- res_log_quad_adj[[j]]$predictions[5,4]
      IMR_adj[j] <- res_log_quad_adj[[j]]$predictions[16,4]
      U5MR_adj[j] <- res_log_quad_adj[[j]]$predictions[22,4]
      # NMR_lower_adj[j] <- res_log_quad_adj[[j]]$predictions[5,5]
      # NMR_upper_adj[j] <- res_log_quad_adj[[j]]$predictions[5,6]
      # IMR_lower_adj[j] <- res_log_quad_adj[[j]]$predictions[16,5]
      # IMR_upper_adj[j] <- res_log_quad_adj[[j]]$predictions[16,6]
      # U5MR_lower_adj[j] <- res_log_quad_adj[[j]]$predictions[22,5]
      # U5MR_upper_adj[j] <- res_log_quad_adj[[j]]$predictions[22,6]
    },
    error = function(e){
      message("log quad failed. Model equals NA.")
      print(e)
      NMR_adj[j] <- NA
      IMR_adj[j] <- NA
      U5MR_adj[j] <- NA
      # NMR_lower_adj[j] <- NA
      # NMR_upper_adj[j] <- NA
      # IMR_lower_adj[j] <- NA
      # IMR_upper_adj[j] <- NA
      # U5MR_lower_adj[j] <- NA
      # U5MR_upper_adj[j] <- NA
      log_quad_fail_count_adj <- log_quad_fail_count_adj + 1
    }
    )
    
  }
  
  summary_lquad_adj <- data.frame(model = model_adj, period = period_adj, NMR = NMR_adj, IMR = IMR_adj, U5MR = U5MR_adj,
                              NMR_lower = NMR_lower_adj, NMR_upper = NMR_upper_adj, IMR_lower = IMR_lower_adj, 
                              IMR_upper = IMR_upper_adj, U5MR_lower = U5MR_lower_adj, U5MR_upper = U5MR_upper_adj)
  
  summary_results <- rbind(summary_results, summary_lquad_adj)
  
  names(res_log_quad_adj) <- c("period_1", "period_2", "period_3", "period_4", "period_5")
  
  
  #------------------------------------------------
  
  list_res <- list(surv_lnorm = res_surv_lnorm, surv_lnorm_adjusted = res_surv_lnorm_adjusted, 
                   surv_weibull = res_surv_weibull, discrete_hazard = res_discrete_hazard, log_quad = res_log_quad,
                   log_quad_adjusted = res_log_quad_adj)
  
  
  
  #Save the list somewhere (Probably new folder) STILL NEED RELATIVE PATH
  print(i)
  saveRDS(list_res, str_c("../models/", str_replace(all_files[i], "sim", "res")))
  
  #Save the summary results somewhere (Probably new folder) STILL NEED RELATIVE PATH
  
  saveRDS(summary_results, str_c("../summaries/", str_replace(all_files[i], "sim", "summary")))
}

#Run the for loop above.

#Then use warnings() to look through the warnings output from logquad. 

#inputs_row_51 is a list of each periods' input data, which comes from `format_data` and is used in `lagrange5q0`

#res_row_51 is a list of each periods' results from `lagrange5q0`. Period 3 is the only one with uncertainty measures.

res_row_51[[3]]
