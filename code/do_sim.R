source("sim_functions.R")
library(tibble)

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_range1_periodLength_periods_proportion1_seed

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
# "sim_ROW=1_numSIM=1_lnormmean=15_lnormsd=2_samplesize=1000_age1=12_range1=921_periodLength=12_periods=5_proportion1=01_seed=1"


#Getting all generated data frames' file paths.
all_files <- list.files(path = "/Users/kylesuelflow/Macalester-Stuff/Research-Taylor/age_heaping2024/data")

#IMPORTANT: this data line is used to help fit logquad model.
data("fin1933")

is_bad <- rep(0, length(all_files))

for (i in 3:length(all_files)){
  
  print(all_files[i])
  
  # if (i != 2){
  #   next
  # }
  # 
  sim <- readRDS(str_c("../data/", all_files[i]))
  
  #fit surv_synthetic() with lognormal distribution WITHOUT adjusting for age heaping.
  
  # res_surv_lnorm <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
  #                                  I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
  #                                  t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
  # 
  # #Table of summary results. We will rbind to this data frame throughout the rest of the models. 
  # summary_results <- get_uncertainty_surv_synthetic(res_surv_lnorm)%>%
  #   mutate(model = "surv_synth_lnorm", .before = period)
  # 
  # print("fit lnorm model")
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
  
  #fit surv_synthetic() with Weibull distribution
  
  res_surv_weibull <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "weibull")
  
  summary_results <- rbind(summary_results, get_uncertainty_surv_synthetic(res_surv_weibull, distribution = "weibull")%>%
                             mutate(model = "surv_synth_weibull", .before = period))
  
  print("fit weibull model")
  #------------------------------------------------
  
  # reformat data to work with discrete hazards model - 6 age periods
  
  dh_data <- reformat_sims_disc_haz(sims = sim,
                         agegroup_splits = c(0,1,12,24,36,48,60))
  
  res_discrete_hazard <- fit_disc_haz(dh_data)
  
  summary_results <- rbind(summary_results, summarize_disc_haz(res_discrete_hazard)%>%
                             mutate(model = "discrete_hazards", .before = period))

  print("fit discrete hazards")
  #------------------------------------------------
  
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
  for (i in 1:5){
    input <- format_data(
      lower_age = log_quad_pars$lower_age[log_quad_pars$period == i],
      upper_age = log_quad_pars$upper_age[log_quad_pars$period == i],
      rate      = log_quad_pars$rate[log_quad_pars$period == i],
      type      = log_quad_pars$type[log_quad_pars$period == i],
      sex       = log_quad_pars$sex[log_quad_pars$period == i],
      fit       = log_quad_pars$fit[log_quad_pars$period == i],
      weight    = log_quad_pars$weight[log_quad_pars$period == i])
    
    tryCatch({
      res_log_quad[[i]] <- lagrange5q0(data = input)
      print(str_c("fit log quad on period ", i))
      
      NMR[i] <- res_log_quad[[i]]$predictions[5,4]
      IMR[i] <- res_log_quad[[i]]$predictions[16,4]
      U5MR[i] <- res_log_quad[[i]]$predictions[22,4]
      NMR_lower[i] <- res_log_quad[[i]]$predictions[5,5]
      NMR_upper[i] <- res_log_quad[[i]]$predictions[5,6]
      IMR_lower[i] <- res_log_quad[[i]]$predictions[16,5]
      IMR_upper[i] <- res_log_quad[[i]]$predictions[16,6]
      U5MR_lower[i] <- res_log_quad[[i]]$predictions[22,5]
      U5MR_upper[i] <- res_log_quad[[i]]$predictions[22,6]
    },
      error = function(e){
        message("log quad failed. Model equals NA.")
        print(e)
        NMR[i] <- NA
        IMR[i] <- NA
        U5MR[i] <- NA
        NMR_lower[i] <- NA
        NMR_upper[i] <- NA
        IMR_lower[i] <- NA
        IMR_upper[i] <- NA
        U5MR_lower[i] <- NA
        U5MR_upper[i] <- NA
      }
    )
    
  }
  
  summary_lquad <- data.frame(model = model, period = period, NMR = NMR, IMR = IMR, U5MR = U5MR,
                              NMR_lower = NMR_lower, NMR_upper = NMR_upper, IMR_lower = IMR_lower, 
                              IMR_upper = IMR_upper, U5MR_lower = U5MR_lower, U5MR_upper = U5MR_upper)
  
  summary_results <- rbind(summary_results, summary_lquad)
  
  names(res_log_quad) <- c("period_1", "period_2", "period_3", "period_4", "period_5")
  
  # input <- format_data(
  #   lower_age = log_quad_pars$lower_age,
  #   upper_age = log_quad_pars$upper_age,
  #   rate      = log_quad_pars$rate,
  #   type      = log_quad_pars$type,
  #   sex       = log_quad_pars$sex,
  #   fit       = log_quad_pars$fit,
  #   weight    = log_quad_pars$weight)
  # 
  # res_log_quad <- lagrange5q0(data = input)
  
  #------------------------------------------------
  
  list_res <- list(surv_lnorm = res_surv_lnorm, surv_lnorm_adjusted = res_surv_lnorm_adjusted, 
                   surv_weibull = res_surv_weibull, discrete_hazard = res_discrete_hazard, log_quad = res_log_quad)
  
  
  
  #Save the list somewhere (Probably new folder) STILL NEED RELATIVE PATH
  
  saveRDS(list_res, str_c("../models/", str_replace(all_files[i], "sim", "res")))
  
  #Save the summary results somewhere (Probably new folder) STILL NEED RELATIVE PATH
  
  saveRDS(list_res, str_c("../summaries/", str_replace(all_files[i], "sim", "summary")))
}

period_5 <- sim%>%
  filter(period == 5)
view(period_5)

view(log_quad_pars)
input_broken <- input
view(input_broken)
view(input)


saveRDS(input, file = "../data/input_example.rds")


#WARNINGS FROM SUCCESS

# 6: In lagrange5q0(data = input) :
#   Predicted value of k extrapolated. k < -1.1 or k > 1.5.
# 7: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.
# 8: In lagrange5q0(data = input) :
#   Predicted value of k extrapolated. k < -1.1 or k > 1.5.
# 9: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.
# 10: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.
# 11: In lagrange5q0(data = input) :
#   Predicted value of k extrapolated. k < -1.1 or k > 1.5.
# 12: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.
# 13: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.

#WARNINGS FROM FAIL

# 21: In lagrange5q0(data = input) :
#   Predicted value of k extrapolated. k < -1.1 or k > 1.5.
# 22: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.
# 23: In lagrange5q0(data = input) :
#   Predicted value of k extrapolated. k < -1.1 or k > 1.5.
# 24: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx). Prediction extrapolated.
# 25: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.
# 26: In lagrange5q0(data = input) :
#   Predicted value of k extrapolated. k < -1.1 or k > 1.5.
# 27: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx). Prediction extrapolated.
# 28: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.
# 29: In lagrange5q0(data = input) :
#   Predicted value of k extrapolated. k < -1.1 or k > 1.5.
# 30: In lagrange5q0(data = input) :
#   Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.


#DEBUGGING ERROR: Error in optim(par = init_vals, fn = optim_fn, data = df[, c("I_i", "A_i",  : 
# initial value in 'vmmin' is not finite



#row 11 had an error, so I filtered to only have rows with the same distribution.

all_files_filter <- all_files[str_detect(string = all_files, pattern = "lnormmean=15_lnormsd=22")]

#sim_bad gets an error (use all_files_filter to see file path)
sim_bad<- readRDS(str_c("../data/", all_files_filter[18]))

#sim_good runs no problem.
sim_good <- readRDS(str_c("../data/", all_files_filter[16]))
view(sim_bad)
view(sim_good)

sim_bad%>%
  group_by(period)%>%
  summarize(deaths = sum(event))


#fit surv_synthetic() with lognormal distribution WITHOUT adjusting for age heaping.

#Using the below code to test run different dfs to see if they work. The issue comes when using lognormal WITHOUT adjusting
#for age heaping.
res_surv_lnorm <- surv_synthetic(df = sim_prac, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                 I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                 t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")

#The following dfs, out of 18 (length of all_files_filter) encounter an error when running surv_synthetic.
all_files_filter[1]
all_files_filter[12]
all_files_filter[15]
all_files_filter[17]
all_files_filter[18]

all_files[5]
all_files[26]ds


sim <- readRDS(str_c("../data/", all_files[26]))
sim_5 <- readRDS(str_c("../data/", all_files[5]))

sim_5%>%
  filter(across_boundary == 1)


res_surv_lnorm <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                 I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                 t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")

add_intervals <- sim%>%
  mutate(left_interval = if_else(right_interval > 9 & right_interval <= 21, 9, left_interval),
         right_interval = if_else(right_interval > 9 & right_interval <= 21, 21, right_interval),
         across_boundary = if_else(age_at_begin > left_interval & age_at_begin < right_interval, 1, 0))%>%
  group_by(id)%>%
  mutate(across_boundary = if_else(sum(across_boundary) > 0, 1, 0))

res_surv_lnorm_adjusted <- surv_synthetic(df = add_intervals, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                          I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                          t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
res_surv_lnorm_adjusted$result
res_surv_lnorm$result

res_surv_weibull <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "weibull")




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


library(tidyr)
    
new_generate <- function(i, seed){
  
  #Setting seed for reproducibility
  seed <- seed
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
  
  # add_intervals <- cleaned_data%>%
  #   mutate(left_interval = if_else(right_interval > 9 & right_interval <= 21, 9, left_interval),
  #          right_interval = if_else(right_interval > 9 & right_interval <= 21, 21, right_interval),
  #          across_boundary = if_else(age_at_begin > left_interval & age_at_begin < right_interval, 1, 0))%>%
  #   group_by(id)%>%
  #   mutate(across_boundary = if_else(sum(across_boundary) > 0, 1, 0))
  # 
  # res_surv_lnorm_adjusted <- surv_synthetic(df = add_intervals, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
  #                                           I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
  #                                           t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
  # 
  res_surv_lnorm <- surv_synthetic(df = cleaned_data, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal", init_vals = NA)
}

for (j in 1:20){
  print(j)
  new_generate(11, j)
}


sim_5 <- readRDS(str_c("../data/", all_files[5]))

res_surv_lnorm <- surv_synthetic(df = sim_5, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                 I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                 t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal", init_vals = c(rep(exp(2), 5),15,15,15,15,15))




new_sim_row_11 <- new_generate(11)

all_files[5]




res_surv_lnorm <- surv_synthetic(df = sim_5, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                 I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                 t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")

sim_5

view(new_sim_row_11)









