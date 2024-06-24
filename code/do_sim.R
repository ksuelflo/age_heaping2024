source("sim_functions.R")

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_range1_periodLength_periods_proportion1_seed

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
# "sim_ROW=1_numSIM=1_lnormmean=15_lnormsd=2_samplesize=1000_age1=12_range1=921_periodLength=12_periods=5_proportion1=01_seed=1"



all_files <- list.files(path = "/Users/kylesuelflow/Macalester-Stuff/Research-Taylor/age_heaping2024/data")

sim <- readRDS(str_c("../data/", all_files[2]))
view(sim)

a <- sim%>%
  filter(t %% 1 != 0 & t>= 9 & t <= 21)

view(a)

for (i in seq_along(all_files)){
  sim <- readRDS(str_c("../data/", all_files[i]))
  
  #fit surv_synthetic() with lognormal distribution WITHOUT adjusting for age heaping.
  
  res_surv_lnorm <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
  #------------------------------------------------
  
  #fit surv_synthetic() with lognormal distribution ADJUSTING for age heaping (interval censored at 9-21 months)
  
  
  #A_i needs to be updated here.
  
  add_intervals <- sim%>%
    mutate(left_interval = if_else(right_interval > 9 & right_interval <= 21, 9, left_interval),
           right_interval = if_else(right_interval > 9 & right_interval <= 21, 21, right_interval),
           across_boundary = if_else(age_at_begin > left_interval & age_at_begin < right_interval, 1, 0))
  
  res_surv_lnorm_adjusted <- surv_synthetic(df = add_intervals, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
  #------------------------------------------------
  
  #fit surv_synthetic() with Weibull distribution
  
  res_surv_weibull <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "weibull")
  
  #------------------------------------------------
  
  # reformat data to work with discrete hazards model - 6 age periods
  
  dh_data <- reformat_sims_disc_haz(sims = sim,
                         agegroup_splits = c(0,1,12,24,36,48,60))
  
  res_discrete_hazard <- fit_disc_haz(dh_data)
  
  #------------------------------------------------
  
  # get 22 age inputs for logquad model from monthly discrete hazards fit
  dh_data_monthly <- reformat_sims_disc_haz(sims = sim,
                                    agegroup_splits = 0:60)
  
  res_discrete_hazard_monthly <- fit_disc_haz(dh_data_monthly)
  
  log_quad_pars <- get_log_quad_params(res_discrete_hazard_monthly)
  
  # TO DO: Kyle - take it from here!
  
  #log quad fitting
  
  #Example format of data: 22 q(x)'s.
  data(fin1933)
  view(fin1933)
  
  #Where lower age always is 0, upper age is xq0, where x is upper age. n = upperage[i] - upperage[i-1]. rate is nqx. type is always qx
  #sex is always total. Fit is always min, besides last observation (5q0), which is match.
  
  #Calculating weights
  fin1933$weight <- c(fin1933$n[1:22]/(365.25*5), NA)
  
  #Need to have all 7 columns of `fin1933`
  input <- format_data(
    lower_age = fin1933$lower_age,
    upper_age = fin1933$upper_age,
    rate      = fin1933$rate,
    type      = fin1933$type,
    sex       = fin1933$sex,
    fit       = fin1933$fit,
    weight    = fin1933$weight)
  
  res_log_quad <- lagrange5q0(data = input)
  
  #------------------------------------------------
  
  list_res <- list(res_surv_lnorm, res_surv_lnorm_adjusted, res_surv_weibull, res_log_quad, res_discrete_hazard)
  
  saveRDS(list_res, "file_location")
}
