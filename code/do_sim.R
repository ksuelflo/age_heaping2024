source("sim_functions.R")

#General format of files:
# sim_ROW_numSIM_lnormmean_lnormsd_samplesize_age1_age2_range1_range2_periodLength_periods_proportion1_proportion2_seed

#For the first row of clean_params, and the the first simulation (out of 1000), the file looks like this:
#           "sim_1_1_15_2_1000_12_NA_9,21_NA_12_5_.01_NA_1.rds"

all_files <- list.files(path = "/Users/kylesuelflow/Macalester-Stuff/Research-Taylor/age_heaping2024/data")

sim <- readRDS(str_c("../data/", all_files[1]))
view(sim)
for (i in seq_along(all_files)){
  sim <- readRDS(str_c("../data/", all_files[i]))
  
  #fit surv_synthetic() with lognormal distribution WITHOUT adjusting for age heaping.
  
  res_surv_lnorm <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
  #------------------------------------------------
  
  #fit surv_synthetic() with lognormal distribution ADJUSTING for age heaping (interval censored at 9-21 months)
  
  add_intervals <- sim%>%
    mutate(left_interval = if_else(right_interval > 9 & right_interval <= 21, 9, left_interval),
           right_interval = if_else(right_interval > 9 & right_interval <= 21, 21, right_interval))
  
  res_surv_lnorm_adjusted <- surv_synthetic(df = add_intervals, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "lognormal")
  #------------------------------------------------
  
  #fit surv_synthetic() with Weibull distribution
  
  res_surv_weibull <- surv_synthetic(df = sim, individual = "id", survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = "period_length",
                                   I_i = "interval_indicator", A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval",
                                   t_1i = "right_interval", numerical_grad = TRUE, dist = "weibull")
  
  #------------------------------------------------
  
  #log quad fitting
  
  #Example format of data: 22 q(x)'s.
  data(fin1933)
  view(fin1933)
  
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
  
  #fit discrete hazards
  
  res_discrete_hazard <- NA
  
  list_res <- list(res_surv_lnorm, res_surv_lnorm_adjusted, res_surv_weibull, res_log_quad, res_discrete_hazard)
  
  saveRDS(list_res, "file_location")
}
