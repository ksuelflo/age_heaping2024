# Project: Age-Heaping Simulation
# Authors: Kyle Suelflow & Taylor Okonek
# Date: Summer 2024


## Load necessary packages up front
library(survival)
library(dplyr)
library(flexsurv)
library(knitr)
library(pssst)
library(logquad5q0)
library(SUMMER)


# Functions ---------------------------------------------------------------

#' @description
#' This function is used to generate the probability of a child dying, given the parameters provided. It is used within the 
#' `general_sim()` function.
#' 
#' @inheritParams general_sim()
#' @param max_age is a vector of the maximum age a child can live for each child at each period.
#' @param age_at_begin is a vector of the age at the beginning of each period for each child.
#' 
#' @return a vector of probabilities of death for each child for each period.
#' @author Kyle Suelflow
helper_fill_p <- function(distribution, param_matrix, max_age, age_at_begin){
  
  # p is the vector we will populate. It is the same length as max_age, age_at_begin, etc. (We will be putting it in a data frame with them!)
  periods <- nrow(param_matrix)
  curr_period <- 1
  p <- rep(0, length(max_age))
  
  for (i in seq_along(age_at_begin)){
    
    if (distribution == "weibull"){
      #Conditional probability p: p(death now till max_age in time period) - p(death from 0 - age_at_begin) / p(survived up until age_at_begin.)
      
      p[i] <- (pweibull(q = max_age[i], shape = param_matrix[curr_period,2], scale = param_matrix[curr_period,1]) - pweibull(q = age_at_begin[i], shape = param_matrix[curr_period,2], scale = param_matrix[curr_period,1]))/
        (1 - pweibull(q = age_at_begin[i], shape = param_matrix[curr_period,2], scale = param_matrix[curr_period,1]))
      
    }
    
    else if (distribution == "lognormal"){
      #Conditional probability p: p(death now till max_age in time period) - p(death from 0 - age_at_begin) / p(survived up until age_at_begin.)
      
      p[i] <- (plnorm(q = max_age[i], meanlog = param_matrix[curr_period, 1], sdlog = param_matrix[curr_period, 2]) - plnorm(q = age_at_begin[i], meanlog = param_matrix[curr_period, 1], sdlog = param_matrix[curr_period, 2]))/
        (1 - plnorm(q = age_at_begin[i], meanlog = param_matrix[curr_period, 1], sdlog = param_matrix[curr_period, 2]))
    }
    
    if (curr_period == periods){
      curr_period <- 0
    }
    curr_period <- curr_period + 1
    
  }
  return (p)
  
}



#' @description
#' This function simulates the mortality of children starting at age 0.
#' 
#' @param distribution Takes values of either "weibull" or "lognormal".
#' @param param_matrix Takes in a matrix, which has 3 columns. The first column is for `mulog` or `shape`, and the second column is for
#' `sdlog` or `scale`. These are the parameters for the "lognormal" or "weibull" distributions, respectively. The 3rd column is for the
#' number of months in each period. Each row represents a period. 
#' @param num_child Is the number of children that are born on the first day of each month within each period. 
#' 
#' @return A data frame where each row corresponds to a child, with information on whether or not they died in the time period, amongst
#' other things.
#' 
#' @author Kyle Suelflow

general_sim <- function(num_child, param_matrix, distribution = "weibull"){
  
  #NOTE: ASSUMING CONSTANT LENGTH OF PERIODS.
  
  
  #Setting up data frame. The first five plus `age_at_begin` do not iteratively change. `t` and `event` get populated using for loop below.
  
  months <- param_matrix[,3][1]
  
  periods <- nrow(param_matrix)
  # ids <- rep(1:(num_child*periods*months), each = periods)
  ids <- rep(1:(num_child*months), each = periods)
  # birthdates <- rep(1:(months*periods), each = num_child*periods)
  birthdates <- rep(1:(months*periods), each = num_child)
  # period <- rep(1:periods, times = num_child*periods*months)
  period <- rep(1:periods, times = num_child*months)
  # max_age <- (months*period - birthdates) + 1
  # max_age <- rep(0, length(period))
  max_age <- (months*period) - (birthdates) + 1
  
  # for (i in seq_along(period)){
  #   max_age[i] <- (months*period) - (birthdates[i] * period[i]) + 1
  #   total_months <- sum(months[1:period[i]])
  #   max_age[i] <- total_months - birthdates[i] + 1
  # }
  
  
  max_age[max_age < 0] <- 0
  age_at_begin <- c(0,max_age)
  age_at_begin[seq(periods+1, length(age_at_begin), periods)] <- 0
  age_at_begin <- age_at_begin[-length(age_at_begin)]
  
  #FOR NO CENSORING: This sets the max age for the last time period to be infinite (100000000).
  # max_age[period == periods] <- Inf
  
  # Note from Taylor - put in a dataframe, easier to debug
  sim_df <- data.frame(id = ids,
                       birthdate = birthdates,
                       period = period,
                       max_age = max_age,
                       age_at_begin = age_at_begin)
  
  # Filling p based on distribution.
  p <-helper_fill_p(distribution = distribution, param_matrix = param_matrix, max_age = max_age, age_at_begin = age_at_begin)
  
  # Note from Taylor - add to dataframe, easier to debug
  sim_df$p <- p
  
  # This column (age_at_begin) is useful for the for loop: Essentially this allows me to "group_by" 
  # by checking if the age is 0. If it is, then I know it is a new child and to reset the boolean `already_died`
  
  # Note from Taylor - should this be all zeros? Could you just do this the same way as you 
  # define t and event?
  
  
  # these two columns will be iteratively populated in the for loop. T is time of event within period, 
  # event is a binary outcome, either 1 or 0.
  
  t <- rep(0, num_child*sum(months)*periods)
  event <- rep(0, num_child*sum(months)*periods)
  
  # for loop that runs binomial using probability p, and if they die, then simulates an appropriate time. 
  # If not, t is equal to max age and event is equal to 0. 
  
  already_died <- FALSE
  curr_period <- 1
  
  start <- Sys.time()
  
  for (i in seq_along(ids)){
    
    # If the child already died in an earlier period, skip to the else statement. 
    # Otherwise, go through this if statement.
    
    if (!already_died){
      
      died <- ifelse(p[i] <= 0, 0, rbinom(n = 1,size = 1, p = p[i]))
      
      #If it is the last period, then everyone has to die! FOR NO CENSORING
      # if (curr_period == periods){
      #   died <- 1
      # }
      
      # If they died, sample the time of death until it falls within their age at the 
      # beginning of the period and the age at the end of the period. 
      
      if (died == 1) {
        event[i] <- died
        
        # The following series of if, else if, else statements handle which distribution we choose. 
        # For now, only "weibull" works. We randomly sample a death. 
        
        if (distribution == "weibull"){
          time <- rweibull(n = 1, scale = param_matrix[curr_period, 1], shape = param_matrix[curr_period,2])
        }
        
        else if (distribution == "lognormal"){
          time <- rlnorm(n = 1, meanlog = param_matrix[curr_period, 1], sdlog = param_matrix[curr_period,2])
        }
        
        #We continue randomly sampling deaths until the time of death falls within the age at beginning and max_age of the time period. 
        while (time >= max_age[i] | time < age_at_begin[i]){
          
          if (distribution == "weibull"){
            time <- rweibull(n = 1, scale = param_matrix[curr_period, 1], shape = param_matrix[curr_period,2])
          }
          
          
          else if (distribution == "lognormal"){
            time <- rlnorm(n = 1, meanlog = param_matrix[curr_period, 1], sdlog = param_matrix[curr_period,2])
          }
          
        }
        
        #Set already_died to true, so we know not to use rbinom to see if the child died in a future period, since they are already dead.
        already_died <- TRUE
      }
      
      #If they did not die, set the time of event to their maximum age within the period, and set the event to 0.
      else{
        time <- max_age[i]
        event[i] <- 0
      }
      
    }
    
    #If they already died, set their time of event to their maximum age within the period, 
    # and set the event to NA. This allows us to filter out the NA afterwards, 
    # because these observations aren't possible! (the child is dead!)
    else{
      time <- max_age[i]
      event[i] <- NA
      # event[i] <- 1
    }
    
    t[i] <- time
    curr_period <- curr_period + 1
    
    # In the if statement below this one, we check if the next observation has an age_at_begin of 0, 
    # essentially checking if the next observation is a new child. This if statement is included so there 
    # is not an index out of bounds error on the final observation.
    # if (i == length(max_age)){
    #   break
    # }
    
    # Checks if the next observation is a new child. If it is, the boolean already_died is reset.
    if (period[i] == periods){
      already_died <- FALSE
      curr_period <- 1
    }
    
  }
  end <- Sys.time()
  print(start-end)
  
  sim_df$t <- t
  sim_df$event <- event
  
  return(sim_df)
}

#' @description
#' This function simulates the phenomenon of age heaping in child mortality. It takes in a data frame of child mortality data, generated
#' from the `general_sim` function. It then randomly chooses deaths around specific ages of death to "heap" at a round number. This is
#' usually at 3, 6, 12, or 24 months. The interval surrounding the age, the age itself, and the probability to be "chosen" to be included in
#' the heap are all parameters, `range_heap`, `ages_at_heaping`, `proportion_heap` respectively.
#' 
#' @param ages_at_heaping a vector of ages, in months, that represent ages where deaths will be artificially clumped at. 
#' @param proportion_heap a vector of proportions. Each index is the probability that a death will be heaped at the `ages_at_heaping` age.
#' @param range_heap a matrix where each row corresponds to one age where deaths will be heaped. The 1st column of that row is the bottom
#' of the range, and 2nd column is the top of the range. So if the row is (2,3) and the age we heap at is 6, then the full range is 4,9. 
#' Children who die in this range will be eligible for the ages to be heaped to 6 months.
#' @param sim_data A data frame which results from calling `general_sim()`.
#' 
#' @returns A data frame, which is a modified version of `sim_data`, where some deaths have been moved to be heaped at certain ages. 
#' @author Kyle Suelflow
sim_age_heap <- function(ages_at_heaping, proportion_heap, range_heap, sim_data){
  
  #Time and event for each observation. 
  times <- sim_data%>%pull(t)
  events <- sim_data%>%pull(event)
  #Emptu container with which to put intervals into. 
  intervals <- vector("list", length = length(ages_at_heaping))
  #empty vector, which will be populated with a binomial distribution: Either the observation will have their time at event heaped, or
  #it won't. 
  indicators <- rep(0, length(times))
  
  # #Looping through each age to build out the interval at which observations are eligible to be heaped. 
  # for (i in seq_along(ages_at_heaping)){
  #   intervals[[i]] <- c(ages_at_heaping[i] - range_heap[i,1], ages_at_heaping[i] + range_heap[i,2])
  # }
  # 
  #looping through each observation.
  for (i in seq_along(times)){
    
    #If they didn't die, then we don't need to worry about age heaping, and can move to the next observation.
    if (is.na(events[i]) | events[i] == 0){
      next
    }
    
    #Looping through each age at which we heap observations. 
    for (j in seq_along(ages_at_heaping)){
      
      if (is.na(ages_at_heaping[j])){
        next
      }

      #If the time of death is within the interval to be eligible to be heaped:
  
      if (range_heap[j,][1] <= times[i] & times[i] <= range_heap[j,][2]){
        
        #randomly sample the chance of being heaped
        indicators[i] <- rbinom(1,1,proportion_heap[j])
        
        #Heap them if they were sampled, otherwise keep their time of death the same.
        sim_data$t[i] <- ifelse(indicators[i] == 1, ages_at_heaping[j], times[i])
        break
      }
      
    }
  }
  sim_data <- sim_data%>%
    mutate(old_times = times,
           binom_heap = indicators)
  return (sim_data)
}

#'@description
#'This function takes in a data frame, coming from `general_sim()` or `sim_age_heap()`, and adds columns so that the returned data frame
#'can be used in the `pssst` package to fit the `surv_synthetic()` function found in that package. One of the main things this function does is to
#'interval censor all observations, following standard DHS guidelines for doing so.
#'
#'@param sims a data frame, which comes from either `general_sim()` or `sim_age_heap()`.
#'
#'@returns A data frame which has new columns `left_interval`, `right_interval`, `interval_indicator`, `right_censor_age`, `period_length`, and
#'`across_boundary`, in addition to all the original columns in `sims`.
#'
#'@author Kyle Suelflow
clean_data <- function(sims){
  
  clean_df <- sims%>%
    mutate(event = if_else(t>60, 0, event),
           t = if_else(t>60, 60, t))%>%
    mutate(left_interval = case_when(t <= 1 & event == 1 ~ t,
                                     t == floor(t) & event == 1 & t<=24 ~ t, #Handle 6 and 12 month heaping
                                     t > 1 & t <= 24 & event == 1 ~ floor(t),
                                     t > 24 & t <= 36 & event == 1 ~ 24,
                                     t > 36 & t <= 48 & event == 1 ~ 36,
                                     t > 48 & event == 1~ 48, .default = 1000),
           right_interval = case_when(t <= 1& event == 1 ~ t,
                                      t == ceiling(t) & event == 1 & t<=24 ~ t+ 1, #Handle 6 and 12 month heaping
                                      t > 1 & t <= 24 & event == 1~ ceiling(t),
                                      t > 24 & t <= 36 & event == 1~ 36,
                                      t > 36 & t <= 48 & event == 1~ 48,
                                      t > 48 & event == 1~ 60, .default = 1000))%>%
    mutate(event = replace_na(event, 1))%>%
    group_by(id)%>%
    mutate(interval_indicator = if_else(sum(event) >= 1, 1, 0),
           right_censor_age = if_else(interval_indicator == 0, max(t), 1000),
           left_interval = if_else(sum(left_interval) < 5000, min(left_interval), 1000),
           right_interval = if_else(sum(right_interval) < 5000, min(right_interval), 1000))%>%
    ungroup()%>%
    mutate(period_length = sims$max_age[1], #Period_length is equal to the max age in first time period for child born at the beginning.
           across_boundary = if_else(age_at_begin > left_interval & age_at_begin < right_interval, 1, 0))%>%#Fix this to be generalizable once figure out how params df works.
    mutate(age_at_begin = (period - 1)*period_length - birthdate + 1)%>%
    group_by(id)%>%
    mutate(across_boundary = if_else(sum(across_boundary) > 0, 1, 0))%>%
    ungroup()%>%
    mutate(left_interval = if_else(across_boundary == 1 & age_at_begin > left_interval & t%%1 != 0, age_at_begin, left_interval))%>%
    mutate(right_interval = if_else(across_boundary == 1 & age_at_begin < left_interval & t%%1 != 0, max_age, right_interval))%>%
    group_by(id)%>%
    mutate(left_interval = if_else(length(unique(left_interval)) > 1, max(left_interval), min(left_interval)))%>%
    mutate(right_interval = if_else(length(unique(right_interval)) > 1, min(right_interval), max(right_interval)))%>%
    mutate(across_boundary = if_else(age_at_begin > left_interval & age_at_begin < right_interval, 1, 0))
  
  return (clean_df)
}

#' @description This function takes in a data frame, coming from `clean_data()` and
#' reformats it so that the discrete hazards model can be fit to it. Much of the code
#' was cribbed from `getBirths()` in the SUMMER package in R.
#'
#' @param sims a data frame, which comes from `clean_data()`
#' @param agegroup_splits a vector, denoting how age groups should be split up. If
#' fitting the traditional discrete hazards model from Mercer et all, agegroup_splits should be 
#' c(0,1,12,24,36,48,60). If you want monthly hazards (which will be more relevant)
#' for inputting this data into the log-quad model, specify 0:60.
#'
#' @returns A data frame which has binomial counts of deaths (conditional on 
#' survival in each time period) by time period and age group. New columns 
#' have names `period`, `agegroup`, `total`, and `died`.
#' @author Taylor Okonek
reformat_sims_disc_haz <- function(sims, 
                                   agegroup_splits = c(0,1,12,24,36,48,60)) {
  
  # Infer number of time periods from rows of observations per id
  n_periods <- sims %>% filter(id == 1) %>% nrow()
  
  # Assumes all periods are the same length
  period_length <- sims$period_length[1]
  
  df <- as.data.frame(sims)
  
  # only need 1 row per person, and all relevant information needed is identical across time periods
  df <- df %>% distinct(id, .keep_all = TRUE)
  
  # if t < 1, interval censor it from 0 to 1 months rather than exactly observed in accordance with DH model
  df$exact_died <- df$left_interval == df$right_interval
  
  df <- df %>%
    mutate(left_interval = ifelse(t < 1 & t >= 0 & exact_died, 0, left_interval),
           right_interval = ifelse(t < 1 & t >= 0 & exact_died, 1, right_interval))
  
  # make obsStart = df$birthdate
  df$obsStart <- df$birthdate
  df$obsStart <- df$obsStart - 1
  
  # make obsStop = right_interval + birthdate if interval censored, otherwise birthdate + right_censor_age
  df$obsStop <- ifelse(df$interval_indicator == 1, df$obsStart + df$right_interval,
                       df$obsStart + df$right_censor_age)
  df$died <- df$interval_indicator
  df$dob <- df$obsStart 
  
  # remove extraneous columns for space
  df_mini <- df %>%
    dplyr::select(id, died, obsStart, obsStop, dob) 
  
  # cribbed from SUMMER::getBirths
  
  test <- survSplit(Surv(time = obsStart, time2=obsStop, 
                         event = died, origin = dob)~dob+died+id, 
                    data = df_mini, 
                    cut = 1:60, 
                    start = "agemonth", 
                    end = "tstop", 
                    event = "interval_indicator")
  
  # assign time periods to each row
  test$age_at_tstop <- test$dob + test$tstop
  
  test$period <- cut(test$age_at_tstop, breaks = seq(0,period_length * (n_periods + 1), by = period_length))
  
  levels(test$period)
  
  test <- test %>%
    mutate(year = case_when(period == levels(test$period)[1] ~ 1,
                            period == levels(test$period)[2] ~ 2,
                            period == levels(test$period)[3] ~ 3,
                            period == levels(test$period)[4] ~ 4,
                            .default = 5))
  
  
  suppressMessages(dh_df <- test %>%
                     mutate(agegroup = cut(agemonth, breaks = agegroup_splits, right = FALSE)) %>%
                     group_by(year, agegroup) %>%
                     summarize(total = n(),
                               died = sum(interval_indicator)))
  
  return(dh_df)
}

#' @description This function takes in a data frame, coming from `reformat_sims_disc_haz()` and
#' fits the discrete hazards model described in Mercer et al. (2015) to the data.
#'
#' @param df a data frame, which comes from `reformat_sims_disc_haz()`
#'
#' @returns A list of `glm` fits containing logistic regression output for each time period
#' @author Taylor Okonek
fit_disc_haz <- function(df) {
  mod_list <- list()
  for (i in 1:length(unique(df$year))) {
    temp <- df %>% filter(year == i)
    
    mod_list[[i]] <- glm(cbind(died, total - died) ~ -1 + agegroup, data = temp, family = binomial(link = "logit"))
    
    # betas <- summary(mod)$coef[, 1]
    # probs <- expit(betas)
    # ns <- c(1,11,12,12,12,12)
    # mean.est <- (1 - prod((1 - probs)^ns, na.rm = TRUE))
    # (1 - cumprod((1 - probs)^ns))
  }
  
  return(mod_list)
}

#' @description This function takes in the output from `fit_disc_haz()` where
#' monthly age groups are specified, and obtains
#' the 22 age inputs needed for the log quad model at the appropriate age groups.
#' Hazard is assumed constant over first month of life, and q(x) for ages x = 7 days,
#' 14 days, and 21 days is computed accordingly.
#'
#' @param disc_haz_res a list, which comes from `fit_disc_haz()`
#'
#' @returns A data frame containing the upper_age and rate inputs needed for `format_data` 
#' and `lagrange5q0` for the log-quad model in each time period
#' @author Taylor Okonek
get_log_quad_params <- function(disc_haz_res) {
  
  upper_age <- c(7, 14, 21, 28, 60.8750 ,  91.3125,  121.7500  ,152.1875,
                 182.6250 , 213.0625,  243.5000 , 273.9375  ,304.3750 , 334.8125,  365.2500 , 456.5625,
                 547.8750 , 639.1875 , 730.5000 ,1095.7500 ,1461.0000 ,1826.2500 ,1826.2500)
  
  # set up return dataframe
  ret_df <- data.frame(period = NA, upper_age = NA, rate = NA)
  
  # loop through time periods
  
  for (i in 1:length(disc_haz_res)) {
    
    mod <- disc_haz_res[[i]]
    betas <- summary(mod)$coef[, 1]
    probs <- expit(betas)
    cumprob <- unname(1 - cumprod((1 - probs)))
    rate <- cbind(c(cumprob[1]/4,
                    cumprob[1]/2,
                    cumprob[1]/4 * 3,
                    cumprob[c(1:12,15,18,21,24,36,48,60,60)]))
    
    ret_df <- rbind(ret_df, data.frame(period = i,
                                       upper_age = upper_age, 
                                       rate = rate))
  }
  
  return(ret_df[-1,])
}

#' @description This function takes in the output from `fit_disc_haz()` where
#' the classical six age groups are specified, and obtains
#' estimates of NMR, IMR, U5MR, and 95% CI's for each
#'
#' @param res_list a list, which comes from `fit_disc_haz()`
#'
#' @returns A data frame containing estimates of NMR, IMR, and U5MR from the 
#' discrete hazards model, along with confidence intervals for each summary measure
#' @author Taylor Okonek
summarize_disc_haz <- function(res_list) {
  
  # get number of periods
  n_periods <- length(res_list)
  
  # Get confidence band for discrete model
  # Code cribbed from getDirect from SUMMER
  
  # set up lims_mat to contain lower and upper bound of confidence intervals 
  lims_mat <- matrix(NA, nrow = length(res_list[[1]]$coefficients), ncol = 2)
  
  # set up return dataframe for NMR, IMR, U5MR 
  ret_df <- data.frame(period = 1:n_periods,
                       NMR = NA,
                       IMR = NA,
                       U5MR = NA,
                       NMR_lower = NA,
                       NMR_upper = NA,
                       IMR_lower = NA,
                       IMR_upper = NA,
                       U5MR_lower = NA,
                       U5MR_upper = NA
  )
  
  # loop through periods
  for (j in 1:n_periods) {
    
    # get point estimates
    betas <- summary(res_list[[j]])$coef[, 1]
    probs <- expit(betas)
    ns <- c(1,11,12,12,12,12)
    mean_ests <- (1 - cumprod((1 - probs)^ns))
    
    # put mean estimates into return dataframe
    ret_df[j,c("NMR","IMR","U5MR")] <- mean_ests[c(1,2,6)]
    
    # loop through age groups (cumulatively)
    for (i in 1:nrow(lims_mat)) {
      which_vals <- 1:i
      V2 <- stats::vcov(res_list[[j]])[which_vals, which_vals]
      V <- V2
      betas2 <- summary(res_list[[j]])$coef[which_vals, 1]
      betas <- betas2
      probs <- expit(betas)
      ns <- c(1,11,12,12,12,12)[which_vals]
      mean.est <- (1 - prod((1 - probs)^ns, na.rm = TRUE))
      gamma <- prod((1 + exp(betas))^ns, na.rm = TRUE)
      derivatives <- (gamma)/(gamma - 1) * ns * expit(betas)
      derivatives[which(is.na(derivatives))] <- 0
      var.est <- t(derivatives) %*% V %*% derivatives
      lims_mat[i,] <- logit(mean.est) + stats::qnorm(c(0.025, 0.975)) * 
        sqrt(c(var.est))
    }
    
    # get confidence intervals
    cis <- expit(lims_mat)
    
    # put CIs into return dataframe
    ret_df[j,c("NMR_lower","NMR_upper")] <- cis[1,]
    ret_df[j,c("IMR_lower","IMR_upper")] <- cis[2,]
    ret_df[j,c("U5MR_lower","U5MR_upper")] <- cis[6,]
    
  }
  
  return(ret_df)
  
}

#' @description This function takes in the output from `surv_synthetic()` where
#' and obtains confidence intervals for
#' estimates of NMR, IMR, U5MR in each period
#'
#' @param res_surv_synthetic a list of output, which comes from `surv_synthetic()`
#' @param distribution a string referring to which distribution was used in
#' `surv_synthetic`. Currently supports "lognormal" and "weibull", and defaults to 
#' "lognormal".
#'
#' @returns A data frame containing estimates of NMR, IMR, and U5MR from the 
#' discrete hazards model, along with confidence intervals for each summary measure
#' @author Taylor Okonek
get_uncertainty_surv_synthetic <- function(res_surv_synthetic,
                                           distribution = "lognormal") {
  
  # extract covariance matrix and parameter estimates to sample from multivariate normal
  means <- res_surv_synthetic$optim$par
  var_est <- res_surv_synthetic$variance
  
  # sample from multivariate normal
  samps <- mvrnorm(n = 1000, mu = means, Sigma = var_est)
  
  # for each sample, compute estimates of NMR, IMR, U5MR 
  # one matrix of samples per summary measures, columns by period
  n_periods <- nrow(res_surv_synthetic$result)
  nmr_samps <- matrix(NA, nrow = 1000, ncol = n_periods)
  imr_samps <- matrix(NA, nrow = 1000, ncol = n_periods)
  u5mr_samps <- matrix(NA, nrow = 1000, ncol = n_periods)
  
  
  if (distribution == "lognormal") {
    
    for (i in 1:nrow(samps)) {
      # for lognormal, c(log_sigma_mean for each period, log_1overmu_mean for each period)
      log_sigma_ests <- samps[i,][1:n_periods]
      log_1overmu_ests <- samps[i,][(n_periods + 1):(n_periods * 2)]
      
      mu_ests <- 1/exp(log_1overmu_ests)
      sigma_ests <- exp(log_sigma_ests)
      
      nmr_samps[i,] <- plnorm(1, meanlog = mu_ests, sdlog = sigma_ests)
      imr_samps[i,] <- plnorm(12, meanlog = mu_ests, sdlog = sigma_ests)
      u5mr_samps[i,] <- plnorm(60, meanlog = mu_ests, sdlog = sigma_ests)
      
    }
    
  } else if (distribution == "weibull") {
    
    for (i in 1:nrow(samps)) {
      # for weibull, c(log_shape_mean for each period, log_scale_mean for each period)
      log_shape_ests <- samps[i,][1:n_periods]
      log_scale_ests <- samps[i,][(n_periods + 1):(n_periods * 2)]
      
      nmr_samps[i,] <- pssst:::p_weibull(1, log_shape = log_shape_ests, log_scale = log_scale_ests)
      imr_samps[i,] <- pssst:::p_weibull(12, log_shape = log_shape_ests, log_scale = log_scale_ests)
      u5mr_samps[i,] <- pssst:::p_weibull(60, log_shape = log_shape_ests, log_scale = log_scale_ests)
    }
    
  } else {
    stop("distribution not supported")
  }
  
  # create return dataframe
  ret_df <- data.frame(period = 1:n_periods,
                       NMR = res_surv_synthetic$result$NMR,
                       IMR = res_surv_synthetic$result$IMR,
                       U5MR = res_surv_synthetic$result$U5MR,
                       NMR_lower = apply(nmr_samps, 2, quantile, 0.025),
                       NMR_upper = apply(nmr_samps, 2, quantile, 0.975),
                       IMR_lower = apply(imr_samps, 2, quantile, 0.025),
                       IMR_upper = apply(imr_samps, 2, quantile, 0.975),
                       U5MR_lower = apply(u5mr_samps, 2, quantile, 0.025),
                       U5MR_upper = apply(u5mr_samps, 2, quantile, 0.975))
  
  return(ret_df)
}


