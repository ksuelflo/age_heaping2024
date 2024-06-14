# Project: Age-Heaping Simulation
# Authors: Kyle Suelflow & Taylor Okonek
# Date: Summer 2024


## Load necessary packages up front
library(survival)
library(tidyverse)
library(flexsurv)
library(devtools)
library(knitr)
library(pssst)



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

general_sim <- function(num_child, param_matrix, distribution = "weibull"){
  
  #Setting up data frame. The first five plus `age_at_begin` do not iteratively change. `t` and `event` get populated using for loop below.
  
  months <- param_matrix[,3]
  periods <- nrow(param_matrix)
  # ids <- rep(1:(num_child*periods*months), each = periods)
  ids <- rep(1:(sum(num_child*months)), each = periods)
  # birthdates <- rep(1:(months*periods), each = num_child*periods)
  birthdates <- rep(1:(sum(months)), each = num_child*periods)
  # period <- rep(1:periods, times = num_child*periods*months)
  period <- rep(1:periods, times = num_child*sum(months))
  # max_age <- (months*period - birthdates) + 1
  max_age <- rep(0, length(period))
  
  for (i in seq_along(period)){
    total_months <- sum(months[1:period[i]])
    max_age[i] <- total_months - birthdates[i] + 1
  }
  
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
sim_age_heap <- function(ages_at_heaping, proportion_heap, range_heap, sim_data){
   
  #Time and event for each observation. 
  times <- sim_data%>%pull(t)
  events <- sim_data%>%pull(event)
  #Emptu container with which to put intervals into. 
  intervals <- vector("list", length = length(ages_at_heaping))
  #empty vector, which will be populated with a binomial distribution: Either the observation will have their time at event heaped, or
  #it won't. 
  indicators <- rep(0, length(times))
  
  #Looping through each age to build out the interval at which observations are eligible to be heaped. 
  for (i in seq_along(ages_at_heaping)){
    intervals[[i]] <- c(ages_at_heaping[i] - range_heap[i,1], ages_at_heaping[i] + range_heap[i,2])
  }
  
  #looping through each observation.
  for (i in seq_along(times)){
    
    #If they didn't die, then we don't need to worry about age heaping, and can move to the next observation.
    if (events[i] == 0){
      next
    }
    
    #Looping through each age at which we heap observations. 
    for (j in seq_along(ages_at_heaping)){
      
      #If the time of death is within the interval to be eligible to be heaped:
      if (intervals[[j]][1] <= times[i] & times[i] <= intervals[[j]][2]){
        
        #randomly sample the chance of being heaped
        indicators[i] <- rbinom(1,1,proportion_heap[j])
        
        #Heap them if they were sampled, otherwise keep their time of death the same.
        times[i] <- ifelse(indicators[i] == 1, ages_at_heaping[j], times[i])
        break
      }
      
    }
  }
  sim_data <- sim_data%>%
    mutate(new_times = times,
           binom_heap = indicators)
  return (sim_data)
}

#'@description
#'Recovers the parameters for a specified distribution using the `survival` package and the `flexsurv` package. This function cleans the 
#'data, and then recovers the parameters.
#'
#'@param sims A data frame generated from the `general_sim()` function. The output from `general_sim()` can be put directly into this
#'function.
#'@param parameters A matrix of the parameters used to generate `sims`. Used in this function to be columns in the returned data frame, 
#'as a means to compare the recovered parameters to the true parameters.
#'@param distribution A distribution: can either be Weibull, lognormal, or gengamma.
#'
#'@returns A data frame with the recovered parameters and the true parameters as columns, and each row representing a unique period. 
recover_params <- function(sims, parameters, distribution){
  
  periods <- nrow(parameters)
  period <- 1:periods
  param_1 <- rep(0, periods)
  param_2 <- rep(0, periods)
  
  initial_clean_all <- sims%>%
    mutate(event = if_else(t>60, 0, event),
           t = if_else(t>60, 60, t))%>%
    mutate(left_interval = case_when(t <= 1 ~ t,
                                     t > 1 & t <= 24 ~ floor(t),
                                     t > 24 & t <= 36 ~ 24,
                                     t > 36 & t <= 48 ~ 36,
                                     t > 48 ~ 48),
           right_interval = case_when(t <= 1 ~ t,
                                      t > 1 & t <= 24 ~ ceiling(t),
                                      t > 24 & t <= 36 ~ 36,
                                      t > 36 & t <= 48 ~ 48,
                                      t > 48 ~ 60))
    # mutate(left_interval = )
  
  inital_clean_survreg <- initial_clean_all%>%
    filter(!is.na(event),
           t != 0, 
           period == 1,
           age_at_begin<60)
  
  for (i in 1:nrow(parameters)){
    #cleaned data: filters to only the current period, removed observations which have already died, and removed observations which
    #haven't been born yet. 
    # curr_period_data <- initial_clean%>%
    #   filter(period == i,
    #          age_at_begin < 60)
    # 
    
    #If it is the first time period, we don't need to use `flexsurv`, instead just using the `survival` package. 
    # if (i == 1){
    #   res <- survreg(formula = Surv(time = left_interval, time2 = right_interval, event = event, type = "interval") ~ 1, data = curr_period_data, dist = distribution)
    #   
    #   if (distribution == "lognormal"){
    #     param_1[i] <- coef(res)
    #     param_2[i] <- res$scale
    #   }
    #   
    #   else if (distribution == "weibull"){
    #     param_1[i] <- log(1/res$scale)
    #     param_2[i] <- coef(res)   
    #   }
    #   
    # }
    
    #If it is not the first time period, then we need to use `flexsurv`. 
    else{
      # res <- surv_synthetic(curr_period_data, survey = FALSE, p = "period", a_pi = "age_at_begin", l_p = 60, I_i = 1, A_i = 1, t_i = 60, 
      #                       t_0i = "left_interval", t_1i = "right_interval", dist = distribution, individual = id)
      res <- flexsurvreg(formula = Surv(time = age_at_begin, time2 = t, event = event) ~ 1, data = curr_period_data, dist = distribution)
      param_1[i] <- res$coefficients[1]
      param_2[i] <- exp(res$coefficients[2])
    }
    
  }
  
  if (distribution == "weibull"){
    return (data.frame(period = period, rec_shape = param_1, shape = log(parameters[,2]), rec_scale = param_2, scale = log(parameters[,1])))
  }
  else if (distribution == "lognormal"){
    return (data.frame(period = period, rec_mu = param_1, mu = parameters[,1], rec_sigma = param_2, sigma = parameters[,2]))
  }

}

#NEW CLEANING/ RECOVER PARAMS FUNCTION
recover_interval <- function(sims, parameters, distribution){
  clean_df <- sims%>%
    mutate(event = if_else(t>60, 0, event),
           t = if_else(t>60, 60, t))%>%
    mutate(left_interval = case_when(t <= 1 & event == 1 ~ t,
                                     t > 1 & t <= 24 & event == 1 ~ floor(t),
                                     t > 24 & t <= 36 & event == 1 ~ 24,
                                     t > 36 & t <= 48 & event == 1 ~ 36,
                                     t > 48 & event == 1~ 48, .default = 1000),
           right_interval = case_when(t <= 1& event == 1 ~ t,
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
    mutate(age_at_begin = max_age- period_length)
  
  return (clean_df)
  
  res <- surv_synthetic(df = clean_df, survey = FALSE, individual = "id", p = "period", a_pi = "age_at_begin", l_p = "period_length", I_i = "interval_indicator",
                        A_i = "across_boundary", t_i = "right_censor_age", t_0i = "left_interval", t_1i = "right_interval", dist = "lognormal", numerical_grad = TRUE)
  
  # return (res)
}



#WHAT IT SHOULD LOOK LIKE
filename <- file.choose()
example <- load("../testdf.rds")
view(example)





#TESTING RECOVER_INTERVAL()

mus <- c(15,16,17,18,19)
sigmas <- exp(c(2.1,2.2,2.3,2.2,2.1))
period_length <- c(12,12,12,12,12)

matrix_lnorm <- cbind(mus, sigmas, period_length)

sim_5_ln <- general_sim(num_child = 100, param_matrix = matrix_lnorm, distribution = "lognormal")
res <- recover_interval(sim_5_ln, matrix_lnorm, "lognormal")
view(res)

#-------------------------------------------------

#sim_age_heap TEST

shapes_k <- c(-1,-1)
scales_lam <- c(12,15)
period_length <- c(60,60)
matrix_params <- cbind(exp(scales_lam), exp(shapes_k), period_length)

example_sim <- general_sim(num_child = 1000, 
                           param_matrix = matrix_params, 
                           distribution = "weibull")

ages <- c(6, 12)
proportions <- c(.5,.5)
ranges <- rbind(c(2,3), c(4,5))
age_heap_df <- sim_age_heap(ages_at_heaping = ages, proportion_heap = proportions, range_heap = ranges, sim_data = example_sim)

#-------------------------------------------------

#five periods TEST

#Weibull test

shapes <- c(-1.2, -1.1,-1,-.9, -.8)
scales <- c(8,9,10,11, 12)
period_length <- c(12,12,12,12,12)
matrix_dist <- cbind(exp(scales), exp(shapes), period_length)

sims_5 <- general_sim(num_child = 100, param_matrix = matrix_dist, distribution = "weibull")
params_recovered_df <- recover_params(sims = sims_5, parameters = matrix_dist, distribution = "weibull")
view(params_recovered_df)

#Lognormal test

mus <- c(15,16,17,18,19)
sigmas <- exp(c(2.1,2.2,2.3,2.2,2.1))
period_length <- c(60,60,60,60,60)

matrix_lnorm <- cbind(mus, sigmas, period_length)

sim_5_ln <- general_sim(num_child = 1000, param_matrix = matrix_lnorm, distribution = "lognormal")

params_5_ln_recover <- recover_params(sims = sim_5_ln, parameters = matrix_lnorm, distribution = "lognormal")
#--------------------------------------------------

#two periods TEST

#Weibull Test

shapes_k <- c(-1,-1)
scales_lam <- c(12,15)
period_length <- c(60,60)
matrix_params <- cbind(exp(scales_lam), exp(shapes_k), period_length)

sims_2 <- general_sim(num_child = 1000, param_matrix = matrix_params, distribution = "weibull")
params_2_recover <- recover_params(sims = sims_2, parameters = matrix_params, distribution = "weibull")
view(params_2_recover)
#lognormal test

mus <- c(11.44, 11.44)
sigmas <- exp(c(1.8,1.8))
period_length <- c(60,60)

matrix_lnorm <- cbind(mus, sigmas, period_length)

sim_2_ln <- general_sim(num_child = 1000, param_matrix = matrix_lnorm, distribution = "lognormal")
view(sim_2_ln)
params_2_ln_recover <- recover_params(sims = sim_2_ln, parameters = matrix_lnorm, distribution = "lognormal")
#--------------------------------------------------

plnorm()
#one period TEST

#Weibull Test

shape_k <- -1
scale_lam <- 12
period_length <- 60
matrix_weib <- cbind(exp(scale_lam), exp(shape_k), period_length)

sim_1 <- general_sim(num_child = 10000, param_matrix = matrix_weib, distribution = "weibull")
params_1_recover <- recover_params(sims = sim_1, parameters = matrix_weib, distribution = "weibull")

#Lognormal test

mu <- 11.44
sigma <- exp(1.8)
period_length <- 60
matrix_lnorm <- cbind(mu, sigma, period_length)

sim_1_ln <- general_sim(num_child = 10000, param_matrix = matrix_lnorm, distribution = "lognormal")

params_1_ln_recover <- recover_params(sims = sim_1_ln, parameters = matrix_lnorm, distribution = "lognormal")
