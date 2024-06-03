# Project: Age-Heaping Simulation
# Authors: Kyle Suelflow & Taylor Okonek
# Date: Summer 2024


## Load necessary packages up front
library(survival)
library(tidyverse)
library(flexsurv)


# Functions ---------------------------------------------------------------

#function to fill the vector of probabilities of death. 

#' Add roxygen documentation here
helper_fill_p <- function(distribution, param_matrix, num_child, months, max_age, age_at_begin, ids){
  
  # p is the vector we will populate. periods and len_period are helpful variables to have for later. 
  # Note from Taylor - what exactly is len_period? Could use more documentation here
  periods <- nrow(param_matrix)
  len_period <- num_child*months*periods
  curr_period <- 1
  p <- rep(0, periods*len_period)
  
  for (i in seq_along(age_at_begin)){
    
    if (curr_period == periods){
      curr_period <- 1
      next
    }
    
    if (distribution == "weibull"){
      #Added this if statement: My thought is, if the max_age is greater than 60, then we should just cap it at 60. Otherwise, 
      #the probability `p` is higher than it should be, as it includes the chance of dying from 60 to max_age. Since we
      #are censoring at 60, this does not make sense to have. HOWEVER, this only exacerbates the difference between the recovered parameters
      #and the actual parameters for periods after period 1. 
      
      #Commented out in order to remove censoring.
      if (max_age[i] > 60){
        max_age[i] <- 60
      }
      
      p[i] <- (pweibull(q = max_age[i], shape = param_matrix[curr_period,2], scale = param_matrix[curr_period,1]) - pweibull(q = age_at_begin[i], shape = param_matrix[curr_period,2], scale = param_matrix[curr_period,1]))/
        (1 - pweibull(q = age_at_begin[i], shape = param_matrix[curr_period,2], scale = param_matrix[curr_period,1]))
      
    }
    
    else if (distribution == "lognormal"){
      #Use plnorm
    }
    
    else{
      #Use pgengamma
    }
    curr_period <- curr_period + 1
    
    # In the if statement below this one, we check if the next observation has an age_at_begin of 0, 
    # essentially checking if the next observation is a new child. This if statement is included so there 
    # is not an index out of bounds error on the final observation.
    # if (i == length(age_at_begin)){
    #   break
    # }
    # 
    # # Checks if the next observation is a new child. If it is, reset the curr_period to 1.
    # if (ids[i+1] != ids[i]){
    #   curr_period <- 1
    # }
    
  }
  return (p)

}

#' Add roxygen documentation here
general_sim <- function(num_child, param_matrix, distribution = "weibull", months){
  
  #Setting up data frame. The first five plus `age_at_begin` do not iteratively change. `t` and `event` get populated using for loop below.
  
  periods <- nrow(param_matrix)
  ids <- rep(1:(num_child*periods*months), each = periods)
  birthdates <- rep(1:(months*periods), each = num_child*periods)
  period <- rep(1:periods, times = num_child*periods*months)
  max_age <- (months*period - birthdates) + 1
  max_age[max_age < 0] <- 0
  age_at_begin <- c(0,max_age)
  age_at_begin[seq(periods+1, length(age_at_begin), periods)] <- 0
  age_at_begin <- age_at_begin[-length(age_at_begin)]
  
  #FOR NO CENSORING: This sets the max age for the last time period to be infinite (100000000).
  max_age[period == periods] <- Inf
  
  # Note from Taylor - put in a dataframe, easier to debug
  sim_df <- data.frame(id = ids,
             birthdate = birthdates,
             period = period,
             max_age = max_age,
             age_at_begin = age_at_begin)

  # Filling p based on distribution.
  p <-helper_fill_p(distribution = distribution, param_matrix = param_matrix, num_child = num_child, months = months, max_age = max_age, age_at_begin = age_at_begin, ids = ids)
  
  # Note from Taylor - add to dataframe, easier to debug
  sim_df$p <- p
  
  # returning this for debugging at this point
  # return(sim_df)
  
  # This column (age_at_begin) is useful for the for loop: Essentially this allows me to "group_by" 
  # by checking if the age is 0. If it is, then I know it is a new child and to reset the boolean `already_died`
  
  # Note from Taylor - should this be all zeros? Could you just do this the same way as you 
  # define t and event?

  
  # these two columns will be iteratively populated in the for loop. T is time of event within period, 
  # event is a binary outcome, either 1 or 0.
  
  t <- rep(0, num_child*periods*months*periods)
  event <- rep(0, num_child*periods*months*periods)

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
      if (curr_period == periods){
        died <- 1
      }
      
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
          time <- 0
        }
        
        else{
          time <- 0
        }
        
        #We continue randomly sampling deaths until the time of death falls within the age at beginning and max_age of the time period. 
        while (time >= max_age[i] | time < age_at_begin[i]){
          
          if (distribution == "weibull"){
            time <- rweibull(n = 1, scale = param_matrix[curr_period, 1], shape = param_matrix[curr_period,2])
          }
          
          else if (distribution == "lognormal"){
            time <- 0
          }
          
          else{
            time <- 0
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

#' Add roxygen documentation here
weibull_params <- function(low, high, periods){
  if (periods == 1){
    return (exp(mean(c(low,high))))
  }
  else if (low == high){
    return (rep(exp(low), periods))
  }
  else{
    return (exp(seq(low, high, (high-low)/(periods -1 ))))
  }
}

#' @description
#' This function simulates the phenomenon of age heaping in child mortality. It takes in a data frame of child mortality data, generated
#' from the `general_sim` function. It then randomly chooses deaths around specific ages of death to "heap" at a round number. This is
#' usually at 3, 6, 12, or 24 months. The interval surrounding the age, the age itself, and the probability to be "chosen" to be included in
#' the heap are all parameters, `range_heap`, `ages_at_heaping`, `proportion_heap` respectively.
#' 
#' @param ages_at_heaping a vector of ages, in months, that represent ages where deaths will be artificially clumped at. 
#' @param proportion_heap a vector of proportions. Each index is the probability of being chosen to 
#' @param range_heap
#' @param sim_data 
#' 
#' @returns description
sim_age_heap <- function(ages_at_heaping, proportion_heap, range_heap, sim_data){
  
  times <- sim_data%>%pull(t)
  events <- sim_data%>%pull(event)
  intervals <- vector("list", length = length(ages_at_heaping))
  indicators <- rep(0, length(times))
  
  for (i in seq_along(ages_at_heaping)){
    intervals[[i]] <- c(ages_at_heaping[i] - range_heap[i,1], ages_at_heaping[i] + range_heap[i,2])
  }
  
  for (i in seq_along(times)){
    
    if (events[i] == 0){
      next
    }
    
    for (j in seq_along(ages_at_heaping)){
      if (intervals[[j]][1] <= times[i] & times[i] <= intervals[[j]][2]){
        indicators[i] <- rbinom(1,1,proportion_heap[j])
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

example_sim <- general_sim(num_child = 1000, # Note from Taylor - it'll take a while, but I want you to try this with an even bigger number. Add another zero, and see if you get closer to the "truth"
                           param_matrix = matrix_params, 
                           distribution = "weibull", 
                           months = 60)

ages <- c(6, 12)
proportions <- c(.5,.5)
ranges <- rbind(c(2,3), c(4,5))
age_heap_df <- sim_age_heap(ages_at_heaping = ages, proportion_heap = proportions, range_heap = ranges, sim_data = example_sim)

recover_params <- function(sims, parameters, distribution){
  periods <- nrow(parameters)
  period <- 1:periods
  rec_shapes <- rep(0, periods)
  rec_scales <- rep(0, periods)
  
  for (i in 1:nrow(parameters)){
    
    curr_period_data <- sims%>%
      filter(period == i,
             !is.na(event),
             t != 0)
    
    if (i == 1){
      res <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = curr_period_data, dist = "weibull")
      rec_shapes[i] <- log(1/res$scale)
      rec_scales[i] <- coef(res)
    }
    
    else{
      res <- flexsurvreg(formula = Surv(time = age_at_begin, time2 = t, event = event) ~ 1, data = curr_period_data, dist = "weibull")
      rec_shapes[i] <- res$coefficients[1]
      rec_scales[i] <- res$coefficients[2]
    }
    
  }
  
  return (data.frame(period = period, rec_shape = rec_shapes, shape = log(parameters[,2]), rec_scale = rec_scales, scale = log(parameters[,1])))

}

#five periods

shapes <- c(-1.2, -1.1,-1,-.9, -.8)
scales <- c(8,9,10,11, 12)
matrix_dist <- cbind(exp(scales), exp(shapes))

sims_5 <- general_sim(num_child = 1000, param_matrix = matrix_dist, distribution = "weibull", months = 60)

params_recovered_df <- recover_params(sims = sims_5, parameters = matrix_dist, distribution == "weibull")

view(params_recovered_df)

#two periods

shapes_k <- c(-1,-1)
scales_lam <- c(12,15)
matrix_params <- cbind(exp(scales_lam), exp(shapes_k))

sims_2 <- general_sim(num_child = 1000, param_matrix = matrix_params, distribution = "weibull", months = 60)

params_2_recover <- recover_params(sims = sims_2, parameters = matrix_params, distribution == "weibull")
view(params_2_recover)

# Code Testing - One time period ------------------------------------------

# Building out parameter matrix. Eventually this needs to be a function.
shapes_k <- weibull_params(-1.5, -.5,1)
scales_lam <- weibull_params(15, 17, 1)

# matrix_params is a n x m matrix, where n is the number of periods, and m is the 
# number of parameters in the distribution. 
matrix_params <- cbind(scales_lam, shapes_k)

# Simulate data
sims <- general_sim(num_child = 100000, # Note from Taylor - it'll take a while, but I want you to try this with an even bigger number. Add another zero, and see if you get closer to the "truth"
                    param_matrix = matrix_params, 
                    distribution = "weibull", 
                    months = 60)

# Testing period 1 parameters.
res_1 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims, dist = "weibull")
shape_1 <- log(1/res_1$scale) #Should be exp(-1) or .368 || LOGGED should be -1
scale_1 <- log(exp(coef(res_1))) #Should be exp(16) or 8886111 || LOGGED should be 16

# Code Testing - Two time periods -----------------------------------------

shapes_k_2 <- weibull_params(low = -1.5, high = -.5, periods = 2)
scales_lam_2 <- weibull_params(low = 15, high = 17, periods = 2)
matrix_params_2 <- cbind(scales_lam_2, shapes_k_2)

sims_2 <- general_sim(num_child = 10000, 
                      param_matrix = matrix_params_2, 
                      distribution = "weibull", 
                      months = 60)

sims_2_clean <- sims_2%>%
  mutate(keep_row = if_else((period == 2 & birthdate < 121) |(period == 1 & birthdate < 61), TRUE, FALSE))%>%
  filter(keep_row)%>%
  na.omit(event)%>%
  filter(age_at_begin < 60)%>%
  mutate(event = if_else(t >= 60, 0, event),
         t = if_else(t>60, 60, t))%>%
  select(-keep_row)

sims_2_1 <- sims_2_clean%>%
  filter(period == 1)

sims_2_2 <- sims_2_clean%>%
  filter(period == 2)

res_1 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims_2_1, dist = "weibull")
shape_1 <- log(1/res_1$scale) #Should be exp(-1) or .368 || LOGGED should be -1.5
scale_1 <- log(exp(coef(res_1))) #Should be exp(16) or 8886111 || LOGGED should be 15

res_2 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims_2_2, dist = "weibull")
shape_2 <- log(1/res_2$scale) #Should be exp(-1) or .368 || LOGGED should be -.5
scale_2 <- log(exp(coef(res_2))) #Should be exp(16) or 8886111 || LOGGED should be 17

res_flex <- flexsurvreg(formula = Surv(time = age_at_begin, time2 = t, event = event, type = "counting") ~ 1, data = sims_2_2, dist = "weibull")
res_flex$coefficients
# Code Testing - Five time periods ----------------------------------------

#Building out parameter matrix. Eventually this needs to be a function.
shapes_k_5 <- weibull_params(low = -2, high = -1, periods = 5)
scales_lam_5 <- weibull_params(low = 14, high = 16, periods = 5)
#matrix_params is a nxm matrix, where n is the number of periods, and m is the number of parameters in the distribution. 
matrix_params_5 <- cbind(scales_lam_5, shapes_k_5)

sims_5 <- general_sim(num_child = 1000, param_matrix = matrix_params_5, distribution = "weibull", months = 60)

#Removing all observations where the child is already dead (event == NA). Also removes observations where the child was born in a later period. (For example, if the child was born in period 5, the observations in periods 1 through 4 are removed.)
sims_clean <- sims_5%>%
  mutate(keep_row = if_else(period == 5 |
                              (period == 4 & birthdate < 241) |
                              (period == 3 & birthdate < 181) |
                              (period == 2 & birthdate < 121) |
                              (period == 1 & birthdate < 61), TRUE, FALSE))%>%
  filter(keep_row)%>%
  na.omit(event)%>%
  #Code to omit children who are already past the 5 year mark when a time period began. Also, censors both t and event at 60 months.
  filter(age_at_begin < 60)%>%
  mutate(event = if_else(t >= 60, 0, event),
         t = if_else(t>60, 60, t))%>%
  select(-keep_row)

#Data frames for each period, in order to test parameters.

sims_clean_1 <- sims_clean%>%
  filter(period == 1)

sims_clean_2 <- sims_clean%>%
  filter(period == 2)

sims_clean_3 <- sims_clean%>%
  filter(period == 3)

sims_clean_4 <- sims_clean%>%
  filter(period == 4)

sims_clean_5 <- sims_clean%>%
  filter(period == 5)

#Testing period 1 parameters.
res_1 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims_clean_1, dist = "weibull")
shape_1 <- log(1/res_1$scale) #Should be -2
scale_1 <- log(exp(coef(res_1))) #Should 14

#Testing period 2 parameters.
res_2 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims_clean_2, dist = "weibull")
shape_2 <- log(1/res_2$scale) #Should be -1.75
scale_2 <- log(exp(coef(res_2))) #Should be 14.5

#Testing period 3 parameters.
res_3 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims_clean_3, dist = "weibull")
shape_3 <- log(1/res_3$scale) #Should be -1.5
scale_3 <- log(exp(coef(res_3))) #Should be 15

#Testing period 4 parameters.
res_4 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims_clean_4, dist = "weibull")
shape_4 <- log(1/res_4$scale) #Should be -1.25
scale_4 <- log(exp(coef(res_4))) #Should be 15.5

#Testing period 5 parameters.
res_5 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims_clean_5, dist = "weibull")
shape_5 <- log(1/res_5$scale) #Should be -1
scale_5 <- log(exp(coef(res_5))) #Should be 16

res <- flexsurvreg(formula = Surv(time = age_at_begin, time2 = t, event = event) ~ 1, data = sims_clean_2, dist = "weibull")
res$coefficients

res5 <- survreg(formula = Surv(time = t, event = event, type = "right") ~ 1, data = sims_clean_5, dist = "weibull")
