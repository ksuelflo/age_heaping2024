source("sim_functions.R")
library(ggplot2)

#Build folders

lnorm_mean_vec <- c(12,15,20)
lnorm_sd_vec <- c(1.9,2.1,2.3)
sample_size_vec_old <- c(100,500)
age_1_vec <- 12
range_1_vec <- c("9,21", "6,18", "8,24")
period_length_vec <- c(60)
periods_vec <- 5
proportion_1_vec <- c(0,.10,.20,.50)

#Every possible combination of sim settings found using expand.grid. Even ones that aren't plausible.

clean_params_old <- expand.grid(lnorm_mean = lnorm_mean_vec, lnorm_sd = lnorm_sd_vec, sample_size = sample_size_vec_old, age_1 = age_1_vec, 
                            range_1 = range_1_vec, period_length = period_length_vec, 
                            periods = periods_vec, proportion_1 = proportion_1_vec, stringsAsFactors = FALSE)

#Filtering out combos to only include plausible combinations.
clean_params_old <- clean_params_old%>%
  mutate(double = lnorm_sd + lnorm_mean)

#Filtering out combos to only include plausible combinations.
clean_params_old <- clean_params_old%>%
  filter(double %in% c(13.9, 17.3, 17.1, 22.3))%>%
  dplyr::select(-double)

#SAME THING, just with samplesize= 300,1000

sample_size_vec_new <- c(300,1000)

clean_params_new <- expand.grid(lnorm_mean = lnorm_mean_vec, lnorm_sd = lnorm_sd_vec, sample_size = sample_size_vec_new, age_1 = age_1_vec, 
                                range_1 = range_1_vec, period_length = period_length_vec, 
                                periods = periods_vec, proportion_1 = proportion_1_vec, stringsAsFactors = FALSE)

#Filtering out combos to only include plausible combinations.
clean_params_new <- clean_params_new%>%
  mutate(double = lnorm_sd + lnorm_mean)

#Filtering out combos to only include plausible combinations.
clean_params_new <- clean_params_new%>%
  filter(double %in% c(13.9, 17.3, 17.1, 22.3))%>%
  dplyr::select(-double)

#Combining clean_params_old with clean_params_new to get one 196 length sim setting df

clean_params <- clean_params_old%>%
  rbind(clean_params_new)

#Creating directories in "../summaries_100"

for (i in 1:nrow(clean_params)){
  foldername <- str_c("sim_ROW=", i, "_lnormmean=", clean_params$lnorm_mean[i], "_lnormsd=", clean_params$lnorm_sd[i],
                      "_samplesize=", clean_params$sample_size[i], "_age=", clean_params$age_1[i], "_range=", clean_params$range_1[i],
                      "_periodlength=", clean_params$period_length[i], "_periods=", clean_params$periods[i],
                      "_proportion=", clean_params$proportion_1[i])
  dir.create(path = str_c("../summaries_100/", foldername))
}


#------------------------------------------------------------

#100 dfs per sim setting Analysis

#Reading them in 

#vector of folders: Each folder represents a different sim setting.
summary_folders <- list.files(path = "../summaries_100")

#making big container to store all the `container` into.
big_container <- vector("list", length = length(summary_folders))

#Looping through each folder (each sim setting)
for (i in seq_along(summary_folders)){
  
  print(i)
  #getting list of summaries for the ith sim setting
  summaries <- list.files(path = str_c("../summaries_100/", summary_folders[i]))
  
  #Create a container for each sim setting
  container <- vector("list", length = length(summaries))
  
  #Looping through each summary (all the same sim setting, different seeds).
  for (j in seq_along(summaries)){
    
    summary <- readRDS(str_c("../summaries_100/", summary_folders[i], "/",summaries[j]))%>%
      mutate(seed = str_extract(summaries[j], pattern = "(?<=seed=)\\d*"),
             lnorm_mean = as.numeric(str_extract(summaries[j], pattern = "(?<=lnormmean=)\\d*")),
             lnorm_sd = as.numeric(str_extract(summaries[j], pattern = "(?<=lnormsd=)\\d*"))/10,
             sample_size = as.numeric(str_extract(summaries[j], pattern = "(?<=samplesize=)\\d*")),
             low_range = as.numeric(str_extract(summaries[j], pattern = "(?<=range=)\\d")),
             high_range = as.numeric(str_extract(summaries[j], pattern = "(?<=range=\\d)\\d*")),
             proportion = as.numeric(str_extract(summaries[j], pattern = "(?<=proportion=)\\d*"))/10,
             sim_setting = i)
    summary$model[21:25] <- "logquad_adj"
    container[[j]] <- summary
  }
  #Combining all the seeds into one df
  container <- bind_rows(container)
  
  #Adding the list of dfs for the ith sim setting to the master container.
  big_container[[i]] <- container
  
}

saveRDS(big_container, "../data/summary_192_100dfs_list.rds")




