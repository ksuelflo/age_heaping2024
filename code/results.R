source("sim_functions.R")

#Getting all generated data frames' file paths. Setting up directory structure
# all_folders <- list.files(path = "../models")
# summaries_folders <- str_replace(all_folders, "res", "summary")

# for (i in seq_along(all_folders)){
#   dir.create(path = str_c("../summaries/", summaries_folders[i]))
# }

#vector of folders: Each folder represents a different sim setting.
summary_folders <- list.files(path = "../summaries2")

#making big container to store all the `container` into.
big_container <- vector("list", length = length(summary_folders))

#Looping through each folder (each sim setting)
for (i in seq_along(summary_folders)){
  print(i)
  #getting list of summaries for the ith sim setting
  summaries <- list.files(path = str_c("../summaries2/", summary_folders[i]))
  
  #Create a container for each sim setting
  container <- vector("list", length = length(summaries))
  
  #Looping through each summary (all the same sim setting, different seeds).
  for (j in seq_along(summaries)){
    
    summary <- readRDS(str_c("../summaries2/", summary_folders[i], "/",summaries[j]))%>%
      mutate(seed = str_extract(summaries[j], pattern = "(?<=seed=)\\d*"),
             lnorm_mean = as.numeric(str_extract(summaries[j], pattern = "(?<=lnormmean=)\\d*")),
             lnorm_sd = as.numeric(str_extract(summaries[j], pattern = "(?<=lnormsd=)\\d*"))/10,
             sample_size = as.numeric(str_extract(summaries[j], pattern = "(?<=samplesize=)\\d*")),
             low_range = as.numeric(str_extract(summaries[j], pattern = "(?<=range=)\\d")),
             high_range = as.numeric(str_extract(summaries[j], pattern = "(?<=range=\\d)\\d*")),
             proportion = as.numeric(str_extract(summaries[j], pattern = "(?<=proportion=)\\d*"))/10)
    summary$model[21:25] <- "logquad_adj"
    container[[j]] <- summary
  }
  #Combining all the seeds into one df
  container <- bind_rows(container)
  
  #Adding the list of dfs for the ith sim setting to the master container.
  big_container[[i]] <- container
  
}

saveRDS(big_container, "../Results/in_progress_summary_list2.rds")

