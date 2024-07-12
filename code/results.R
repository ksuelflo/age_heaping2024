source("sim_functions.R")

#vector of folders: Each folder represents a different sim setting.
summary_folders <- list.files(path = "../summaries")

#making big container to store all the `container` into.
big_container <- vector("list", length = length(summary_folders))

#Looping through each folder (each sim setting)
for (i in seq_along(summary_folders)){
  
  #getting list of summaries for the ith sim setting
  summaries <- list.files(path = str_c("../summaries", summary_folders[i]))
  
  #Create a container for each sim setting
  container <- vector("list", length = length(summaries))
  
  #Looping through each summary (all the same sim setting, different seeds).
  for (j in seq_along(summaries)){
    
    summary <- readRDS(str_c("../summaries", summary_folders[i], summaries[j]))%>%
      mutate(seed = str_extract(summaries[j], ))
    container[[j]] <- summary
    
  }
  
  #Adding the list of dfs for the ith sim setting to the master container.
  big_container[[i]] <- container
  
}

#Now we have read in all RDS files and put them in one big list, with sublists for each sim setting.

#looping through each sim setting and combining each seed into one big df. 
for (k in seq_along(big_container)){
  
  big_container[[k]] <- bind_rows(big_container[[k]])
  
}



