library(readstata13)
library(tidyverse)
library(pssst)
library(knitr)
library(SUMMER)

#Randomly generated numbers, to then use in DHS alphabetical list to select countries
random_nums <- sample(1:42, 10, replace = F)
random_countries <- c("Burkina Faso, Cameroon", "Chad", "Eritrea", "Ghana", 
                      "Madagascar", "Malawi", "Mauritania", "Nigeria", "Rwanda")



folders <- list.files(path = "/Users/kylesuelflow/Macalester-Stuff/Research-Taylor/data")
years <- as.numeric(str_extract(folders, "\\d+"))
dfs <- vector(mode = "list", length = length(folders))
countries <- str_extract(folders, ".*(?=_(\\d)*)")
intervals <- vector("list", length(years))

for (i in seq_along(years)){
  intervals[[i]] <- c(years[i] - 5, years[i])
}

intervals[[17]] <- c(2006, 2013)
intervals[[18]] <- c(2011, 2018)

mus <- rep(0,length(folders))
sigmas <- rep(0, length(folders))

for (i in 18:length(folders)){
  print(str_c(countries[i], " ", years[i], " begins."))
  path <- str_subset(string = list.files(path = str_c("../data/", folders[i])), pattern = ".DTA") 
  raw_data <- read.dta13(str_c("../data/", folders[i], "/", path))
  cleaned_data <- format_dhs(df = raw_data, survey_year = years[i], period_boundaries = intervals[[i]])
  res <- surv_synthetic(cleaned_data, household = "household", strata = "strata", cluster = "cluster", weights = "weights", dist = "lognormal", numerical_grad = TRUE)
  dfs[[i]] <- res$result
  print(str_c(countries[i], " ", years[i], " was successful."))
}

print(str_c(countries[17], " ", years[17], " begins."))
path <- str_subset(string = list.files(path = str_c("../data/", folders[17])), pattern = ".DTA") 
raw_data <- read.dta13(str_c("../data/", folders[17], "/", path))

raw_data%>%
  count(b2)
view(raw_data)

for (i in 1:16){
  mus[i] <- dfs[[i]]$mu_mean
  sigmas[i] <- dfs[[i]]$log_sigma_mean
}
for (i in 18:21){
  mus[i] <- dfs[[i]]$mu_mean
  sigmas[i] <- dfs[[i]]$log_sigma_mean
}

view(raw_data)

params_DHS <- data.frame(Country = countries, Year = years, mu = mus, sigma = sigmas)%>%
  arrange(Year)
view(params_DHS[-10,])
kable(params_DHS, "latex")
