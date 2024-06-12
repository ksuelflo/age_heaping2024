library(readstata13)
library(tidyverse)
library(pssst)
library(knitr)

folders <- list.files(path = "/Users/kylesuelflow/Macalester-Stuff/Research-Taylor/data")
paths <- str_c("data/", folders, "/", str_subset(string = list.files(path = str_c(getwd(), "../data/", folders)), pattern = "DTA"))
years <- as.numeric(str_extract(folders, "\\d+"))
dfs <- vector(mode = "list", length = length(folders))
countries <- str_extract(folders, ".*(?=_(\\d)*)")
intervals <- list(c(2005,2010),c(2016,2021),c(1998,2003),c(2003,2008),c(2003,2008),c(2009,2014),c(2017,2022),c(2005,2010),c(2010,2015))
mus <- rep(0,length(folders))
sigmas <- rep(0, length(folders))



for (i in seq_along(paths)){
  print(str_c(countries[i], " ", years[i], " begins."))
  raw_data <- read.dta13(paths[i])
  cleaned_data <- format_dhs(df = raw_data, survey_year = years[i], period_boundaries = intervals[[i]])
  res <- surv_synthetic(cleaned_data, household = "household", strata = "strata", cluster = "cluster", weights = "weights", dist = "lognormal", numerical_grad = TRUE)
  dfs[[i]] <- res$result
  print(str_c(countries[i], " ", years[i], " was successful."))
}

raw_data <- read.dta13(paths[1])
cleaned_data <- format_dhs(df = raw_data, survey_year = years[1], period_boundaries = intervals[[1]])
res <- surv_synthetic(cleaned_data, household = "household", strata = "strata", cluster = "cluster", weights = "weights", dist = "lognormal", numerical_grad = TRUE)
dfs[[1]] <- res$result

raw <- read.dta13("../data/Burkina_Faso_2010/BFBR62FL.DTA")
clean <- format_dhs(df = raw, survey_year = 2010, period_boundaries = c(2005,2010))
res <- surv_synthetic(clean, household = "household", strata = "strata", cluster = "cluster", weights = "weights", dist = "lognormal", numerical_grad = TRUE)

res$result

for (i in seq_along(dfs)){
  mus[i] <- dfs[[i]]$mu_mean
  sigmas[i] <- dfs[[i]]$log_sigma_mean
}



params_DHS <- data.frame(Country = countries, Year = years, mu = mus, sigma = sigmas)
view(params_DHS)
kable(params_DHS, "latex")
