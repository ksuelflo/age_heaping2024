library(readstata13)
library(tidyverse)
library(pssst)
library(knitr)
library(SUMMER)
library(survey)

folders <- list.files(path = "/Users/kylesuelflow/Macalester-Stuff/Research-Taylor/data")
paths <- str_c("../../data/", folders, "/", str_subset(string = list.files(path = str_c("../../data/", folders)), pattern = "DTA"))
dfs <- vector(mode = "list", length = length(folders))
years <- as.numeric(str_extract(folders, "\\d+"))
countries <- str_extract(folders, ".*(?=_(\\d)*)")
intervals <- list(c(2005,2010),c(2016,2021),c(1998,2003),c(2003,2008),c(2003,2008),c(2009,2014),c(2017,2022),c(2005,2010),c(2010,2015))

for (i in seq_along(folders)){
  
  #Reading in data
  births <- read.dta13(paths[i], generate.factors = TRUE)
  
  #Formatting using getBirths, also does other cool stuff.
  dat <- getBirths(data = births, strata = c("v023"), surveyyear = years[i], year.cut = intervals[[i]])
  
  #Renaming columns to be more clear.
  dat <- dat[, c("v001", "v002", "v024", "time", "age", "v005", "strata", "died")]
  colnames(dat) <- c("clustid", "id", "region", "time", "age", "weights", "strata",
                     "died")
  dfs[[i]] <- dat
  
}

NMRs <- rep(0,length(dfs))
IMRs <- rep(0,length(dfs))
U5MRs <- rep(0,length(dfs))

#Formatting `intervals` do be a vector of strings like "2005-2010", instead of a list of numeric vectors.
formatted_intervals <- sapply(intervals, function(x) str_c(x, collapse = "-"))


for (i in seq_along(dfs)){
  
  #Incoporating survey design.
  design_births <- svydesign(id = ~clustid, 
                             strata = ~strata,
                             weights = ~weights, 
                             data = dfs[[i]] %>% mutate(weights = weights/1e06))
  mod <- svyglm(died ~ -1 + age, design = design_births, family = quasibinomial(link = "logit"))
  betas <- mod$coefficients
  probs <- expit(betas)
  NMRs[i] <- probs[1]
  IMRs[i] <- 1- prod(1- probs[1:2])
  U5MRs[i] <- 1 - prod(1- probs)

}

results <- data.frame(country = countries, survey_year = years, period = formatted_intervals, NMR = NMRs, IMR = IMRs, U5MR = U5MRs)

results_longer <- pivot_longer(results, cols = c(NMR, IMR, U5MR), names_to = "mortality", values_to = "y")%>%
  mutate(age = rep(c(1,12,60), length(dfs)))

#Getting lnorm distribution.
y <- dlnorm(ages, 15, exp(2))
y_update <- rep(0, length(y))

for (i in seq_along(y)){
  #Applying discrete hazards method to lnorm output. 
  y_update[i] <- 1 - prod(1- y[1:i])
  
}

ages <- 0:60

ggplot()+
  geom_line(aes(x = ages, y = y_update))+
  geom_line(data = results_longer, aes(x = age, y = y, group = interaction(country, survey_year), 
                                       color = interaction(country, survey_year)))+
  labs(x = "age (months)", y = "mortality rate", color = "Country Year")+
  theme_classic()
  

