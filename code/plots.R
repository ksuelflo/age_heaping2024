library(readstata13)
library(tidyverse)
library(pssst)
library(knitr)
library(SUMMER)
library(survey)
library(ggnewscale)

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
year2 <- rep(0, length(dfs))
year3 <- rep(0, length(dfs))
year4 <- rep(0, length(dfs))
U5MRs <- rep(0,length(dfs))
zeros <- rep(0,length(dfs))

#Formatting `intervals` do be a vector of strings like "2005-2010", instead of a list of numeric vectors.
formatted_intervals <- sapply(intervals, function(x) str_c(x, collapse = "-"))


for (i in seq_along(dfs)){
  
  #Incorporating survey design.
  design_births <- svydesign(id = ~clustid, 
                             strata = ~strata,
                             weights = ~weights, 
                             data = dfs[[i]] %>% mutate(weights = weights/1e06))
  mod <- svyglm(died ~ -1 + age, design = design_births, family = quasibinomial(link = "logit"))
  betas <- mod$coefficients
  probs <- expit(betas)
  
  probs <- expit(betas)
  
  #Because each interval may be different lengths, we need to keep that in mind. (0,1), for example, should not get
  #the same weighting as (12,24).
  interval_lengths <- c(1,11,12,12,12,12)
  
  NMRs[i] <- probs[1]
  IMRs[i] <- (1- prod((1- probs[1:2])^interval_lengths[1:2]))
  U5MRs[i] <- (1 - prod((1 - probs)^interval_lengths))
  year2[i] <- (1- prod((1- probs[1:3])^interval_lengths[1:3]))
  year3[i] <- (1- prod((1- probs[1:4])^interval_lengths[1:4]))
  year4[i] <- (1- prod((1- probs[1:5])^interval_lengths[1:5]))
  
}

results <- data.frame(country = countries, survey_year = years, period = formatted_intervals, NMR = NMRs, IMR = IMRs, year2 = year2,
                      year3 = year3, year4 = year4, U5MR = U5MRs, zeros = zeros)

results_longer <- pivot_longer(results, cols = c(zeros,NMR, IMR, year2, year3, year4, U5MR), names_to = "mortality", values_to = "y")%>%
  mutate(age = rep(c(0, 1, 12, 24, 36, 48, 60), length(dfs)))
  
view(results_longer)

#Getting lnorm distribution.

ages <- 0:60
y <- list(dlnorm(ages, 15, exp(2)), dlnorm(ages, 20, exp(2.4)), dlnorm(ages, 15, exp(2.2)), dlnorm(ages, 12, exp(2.2)))
y_update <- list(rep(0, length(y[[1]])), rep(0, length(y[[1]])),rep(0, length(y[[1]])),rep(0, length(y[[1]])))

for (i in seq_along(y)){
  
  for (j in seq_along(y[[i]])){
    #Applying discrete hazards method to lnorm output. 
    y_update[[i]][j] <- 1 - prod(1- y[[i]][1:j])
  }
}

lnorm_data <- data.frame(y = list_c(y_update), age = rep(ages, 4), type = rep(c("lnorm(15, exp(2)", "lnorm(20, exp(2.4)", 
                                                                                "lnorm(15, exp(2.2)", "lnorm(12, exp(2.2)"), each = 61))
view(lnorm_data)
p <- ggplot()+
  geom_line(data = lnorm_data, aes(x = age, y = y, linetype = type))+
  geom_line(data = results_longer, aes(x = age, y = y, group = interaction(country, survey_year), 
                                       color = interaction(country, survey_year)))+
  labs(x = "age (months)", y = "mortality rate", color = "Country Year", linetype = "lognormal curves", caption = "Lognormal curve uses mean of 15, sd of exp(2.2).
       The resulting values are updated using the discrete hazards approach.")+
  scale_linewidth_manual(values = c(.5, .5, .5, .5))+
  theme_classic()
  
p
#Save plot when finished:

# ggsave(filename = "disc_haz_lnorm.png", plot = p, path = "../plots")


