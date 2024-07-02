library(readstata13)
library(tidyverse)
library(pssst)
library(knitr)
library(SUMMER)
library(survey)
library(ggnewscale)

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
for (i in seq_along(folders)){
  
  #Getting file path
  path <- str_subset(string = list.files(path = str_c("../../data/", folders[i])), pattern = ".DTA") 
  
  #Reading in data
  births <- read.dta13(str_c("../../data/", folders[i], "/", path), generate.factors = TRUE)
  
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
  
  print(i)
  
  #This DF doesn't work
  if (i == 17){
    next
  }
  
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
  
#Getting rid of Nigeria 2013
results_longer <- results_longer[-(113:118),]

saveRDS(results_longer, "../plots/dhs_estimates_df_for_plot.rds")

#Getting lnorm distribution.
#Using the same 4 different lnorm parameter settings as in our simulation: (12, 1.9), (20,2.1), (15, 2.1), 20(2.3)
ages <- 0:60
y <- list(dlnorm(ages, 15, exp(2.1)), dlnorm(ages, 15, exp(2.3)), dlnorm(ages, 12, exp(1.9)), dlnorm(ages, 20, exp(2.3)))
y_update <- list(rep(0, length(y[[1]])), rep(0, length(y[[1]])),rep(0, length(y[[1]])),rep(0, length(y[[1]])))
y
for (i in seq_along(y)){
  
  for (j in seq_along(y[[i]])){
    #Applying discrete hazards method to lnorm output. 
    y_update[[i]][j] <- 1 - prod(1- y[[i]][1:j])
  }
}

#Building lnorm df, to use in plots.
lnorm_data <- data.frame(y = list_c(y_update), age = rep(ages, 4), type = rep(c("lnorm(15, exp(2.1))", "lnorm(15, exp(2.3))", 
                                                                                "lnorm(12, exp(1.9))", "lnorm(20, exp(2.3))"), each = 61))

#Splitting up df into 3. Doing this purely to make 3 different plots, because if we made only one plot, there would be too many lines.
results_subset_1 <- results_longer%>%
  filter(country %in% c("Burkina_Faso", "Cameroon", "Chad"))

results_subset_2 <- results_longer%>%
  filter(country %in% c("Ghana", "Madagascar", "Malawi", "Mauritania"))

results_subset_3 <- results_longer%>%
  filter(country %in% c("Namibia", "Nigeria", "Rwanda"))

#Burkina_Faso, Cameroon, and Chad plot, with lognormal curves overlayed.
p_bf_cm_td <- ggplot()+
  geom_line(data = lnorm_data, aes(x = age, y = y, linetype = type))+
  geom_line(data = results_subset_1, aes(x = age, y = y, group = interaction(country, survey_year), 
                                       color = interaction(country, survey_year)))+
  labs(x = "age (months)", y = "mortality rate", color = "Country Year", linetype = "lognormal curves", caption = "For the curve 'lnorm(12, exp(1.9))',
  12 is equal to the parameter mu, and exp(1.9) is equal to the parameter sigma. The resulting values are updated using the discrete hazards method.")+
  scale_linewidth_manual(values = c(.5, .5, .5, .5))+
  scale_color_brewer(palette = "Set1")+
  theme_classic()

p_bf_cm_td
ggsave(filename = "p_bf_cm_td.png", plot = p_bf_cm_td, path = "../plots")

#Ghana, Madagascar, Malawi, and Mauritania plot, with lognormal curves overlayed.
p_gh_md_mw_ma <- ggplot()+
  geom_line(data = lnorm_data, aes(x = age, y = y, linetype = type))+
  geom_line(data = results_subset_2, aes(x = age, y = y, group = interaction(country, survey_year), 
                                         color = interaction(country, survey_year)))+
  labs(x = "age (months)", y = "mortality rate", color = "Country Year", linetype = "lognormal curves", caption = "Lognormal curve uses mean of 15, sd of exp(2.2).
       The resulting values are updated using the discrete hazards method.")+
  scale_linewidth_manual(values = c(.5, .5, .5, .5))+
  scale_color_brewer(palette = "Set2")+
  theme_classic()

p_gh_md_mw_ma
ggsave(filename = "p_gh_md_mw_ma.png", plot = p_gh_md_mw_ma, path = "../plots")

#Namibia, Nigeria, and Rwanda plot, with lognormal curves overlayed.
p_na_ng_rw <- ggplot()+
  geom_line(data = lnorm_data, aes(x = age, y = y, linetype = type))+
  geom_line(data = results_subset_3, aes(x = age, y = y, group = interaction(country, survey_year), 
                                         color = interaction(country, survey_year)))+
  labs(x = "age (months)", y = "mortality rate", color = "Country Year", linetype = "lognormal curves", caption = "Lognormal curve uses mean of 15, sd of exp(2.2).
       The resulting values are updated using the discrete hazards method.")+
  scale_linewidth_manual(values = c(.5, .5, .5, .5))+
  scale_color_brewer(palette = "Set3")+
  theme_classic()

p_na_ng_rw
ggsave(filename = "p_na_ng_rw.png", plot = p_na_ng_rw, path = "../plots")






