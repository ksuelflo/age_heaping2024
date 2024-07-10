library(readstata13)
library(tidyverse)
library(pssst)
library(knitr)
library(SUMMER)
library(survey)
library(ggnewscale)
library(purrr)
library(latex2exp)

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



#read rds
results_loaded <- readRDS("../plots/dhs_estimates_df_for_plot.rds")

# view(results_loaded)
#Getting lnorm distribution.
#Using the same 4 different lnorm parameter settings as in our simulation: (12, 1.9), (20,2.1), (15, 2.1), 20(2.3)
ages <- 0:60
y <- list(plnorm(ages, 15, exp(2.1)), plnorm(ages, 15, exp(2.3)), plnorm(ages, 12, exp(1.9)), plnorm(ages, 20, exp(2.3)))

#Building lnorm df, to use in plots.
lnorm_data <- data.frame(y = list_c(y), age = rep(ages, 4), type = rep(c("15, 2.1", "15, 2.3", "12, 1.9", "20, 2.3"), each = 61))

#Splitting up df into 3. Doing this purely to make 3 different plots, because if we made only one plot, there would be too many lines.
results_subset_1 <- results_loaded%>%
  filter(survey_year < 2011)

results_subset_2 <- results_loaded%>%
  filter(survey_year < 2018 & survey_year >2010)

results_subset_3 <- results_loaded%>%
  filter(survey_year > 2017)

#Survey years less than 2011 plot, with lognormal curves overlayed.
p_2010 <- ggplot()+
  geom_line(data = lnorm_data, aes(x = age, y = y, linetype = type, linewidth = type))+
  geom_line(data = results_subset_1, aes(x = age, y = y, group = interaction(country, survey_year), 
                                       color = interaction(country, survey_year)))+
  labs(x = TeX("$x$"), y = TeX("$_xq_0$"), color = "Country (Survey Year)", linetype = TeX("Lognormal Parameters ($\\mu, \\sigma$)"))+
  scale_linewidth_manual(values = c(1, 1, 1, 1))+
  guides(linewidth = "none")+
  scale_color_viridis_d(labels = rename_color_legend(country = results_subset_1$country, survey_year = results_subset_1$survey_year))+
  scale_linetype(labels = c(TeX("(15, $e^{2.1}$)"), TeX("(15, $e^{2.3}$)"),TeX("(12, $e^{1.9}$)"),TeX("(20, $e^{2.3}$)")))+
  theme_classic()+
  theme(axis.title.y = element_text(size = 15))

p_2010
ggsave(filename = "p_2010.png", plot = p_2010, path = "../plots")

#Ghana, Madagascar, Malawi, and Mauritania plot, with lognormal curves overlayed.
p_2015 <- ggplot()+
  geom_line(data = lnorm_data, aes(x = age, y = y, linetype = type))+
  geom_line(data = results_subset_2, aes(x = age, y = y, group = interaction(country, survey_year), 
                                         color = interaction(country, survey_year)))+
  labs(x = TeX("$x$"), y = TeX("$_xq_0$"), color = "Country (Survey Year)", linetype = TeX("Lognormal Parameters ($\\mu, \\sigma$)"))+
  scale_linewidth_manual(values = c(.5, .5, .5, .5))+
  scale_linetype(labels = c(TeX("(15, $e^{2.1}$)"), TeX("(15, $e^{2.3}$)"),TeX("(12, $e^{1.9}$)"),TeX("(20, $e^{2.3}$)")))+
  scale_color_viridis_d(labels = rename_color_legend(country = results_subset_2$country, survey_year = results_subset_2$survey_year))+
  theme_classic()+
  theme(axis.title.y = element_text(size = 15))

p_2015
ggsave(filename = "p_2015a.png", plot = p_2015, path = "../plots")

#Namibia, Nigeria, and Rwanda plot, with lognormal curves overlayed.
p_2022 <- ggplot()+
  geom_line(data = lnorm_data, aes(x = age, y = y, linetype = type))+
  geom_line(data = results_subset_3, aes(x = age, y = y, group = interaction(country, survey_year), 
                                         color = interaction(country, survey_year)))+
  labs(x = TeX("$x$"), y = TeX("$_xq_0$"), color = "Country (Survey Year)", linetype = TeX("Lognormal Parameters ($\\mu, \\sigma$)"))+
  scale_linewidth_manual(values = c(.5, .5, .5, .5))+
  scale_linetype(labels = c(TeX("(15, $e^{2.1}$)"), TeX("(15, $e^{2.3}$)"),TeX("(12, $e^{1.9}$)"),TeX("(20, $e^{2.3}$)")))+
  scale_color_viridis_d(labels = rename_color_legend(country = results_subset_3$country, survey_year = results_subset_3$survey_year))+
  theme_classic()+
  theme(axis.title.y = element_text(size = 15))

p_2022
ggsave(filename = "p_2022.png", plot = p_2022, path = "../plots")



#Helper functions for three plots.

rename_color_legend <- function(country, survey_year){
  #Getting all of the country year values from plot
  string <- unique(interaction(country, survey_year))
  #replacing dots with space and parenthesis
  in_progress <- str_replace(string, "\\.", " \\(")
  #adding right parenthesis
  finished <- str_c(in_progress, ")")
  #returning completed vector, with any underscores removed (for Burkina Faso)
  return (str_replace_all(finished, "_", " "))
}


