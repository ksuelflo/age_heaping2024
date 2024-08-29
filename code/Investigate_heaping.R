library(survival)
library(dplyr)
library(flexsurv)
library(knitr)
library(pssst)
library(logquad5q0)
library(SUMMER)
library(tidyr)
library(stringr)
library(readstata13)
library(rdhs)
library(foreach)
library(doMC)
library(ggplot2)
library(gghighlight)
library(pssst)

#Getting country codes and country names for DHS countries
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
str(ids)

#Getting all surveys beginning in year 2005 (any survey that starts from 2005 or later), only DHS surveys (not AIDs or malaria, etc.) and 
#Only sub-Saharan African countries. 
survs <- dhs_surveys(surveyYearStart = 2005, surveyType = "DHS", countryIds = c("AO", "BJ", "BT", "BF", "BU", "CM", "CF", "TD", "CG", "CD", "CI", 
                                                                       "ER", "KM", "SZ", "ET", "GA", "GM", "GH", "GN", "KE", "LS", "LB", 
                                                                       "MD", "MW", "ML", "MR", "MZ", "NM", "NI", "NG", "OS", "RW", "SN",
                                                                       "SL", "ZA", "SD", "ST", "TZ", "TG", "UG", "ZM", "ZW"))
#Getting dataset info: Filtering so that it only includes Birth recode files as that is the only thing we want.
datasets <- dhs_datasets(surveyIds = survs$SurveyId, 
                         fileFormat = "flat")%>%
  filter(FileType == "Births Recode")

#Logging in (effectively) to DHS
set_rdhs_config(email = "ksuelflo@macalester.edu",
                project = "Evaluation of Existing Statistical Methods Regarding Age Heaping")

# download datasets. This puts them in a cache somewhere (don't quite understand this part). `downloads` is a list of filepaths corresponding
# to each RDS file (df).

downloads <- get_datasets(datasets$FileName)

#Setting up containers
years_df <- rep(0, 2)
#skipping 47 for now: It has an error with both v023 and v024 strata. And skipping 53, 83

df_container <- vector("list", length = 2)
indicator <- 1
for (i in c(52, 81)){
  
  #Reading in data
  print(i)
  births <- readRDS(downloads[[i]])
  
  #Survey year
  year <- min(births$v007)
  years_df[indicator] <- year
  
  tryCatch(
    expr = {
      df_birth <- format_dhs(df = births, survey_year = year, period_boundaries = c(year-5, year-4, year-3, year-2, year-1, year))
      print(str_c("finished number ", i, " download.", " Using strata(v023)"))
    },
    error = function(e){
      df_birth <- format_dhs(df = births, survey_year = year, period_boundaries = c(year-5, year-4, year-3, year-2, year-1, year), strata = c("v024", "v025"))
      print(str_c("finished number ", i, " download.", " Using strata(v024, v025)"))
    }
  )
  df_container[[indicator]] <- df_birth
  indicator <- indicator + 1
}

#Using `ids` and `downloads` to scrape out the country name for each survey.

countries <- rep("", length(downloads))
#All 90 codes for each df
ssa_codes <- str_sub(names(downloads), 1, 2)

#Looping through each survey, and finding corresponding country for the ith code.
for (i in seq_along(countries)){
  code <- ssa_codes[i]
  countries[i] <- ids$CountryName[ids$DHS_CountryCode == code]
}

#Container for summary deaths by month for each survey
summary_deaths_2 <- vector(mode = "list", length = length(df_container))
#This for loop wrangles data into monthly deaths by period. Included are all intervals of death.
for (i in seq_along(df_container)){
  #Checking dfs for NULL: There are many that had error with `format_dhs`, so they are null. We will skip those.
  if (is.null(df_container[[i]])){
    next
  }
  
  print(i)
  summary_heap_df <- df_container[[i]]%>%
    filter(I_i ==1)%>%
    mutate(avg_death = (t_0i + t_1i)/ 2,
           death_date = avg_death*30 + b_i,
           died_in_period = if_else(y_p <= death_date & death_date < (y_p + 365), 1, 0))%>%
    group_by(p)%>%
    summarize(int_0_1 = sum(0 <= avg_death & avg_death < 1 & died_in_period == 1),
              int_1_2 = sum(1 <= avg_death & avg_death < 2 & died_in_period == 1),
              int_2_3 = sum(2 <= avg_death & avg_death < 3 & died_in_period == 1),
              int_3_4 = sum(3 <= avg_death & avg_death < 4 & died_in_period == 1),
              int_4_5 = sum(4 <= avg_death & avg_death < 5 & died_in_period == 1),
              int_5_6 = sum(5 <= avg_death & avg_death < 6 & died_in_period == 1),
              int_6_7 = sum(6 <= avg_death & avg_death < 7 & died_in_period == 1),
              int_7_8 = sum(7 <= avg_death & avg_death < 8 & died_in_period == 1),
              int_8_9 = sum(8 <= avg_death & avg_death < 9 & died_in_period == 1),
              int_9_10 = sum(9 <= avg_death & avg_death < 10 & died_in_period == 1),
              int_10_11 = sum(10 <= avg_death & avg_death < 11 & died_in_period == 1),
              int_11_12 = sum(11 <= avg_death & avg_death < 12 & died_in_period == 1),
              int_12_13 = sum(12 <= avg_death & avg_death < 13 & died_in_period == 1),
              int_13_14 = sum(13 <= avg_death & avg_death < 14 & died_in_period == 1),
              int_14_15 = sum(14 <= avg_death & avg_death < 15 & died_in_period == 1),
              int_15_16 = sum(15 <= avg_death & avg_death < 16 & died_in_period == 1),
              int_16_17 = sum(16 <= avg_death & avg_death < 17 & died_in_period == 1),
              int_17_18 = sum(17 <= avg_death & avg_death < 18 & died_in_period == 1),
              int_18_19 = sum(18 <= avg_death & avg_death < 19 & died_in_period == 1),
              int_19_20 = sum(19 <= avg_death & avg_death < 20 & died_in_period == 1),
              int_20_21 = sum(20 <= avg_death & avg_death < 21 & died_in_period == 1),
              int_21_22 = sum(21 <= avg_death & avg_death < 22 & died_in_period == 1),
              int_22_23 = sum(22 <= avg_death & avg_death < 23 & died_in_period == 1),
              int_23_24 = sum(23 <= avg_death & avg_death < 24 & died_in_period == 1),
              int_24_36 = sum(24 <= avg_death & avg_death < 36 & died_in_period == 1),
              int_36_48 = sum(36 <= avg_death & avg_death < 48 & died_in_period == 1),
              int_48_60 = sum(48 <= avg_death & avg_death < 60 & died_in_period == 1),
              int_60_inf = sum(60 <= avg_death & avg_death < Inf & died_in_period == 1))%>%
    mutate(country = countries[i], #Use `i` instead when in the for loop.
           survey_year = years_df[i], 
           .before = p)%>% #Use `i` instead when in the for loop.
    mutate(p_begin = survey_year - 6 + p,
           p_end = p_begin + 1,
           .after = p)
  
  summary_deaths_2[[i]] <- summary_heap_df
}

saveRDS(summary_deaths, file = "../output/summary.rds")

summary_deaths <- readRDS(file = "../output/summary.rds")

#Now to calculate heaping indexes. Some info about heaping index

#Hill and Choi 2006 index:
# Note that they used 5-9 days, 7 days as the heaped day. This was for early neonatal deaths. However, we are
# extrapolating this method to use at 12 months. 
# Divide the deaths at 12 by 1/5 the total deaths from ages 10-14 months. Can experiment with different ranges here. 

#Pullum and Becker 2014 index:
# Calculate the average deaths from 10-14 months. Then the relative difference between the actual number of deaths
# at 12 months and the average, times 100, is the index. For example, if the average deaths was 50, and the actual
# number was 150, then then index would be 200 (200% higher).

#loop that calculated all the heap indexes for each survey.

#container for new dfs

#START HERE 
# summary_deaths <- readRDS(file = "../output/summary.rds")

indexed_dfs <- vector(mode = "list", length = length(summary_deaths))

#Calculating heap indexes
for (i in seq_along(indexed_dfs)){
  #Checking if the df didn't read in, we skip it.
  # if (is.null(summary_deaths[[i]])){
  #   next
  # }
  # 
  df <- summary_deaths[[i]]%>%
    mutate(sample_size = int_10_11 + int_11_12 + int_12_13 + int_13_14 + int_14_15,
           hillChoi1014 = int_12_13*5/(sample_size),
           pullum1014 = (int_12_13 - ((sample_size)/5))/((sample_size)/5)*100,
           standardIndex = (int_12_13 - (sample_size/5))/sample_size)
  indexed_dfs[[i]] <- df
}
#Setting up folders for plots and tables
#Also making a container for folders to make it easier

all_surveys <- bind_rows(indexed_dfs)

#Fixing Ethiopia DONE
# all_surveys$survey_year[all_surveys$country == "Ethiopia"] <- all_surveys$survey_year[all_surveys$country == "Ethiopia"] + 8
# all_surveys$p_begin[all_surveys$country == "Ethiopia"] <- all_surveys$survey_year[all_surveys$country == "Ethiopia"] + 8
# all_surveys$p_end[all_surveys$country == "Ethiopia"] <- all_surveys$survey_year[all_surveys$country == "Ethiopia"] + 8

all_surveys <- all_surveys%>%
  dplyr::select(-c(int_24_36, int_36_48, int_48_60, int_60_inf))


all_surveys_combined <- all_surveys%>%
  group_by(country, survey_year)%>%
  summarize(int_10_11_sum = sum(int_10_11),
            int_11_12_sum = sum(int_11_12),
            int_12_13_sum = sum(int_12_13),
            int_13_14_sum = sum(int_13_14),
            int_14_15_sum = sum(int_14_15),
            sample_size_sum = sum(sample_size))%>%
  ungroup()%>%
  mutate(standardIndex = (int_12_13_sum - (sample_size_sum/5))/sample_size_sum,
         hillChoi1014 = int_12_13_sum*5/(sample_size_sum))


latex_table <- all_surveys_combined%>%
  dplyr::select(-c(int_10_11_sum, int_11_12_sum, int_12_13_sum, int_13_14_sum, int_14_15_sum, standardIndex))%>%
  mutate(hillChoi1014 = round(hillChoi1014, 2))


  
kable(latex_table, format = "latex")

#Spaghetti plot: BUST
all_surveys%>%
  ggplot(aes(x = p_begin, y = hillChoi1014))+
  geom_point(aes(color = country))+
  geom_line(aes(group = interaction(country, survey_year), color = country))+
  theme_minimal()

#Plot of index vs sample size
p_1 <- all_surveys%>%
  ggplot(aes(x = sample_size, y = hillChoi1014))+
  geom_point()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_minimal()+
  labs(y = "Heap Index", title = "Heap Index vs Sample Size Across 87 sub-Saharan Africa Surveys",
       x = "Sample Size", subtitle = "Each survey contains five one year long periods, thus there is five data points per survey.")

ggsave(filename = "../plots/p_1.png", bg = "white", plot = p_1, width = 300, height = 300, units = "px")
#Plot of standardized index vs sample size
#Coloring by country is a little crazy
p_2 <- all_surveys%>%
  ggplot(aes(x = sample_size, y = standardIndex))+
  geom_point()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_minimal()+
  labs(y = "Standardized Heap Index", title = "Standardized Heap Index vs Sample Size Across 87 sub-Saharan Africa Surveys",
       x = "Sample Size", subtitle = "Each survey contains five one year long periods, thus there is five data points per survey.")

ggsave(filename = "../plots/p_2.png", bg = "white", plot = p_2, width = 300, height = 300, units = "px")

#Plot of standardized index vs sample size, combining periods
p_3 <- all_surveys_combined%>%
  ggplot(aes(x = sample_size_sum, y = standardIndex))+
  geom_point(aes(color = country))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_minimal()+
  labs(y = "Standardized Heap Index", title = "Standardized Heap Index vs Sample Size Across 87 sub-Saharan Africa Surveys",
       x = "Sample Size", subtitle = "Each survey contains one five year long period.")

ggsave(filename = "../plots/p_3.png", bg = "white", plot = p_3, width = 300, height = 300, units = "px")

#Manual color scale
cols <- c("Niger" = "#32ECC2", "Nigeria" = "#FFC107", "Tanzania" = "#1E88E5", "Senegal" = "#D81B60")

p_4 <- all_surveys_combined%>%
  ggplot(aes(x = survey_year, y = hillChoi1014, size = sample_size_sum, color = country))+
  geom_point()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_minimal()+
  scale_color_manual(values = cols, na.value = "black")+
  labs(x = "Year", y = "Heap Index", size = "Total Deaths\n10-14 months", title = "Heap Indexes in sub-Saharan Africa Surveys since 2005",
       country = "Country")

ggsave(filename = "../plots/p_4.png", bg = "white", plot = p_4, width = 1000, height = 800, units = "px")

#--------------------------------------------------------------------
#Making plots (ALREADY MADE, SKIP)
plot_folders <- rep("", length(indexed_dfs))
for (i in seq_along(indexed_dfs)){
  #dir.create(path = str_c("../plots/", str_replace_all(countries[i], pattern = "\\s", replacement = "-"), "_", years_df[i]))
  plot_folders[i] <- str_c("../plots/", str_replace_all(countries[i], pattern = "\\s", replacement = "-"), "_", unique(indexed_dfs[[i]]$survey_year))
}

#Making a vector of ages that is numeric, to help make plots:
numeric_ages <- rep((1:24), 5)


#Making plots for every period within every survey
for (i in seq_along(indexed_dfs)){
  df <- indexed_dfs[[i]]
  if (is.null(df)){
    next
  }
  for (j in 1:5){
    df_periods <- df%>%
      dplyr::select(-c(int_60_inf, int_24_36, int_36_48, int_48_60))%>%
      pivot_longer(names_to = "interval", values_to = "deaths", cols = !c(country, survey_year, p_begin, p_end, p, hillChoi1014, pullum1014))%>%
      mutate(age = numeric_ages, 
             p_begin = as.factor(p_begin))
    
    p_death_dist <- df_periods%>%
      filter(p == j)%>%
      ggplot(aes(x = age, y = deaths))+
      geom_line()+
      geom_point(size = 1.5)+
      gghighlight(age ==13)+
      labs(title = str_c(countries[i], " ", unique(df_periods$survey_year)), subtitle = str_c("period ", j, " (", unique(df$p_begin), " - ", unique(df$p_end), ")"),
           x = "Age at death (months)", y = "total")+
      theme_minimal()+
      theme(plot.title = element_text(hjust = 0.5, size = 20), 
            plot.subtitle = element_text(hjust = .5), 
            panel.background = element_blank())
    
    ggsave(plot = p_death_dist, filename = str_c(plot_folders[i], "/", "period_", j, "plot.png"), width = 6, height = 6, units = "in", bg = "#ffffff")
  }
  
  p_all_periods <- df_periods%>%
    ggplot(aes(x = age, y = deaths))+
    geom_line(aes(group = p_begin, color = p_begin), alpha = .5)+
    geom_point(aes(color = p_begin), size = 1.5, alpha = .5)+
    geom_point(data = df_periods%>%filter(age == 13))+
    labs(title = str_c(countries[i], " ", unique(df_periods$survey_year)), x = "Age at death (months)", y = "total",
         color = "Year")+
    scale_colour_viridis_d()+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          plot.subtitle = element_text(hjust = .5),
          panel.background = element_blank())
  
  ggsave(plot = p_all_periods, filename = str_c(plot_folders[i], "/all_periods_plot.png"), width = 6, height = 6, units = "in", bg = "#ffffff")
  
  hill_p <- df%>%
    ggplot(aes(x = p_end, y = hillChoi1014))+
    geom_point(size = 1.5)+
    geom_line()+
    geom_hline(yintercept = 1, linetype = "dashed")+
    labs(title = str_c(countries[i], " ", unique(df_periods$survey_year)), subtitle = "Heap index across five one year periods",
         y = "Heap Index", x = "year")+
    scale_y_continuous(limits = c(0, max(df$hillChoi1014) + .5))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          plot.subtitle = element_text(hjust = .5),
          panel.background = element_blank())
  
  ggsave(plot = hill_p, filename = str_c(plot_folders[i], "/heapindex.png"), width = 6, height = 6, units = "in", bg = "#ffffff")
  
  
}

