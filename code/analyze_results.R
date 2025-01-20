source("sim_functions.R")
library(ggplot2)
#Reading in data
df_100 <- readRDS("../data/summary_192_100dfs_list.rds")

df_100 <- bind_rows(df_100)

# Calculate "TRUE NMR, IMR, U5MR" from lnorm_mean and lnorm_sd
# Example
plnorm(60, mean = 12, sd = exp(1.9))



#Wrangling data to utilize Coverage.
coverage_df <- df %>% 
  filter(model != "surv_synth_weibull") %>%
  mutate(true_u5mr = plnorm(60, mean = lnorm_mean, sd = exp(lnorm_sd)),
         true_imr = plnorm(12, mean = lnorm_mean, sd = exp(lnorm_sd)),
         true_nmr = plnorm(1, mean = lnorm_mean, sd = exp(lnorm_sd))) %>%
  mutate(u5mr_covered = ifelse(U5MR_lower <= true_u5mr & U5MR_upper >= true_u5mr, 1, 0),
         nmr_covered = ifelse(NMR_lower <= true_nmr & NMR_upper >= true_nmr, 1, 0),
         imr_covered = ifelse(IMR_lower <= true_imr & IMR_upper >= true_imr, 1, 0)) %>%
  group_by(model, sample_size, proportion) %>%
  summarize(u5mr_coverage = mean(u5mr_covered, na.rm = TRUE),
            imr_coverage = mean(imr_covered, na.rm = TRUE),
            nmr_coverage = mean(nmr_covered, na.rm = TRUE)) %>%
  pivot_longer(cols = u5mr_coverage:nmr_coverage,
               names_to = "Estimate",
               values_to = "Value")%>%
  mutate(sample_size = as.factor(sample_size))

#----------------------------------------------------------------------------

#PLOTS

coverage_df%>%
  group_by(model, Estimate)%>%
  summarize(mean_cov = mean(Value))
#COVERAGE PLOTS

#Labels for facet wrap
labels_wrap <- c(`imr_coverage` = "IMR", `nmr_coverage` = "NMR", `u5mr_coverage` = "U5MR")

#Main coverage plot
coverage_df%>%
  filter(Value != 0)%>%
  filter(model != "logquad_adj")%>%
  ggplot(aes(x = as.factor(sample_size), y = Value, col = model, shape = factor(proportion))) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.95) +
  facet_wrap(~Estimate,labeller = as_labeller(labels_wrap))+
  theme_minimal()+
  labs(x = "Sample Size", y = "Coverage", shape = "Proportion", col = "Model",
       title = "Coverage of NMR, IMR, and U5MR Estimates")+
  scale_color_viridis_d(labels = c("Discrete Hazards", "Log-quad", "Lognormal", "Lognormal adj"), option = "H") +
  coord_cartesian(ylim = c(0.88, 0.98)) +
  scale_y_continuous(breaks = seq(0.88, 0.98, by = 0.01))+
  theme(text = element_text(size = 24))

coverage_df%>%
  filter(model == "discrete_hazards")%>%
  ggplot(aes(x = proportion, y = Value, col = sample_size))+
  geom_point()+
  geom_hline(yintercept = .95, linetype = "dashed")+
  facet_wrap(~Estimate)+
  theme_minimal()

#------------------------------------------------------------------------

#First layer, generating survival curve using lognormal distribution
x <- seq(0, 60, by = .1)
y_12_19 <- 1 - plnorm(x, meanlog = 12, sdlog = exp(1.9))
y_20_23 <- 1 - plnorm(x, meanlog = 20, sdlog = exp(2.3))
y_15_21 <- 1 - plnorm(x, meanlog = 15, sdlog = exp(2.1))
y_15_23 <- 1 - plnorm(x, meanlog = 15, sdlog = exp(2.3))

survival_curve_df <- data.frame(x = x, y_12_19 = y_12_19, y_20_23 = y_20_23, y_15_21 = y_15_21, y_15_23 = y_15_23)
View(survival_curve_df)

colnames(df)[3] <- "NMR_estimate"
colnames(df)[4] <- "IMR_estimate"
colnames(df)[5] <- "U5MR_estimate"

df_2 <- df%>%
  pivot_longer(cols = NMR_estimate:U5MR_upper, names_to = c("measure", ".value"), 
               names_pattern = "(.*)_(.*)")

df_3 <- df_2%>%
  group_by(lnorm_mean, lnorm_sd, measure, proportion)%>%
  summarize(mean_upper = 1 - mean(upper, na.rm = T),
            mean_lower = 1 - mean(lower, na.rm = T), 
            mean_estimate = 1 - mean(estimate, na.rm=T))%>%
  ungroup()%>%
  mutate(x = rep(c(12,1,60), 4))
#lognormal 12, 1.9
survival_curve_df%>%
  ggplot()+
  geom_line(aes(x = x, y = y_12_19, color = "lnorm(mu = 12, sd = 1.9)"))+
  geom_errorbar(data = df_3%>%filter(lnorm_mean == 12), aes(x=x,y = mean_estimate, ymin = mean_lower, ymax = mean_upper, color = measure))+
  theme_minimal()+
  scale_y_continuous(limits = c(.8,1))+
  labs(x = "Age (months)", y = "Survival Rate", title = "Lognormal Survival Curve", color = "True Estimate")

#lognormal 20, 2.3
survival_curve_df%>%
  ggplot()+
  geom_line(aes(x = x, y = y_20_23))+
  theme_minimal()+
  scale_y_continuous(limits = c(.8,1))+
  labs(x = "Age (months)", y = "Survival Rate", title = "Lognormal Survival Curve", 
       subtitle = "Using parameters mu")

#lognormal 15, 2.1
survival_curve_df%>%
  ggplot()+
  geom_line(aes(x = x, y = y_15_21))+
  theme_minimal()+
  scale_y_continuous(limits = c(.8,1))+
  labs(x = "Age (months)", y = "Survival Rate", title = "Lognormal Survival Curve")

#lognormal 15, 2.3
survival_curve_df%>%
  ggplot()+
  geom_line(aes(x = x, y = y_15_23))+
  theme_minimal()+
  scale_y_continuous(limits = c(.8,1))+
  labs(x = "Age (months)", y = "Survival Rate", title = "Lognormal Survival Curve")

















#NOTES
#Heaping index in regards to the proportions of heaping that we use. 
#Lets say proportion is 0.5. If we assume that, regardless of range of heaping (They all include 10-14 months), 
#deaths at each month between 10 and 14 are roughly equal. (Or linear, that works too).
#If we then give half of the deaths at 10,11,13,14 months to month 12, the heap index works out to be 3.
#For proportion of .2, it is 1.8, and for .1, it is 1.4. These numbers are maybe a little lower than what our analysis of
#DHS surveys suggests, although I think it is important to keep these lower values in. However, it might make sense to 
#Add some more proportions, maybe .4, .35, .6, .75, etc. We should also calculate heap indexes for all the simulated dfs,
#to see what the averages are and see if that lines up with the math/assumptions that I make above. 




