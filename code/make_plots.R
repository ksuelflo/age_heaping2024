source("sim_functions.R")
library(ggplot2)
#Reading in data
df_100 <- readRDS("../data/summary_192_100dfs_list.rds")

df_100 <- bind_rows(df_100)

df_100_imr <- df_100%>%
  mutate(true_u5mr = plnorm(60, mean = lnorm_mean, sd = exp(lnorm_sd)),
         true_imr = plnorm(12, mean = lnorm_mean, sd = exp(lnorm_sd)),
         true_nmr = plnorm(1, mean = lnorm_mean, sd = exp(lnorm_sd)),
         sim_setting = as.factor(sim_setting),
         range_both = str_c(low_range, ", ", high_range))%>%
  group_by(sim_setting, true_imr, proportion, range_both, sample_size)%>%
  summarise(avg_IMR = mean(IMR))%>%
  mutate(diff_IMR = avg_IMR - true_imr)

df_100_imr%>%
  ggplot(aes(x = diff_IMR, fill = as.factor(proportion)))+
  geom_density(alpha = .3)+
  scale_fill_viridis_d()+
  facet_grid(~sample_size)+
  theme_classic

