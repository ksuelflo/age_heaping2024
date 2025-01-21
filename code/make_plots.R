source("sim_functions.R")
library(ggplot2)
#Reading in data
df_100 <- readRDS("../data/summary_192_100dfs_list.rds")

df_100 <- bind_rows(df_100)%>%
  mutate(true_u5mr = plnorm(60, mean = lnorm_mean, sd = exp(lnorm_sd)),
         true_imr = plnorm(12, mean = lnorm_mean, sd = exp(lnorm_sd)),
         true_nmr = plnorm(1, mean = lnorm_mean, sd = exp(lnorm_sd)),
         sim_setting = as.factor(sim_setting),
         range_both = str_c(low_range, ", ", high_range))

df_100_imr <- df_100%>%
  group_by(sim_setting, true_imr, proportion, range_both, sample_size)%>%
  summarise(avg_IMR = mean(IMR))%>%
  mutate(diff_IMR = avg_IMR - true_imr)

df_100_imr%>%
  ggplot(aes(x = diff_IMR, fill = as.factor(proportion)))+
  geom_density(alpha = .3)+
  scale_fill_viridis_d()+
  facet_grid(~sample_size)+
  theme_classic

# New facet label names for proportion variable
prop.labs <- c("0% heaped", "10% heaped", "20% heaped", "50% heaped")
names(prop.labs) <- c("0", "0.1", "0.2", "0.5")

# New facet label names for range variable
range.labs <- c("6-18 month range", "8-24 month range", "9-21 month range")
names(range.labs) <- c("6, 18", "8, 24", "9, 21")

df_100 %>%
  filter(lnorm_mean == 20, lnorm_sd == 2.3) %>%
  mutate(bias_imr = IMR - true_imr) %>%
  group_by(sim_setting, true_imr, proportion, range_both, sample_size, model) %>%
  summarise(bias_imr = mean(bias_imr)) %>%
  ggplot(aes(x = sample_size, y = bias_imr, color = model)) +
  geom_point() +
  facet_grid(range_both~proportion,
             labeller = labeller(range_both = range.labs, proportion = prop.labs)) +
  geom_hline(yintercept = 0)+
  scale_color_discrete(labels = c("Discrete Hazards", "Log-quad", "Log-quad adj.", "Lognormal", "Lognormal adj."))+
  # theme_minimal()+
  labs(x = "# of Children Born Each Month", y = "Bias in IMR (Infant Mortality Rate)",
       color = "Model Type")
