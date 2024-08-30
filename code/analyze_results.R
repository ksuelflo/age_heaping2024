
df <- readRDS("../Results/in_progress_summary_list2.rds")

df <- bind_rows(df)

# Calculate "TRUE NMR, IMR, U5MR" from lnorm_mean and lnorm_sd
plnorm(60, mean = 12, sd = exp(1.9))


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


#Labels for facet wrap
labels_wrap <- c(`imr_coverage` = "IMR", `nmr_coverage` = "NMR", `u5mr_coverage` = "U5MR")

#Main coverage plot
coverage_df%>%
  filter(Value != 0)%>%
  # filter(model == "logquad_adj")%>%
  ggplot(aes(x = as.factor(sample_size), y = Value, col = model, shape = factor(proportion))) +
  geom_point() +
  geom_hline(yintercept = 0.95) +
  facet_wrap(~Estimate, scales = "free", labeller = as_labeller(labels_wrap))+
  theme_minimal()+
  labs(x = "Sample Size", y = "Coverage", shape = "Proportion", col = "Model",
       title = "Coverage of NMR, IMR, and U5MR Estimates Compared to the True Value",
       subtitle = "Shown are various proportions of heaping at 12 months, and models used to estimate NMR, IMR, and U5MR")

coverage_df%>%
  filter(model == "discrete_hazards")%>%
  ggplot(aes(x = proportion, y = Value, col = sample_size))+
  geom_point()+
  geom_hline(yintercept = .95, linetype = "dashed")+
  facet_wrap(~Estimate)+
  theme_minimal()

y_12_19 <- 1 - plnorm(1:60, meanlog = 12, sdlog = exp(1.9))
y_20_23 <- 1 - plnorm(1:60, meanlog = 20, sdlog = exp(2.3))
y_15_21 <- 1 - plnorm(1:60, meanlog = 15, sdlog = exp(2.1))
y_15_23 <- 1 - plnorm(1:60, meanlog = 15, sdlog = exp(2.3))


df2 <- df%>%
  pivot_longer(cols = NMR:U5MR_upper, names_to = c("estimate", "lower", "upper"), names_pattern = "(^[A-Z0-9]+$)(lower$)(upper$)",
               values_to = "value")
View(df2)

string <- c("NMR", "U5MR", "IMR", "NMR_lower", "U5MR_lower", "IMR_lower", "NMR_upper", "U5MR_upper", "IMR_upper")
str_extract_all(string, "lower$")

# Fix sample sizes
# Remove logquad points from u5mr coverage plot
# make survival curve plot
# figure out how to best incorporate low_range, high_range plot without being overwhelming
# get rid of grouping by period, replace with sample size on x-axis
# make legend, axis labels, etc. look nice


