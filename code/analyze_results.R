
df <- readRDS("../Results/in_progress_summary_list.rds")

df <- bind_rows(df)

# Calculate "TRUE NMR, IMR, U5MR" from lnorm_mean and lnorm_sd
plnorm(60, mean = 12, sd = exp(1.9))

compute_plnorm <- function(q, mean, sd) {
  plnorm(q, mean, sd)
}

df %>% 
  filter(model != "surv_synth_weibull") %>%
  mutate(true_u5mr = plnorm(60, mean = lnorm_mean, sd = exp(lnorm_sd)),
         true_imr = plnorm(12, mean = lnorm_mean, sd = exp(lnorm_sd)),
         true_nmr = plnorm(1, mean = lnorm_mean, sd = exp(lnorm_sd))) %>%
  mutate(u5mr_covered = ifelse(U5MR_lower <= true_u5mr & U5MR_upper >= true_u5mr, 1, 0),
         nmr_covered = ifelse(NMR_lower <= true_nmr & NMR_upper >= true_nmr, 1, 0),
         imr_covered = ifelse(IMR_lower <= true_imr & IMR_upper >= true_imr, 1, 0)) %>%
  group_by(model, sample_size, period, proportion) %>%
  summarize(u5mr_coverage = mean(u5mr_covered, na.rm = TRUE),
            imr_coverage = mean(imr_covered, na.rm = TRUE),
            nmr_coverage = mean(nmr_covered, na.rm = TRUE)) %>%
  pivot_longer(cols = u5mr_coverage:nmr_coverage,
               names_to = "Estimate",
               values_to = "Value") %>%
  ggplot(aes(period, Value, col = model, shape = factor(proportion))) +
  geom_point() +
  geom_hline(yintercept = 0.95) +
  facet_wrap(~Estimate, scales = "free")

# Fix sample sizes
# Remove logquad points from u5mr coverage plot
# make survival curve plot
# figure out how to best incorporate low_range, high_range plot without being overwhelming
# get rid of grouping by period, replace with sample size on x-axis
# make legend, axis labels, etc. look nice
 

