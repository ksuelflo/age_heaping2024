# Simulation -- Weibull, two-parameters

library(survival)
library(dplyr)
library(flexsurv)
library(tidyverse)

log_shapes <- c(-1,-1)
log_scales <- c(12,15)

nchild <- 10
nperiods <- 2
length_period <- 60


# Generate data -----------------------------------------------------------

sim_df <- expand.grid(p = 1:2,
                      child = 1:nchild,
                      birth_month = 1:(length_period * 2))
sim_df$id <- rep(1:(nchild * nperiods * length_period), each = 2)
sim_df$child <- NULL

# add l_p
sim_df$l_p <- length_period

# make a_pi
sim_df <- sim_df %>%
  mutate(a_pi = - (birth_month-1) + l_p * (p-1))
view(sim_df)
# simulate who dies in the first time period (must be born before birth_month 60)
p1_df <- sim_df %>% filter(p == 1, a_pi > -l_p)
view(p1_df)
# calculate probability of dying in this time period (no left truncation because first time period)
p1_df <- p1_df %>% mutate(prob = pweibull(a_pi + l_p, scale = exp(log_scales[1]), shape = exp(log_shapes[1])))
p1_df$died <- sapply(p1_df$prob, function(x) {rbinom(n = 1, size = 1, prob = x)})

# for those who died, sample death time
died_p1 <- p1_df %>% filter(died == 1)
died_p1$t <- NA
view(died_p1)


for (i in 1:nrow(died_p1)) {
  t <- 10000
  while(t > (died_p1$a_pi[i] + died_p1$l_p[i])) {
    t <- rweibull(1, shape = exp(log_shapes[1]), scale = exp(log_scales[1]))
  }
  died_p1$t[i] <- t
}

# for those who died in p1, merge this info back into sim_df
p1_df <- rbind(p1_df %>% filter(died == 0) %>% mutate(t = a_pi + l_p), died_p1)
sim_df <- left_join(sim_df, p1_df)
view(sim_df)

# sample deaths in p2 for those who either
# - are born in p2
# - didn't die in p1
who_died <- died_p1$id
which_p2 <- which(!(sim_df$id %in% who_died) & sim_df$p == 2)

# last period, so just sample death from a_pi onward (pretend it can go to inf)
for (i in 1:length(which_p2)) {
  t <- -1
  while((t < (sim_df$a_pi[which_p2[i]])) | (t < 0)) {
    t <- rweibull(1, shape = exp(log_shapes[2]), scale = exp(log_scales[2]))
  }
  sim_df$t[which_p2[i]] <- t
  sim_df$died[which_p2[i]] <- 1
}
view(sim_df)

# Check params ------------------------------------------------------------

test_df <- sim_df %>% filter(!is.na(died))

# Period 1
res_1 <- survreg(formula = Surv(time = t, event = died, type = "right") ~ 1, data = test_df %>% filter(p == 1), dist = "weibull")
log(1/res_1$scale) #Should be -1
log(exp(coef(res_1))) #Should be 12
view(test_df%>%filter(p==1))
# Period 2
res_2 <- flexsurvreg(formula = Surv(time = a_pi, time2 = t, event = died) ~ 1, data = test_df %>%
                       filter(p == 2) %>%
                       mutate(a_pi = ifelse(a_pi < 0, 0, a_pi)), dist = "weibull")
res_2$coefficients[1]
# Should be -1, 15

nrow(test_df%>%filter(p==1))
nrow(test_df %>%
       filter(p == 2) %>%
       mutate(a_pi = ifelse(a_pi < 0, 0, a_pi)))

nrow(test_df%>%filter(p==1)%>%filter(died == 1))

view(test_df %>%
       filter(p == 2) %>%
       mutate(a_pi = ifelse(a_pi < 0, 0, a_pi)))

per_2 <- test_df %>%
  filter(p == 2) %>%
  mutate(a_pi = ifelse(a_pi < 0, 0, a_pi))

mean(per_2$t)
