# methods from Frost & Thompson (2000) J. R. Statistical Soc. A

library(tidyverse)
library(brms)
library(tidybayes)

# building a dataset with a substudy containing repeat observations in X

N <- 1000 # statistical units in study
n <- 50  # statistical units in substudy

mu1 <- mu2 <- 20 # mean of both variables
sd1 <- sd2 <- 5  # sd of both variables
rho <- .9        # correlation between both variables

data <- MASS::mvrnorm(
  N, 
  mu = c(mu1, mu2),
  Sigma = matrix(c(sd1^2, rep(sd1*sd2*rho, 2), sd2^2), byrow = T, nrow = 2)
)

data.df <- as.data.frame(data) %>%
  rename(Y = V1, Xtrue = V2) %>%
  mutate(Id = paste0("ID_", 1:N), .before = Y)

# standard deviation of measurement error on X
sdx <- 3
data.df$Xobs <- data.df$Xtrue + rnorm(N, 0, sdx)

data.df %>%
  pivot_longer(c(Xtrue, Xobs), names_to = "Xtype", values_to = "X") %>%
  mutate(Xtype = factor(Xtype, c("Xtrue","Xobs"))) %>%
  ggplot(aes(X, Y)) + facet_wrap(~ Xtype) +
  geom_point() +
  geom_smooth(method = "lm")

with(data.df, cor(Y, Xtrue))
with(data.df, cor(Y, Xobs))  # weaker correlation = regression dilution


# models for Y based on Xtrue & Xobs
m_xtrue <- brm(Y ~ Xtrue, family = gaussian(), data = data.df)
summary(m_xtrue); plot(m_xtrue); pp_check(m_xtrue, ndraws = 100)
m_xobs <- brm(Y ~ Xobs, family = gaussian(), data = data.df)
summary(m_xobs); plot(m_xobs); pp_check(m_xobs, ndraws = 100)
# => underestimation of beta with Xobs as predictor


# substudy with measurement error
idx_substudy <- sample(data.df$Id, n, replace = F)

data.df.sub <- data.df %>%
  filter(Id %in% idx_substudy) %>%
  mutate(sdx = sdx)

data.df.sub %>%
  ggplot() +
  geom_smooth(method = "lm", aes(Xobs, Y)) +
  geom_pointrange(aes(Xobs, Y, xmin = Xobs-sdx, xmax = Xobs+sdx), 
                  color = "red", size = .3, alpha = .5)
  
nms <- 2 # measurements per subject
data.df.sub2 <- expand_grid(
  data.df.sub, Id.sub = 1:nms
) %>%
  # overwrite Xobs with nms
  mutate(Xobs = rnorm(n(), Xtrue, sdx))
  

data.df.sub2 %>%
  ggplot(aes(y = Y)) +
  geom_pointrange(aes(x = Xtrue, xmin = Xtrue-sdx, xmax = Xtrue + sdx),
                  color = "red", size = .3, alpha = .5) +
  geom_point(aes(x = Xobs), size = 1)

# ICC for x-measurements
model <- brm(
  Xobs ~ 1 + (1|Id),
  family = gaussian(),
  data = data.df.sub2
)
summary(model)
plot(model)

hyp_icc <- "sd_Id__Intercept^2 / (sd_Id__Intercept^2 + sigma^2) = 0"
(hyp_icc <- hypothesis(model, hyp_icc, class = NULL))
plot(hyp_icc)

samples_icc <- spread_draws(model, sd_Id__Intercept, sigma) %>%
  mutate(
    icc = sd_Id__Intercept^2 / (sd_Id__Intercept^2 + sigma^2),
    lambda = 1 / icc
  )

# regression dilution factor
ggplot(samples_icc, aes(lambda)) +
  geom_density() 


# combining different models seems not to do the trick ... continue here
# # combine samples of beta_x from m_xobs with samples for lambda
# samples_beta_true <- spread_draws(m_xtrue, b_Xtrue)
# samples_beta_corr <- spread_draws(m_xobs, b_Xobs) %>%
#   select(b_Xobs) %>%
#   expand_grid(samples_icc %>% select(lambda)) %>%
#   mutate(b_Xobs_corr = b_Xobs * lambda)
# 
# samples_beta_corr %>% 
#   pivot_longer(c(b_Xobs, b_Xobs_corr), 
#                values_to = "b_X", names_to = "method") %>%
#   ggplot(aes(b_X, fill = method)) +
#   geom_density(alpha = .4) +
#   geom_density(aes(b_Xtrue), data = samples_beta_true, alpha = .1,
#                inherit.aes = F)
#   

