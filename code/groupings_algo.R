library(tidyverse)
library(fixest)
set.seed(1)

N     <- 10000
beta  <- 1
delta <- 1

df <- tibble(i = 1:N, 
             
             # Treatment effect: can make heterogeneous later 
             beta = beta, 
             
             # Randomly assign covariate
             G = sample(c(0, 1), N, replace = TRUE), 
             
             # Assign treatment, random conditional on. G
             D = stats::rbinom(N, 1, 1 / (G + 1)), 
             
             # Calculate the potential outcomes 
             y0 = rnorm(N, mean = delta * G), 
             y1 = y0 + beta, 
             
             # Calculate the observed outcomes. 
             y = D * y1 + (1-D) * y0
             )

feols(y ~ D, df, vcov = "hetero")
feols(y ~ D + G, df, vcov = "hetero")
feols(y ~ -1 + D + G + D*G, df, vcov = "hetero")

df.obs <- df %>% select(y, D)

# B&M ---------------------------------------------------------------------




# trying Gs ---------------------------------------------------------------

# Candidate grouping vector
df.obs$G.t <- sample(c(0, 1), N, replace = TRUE)

# Calculate residuals
m.tst <- feols(y ~ D + G.t, df.obs, vcov = "hetero")
df.obs$resid <- m.tst$residuals

# KS test on conditional distributions 
g1d1 <- df.obs$resid[df.obs$D == 1 & df.obs$G.t == 1]
g1d0 <- df.obs$resid[df.obs$D == 1 & df.obs$G.t == 0]
g0d1 <- df.obs$resid[df.obs$D == 0 & df.obs$G.t == 1]
g0d0 <- df.obs$resid[df.obs$D == 0 & df.obs$G.t == 0]

p1 <- ks.test(g1d1, g1d0)
p2 <- ks.test(g0d1, g0d0)

list(group = df.obs$G.t, 
     p = tibble(p1 = p1, p2 = p2))