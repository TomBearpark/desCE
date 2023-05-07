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

# Check what we would recover if we had full information 
feols(y ~ D, df, vcov = "hetero")
feols(y ~ D | G, df, vcov = "hetero")
feols(y ~ D + G, df, vcov = "hetero")
feols(y ~ -1 + D + G + D*G, df, vcov = "hetero")

# What we would observe in the real world... 
df.obs <- df %>% select(i, y, D)

# B&M ---------------------------------------------------------------------

# Tolerance for convergence
tol <- 0.0001
div <- 100000
# Matrices to store results in 
max.iter <- 1000
beta.mat  <- matrix(nrow = max.iter, ncol = 1)
alpha.mat <- matrix(nrow = max.iter, ncol = 2)

# 1. Initial guess

s = 0
beta.mat[1, ]  <- c(.9)
alpha.mat[1, ] <- c(0,  1)


while(div > tol){
  
  s <- s + 1
  
  beta <- beta.mat[s, ]
  alpha <- alpha.mat[s, ]
  
  # 2. Assignment step
  df.obs$e1 = (df.obs$y - (beta[1]*df.obs$D + alpha[1]))^2
  df.obs$e2 = (df.obs$y - (beta[1]*df.obs$D + alpha[2]))^2
  
  df.obs$g  = ifelse(df.obs$e1 < df.obs$e2, 0, 1)
  
  # 3. Update step 
  m <- feols(y ~ D | g, data = df.obs)
  
  # Store new coefficients
  beta.mat[s+1, ]  <- coef(m)[1]
  alpha.mat[s+1, ] <- drop(fixest::fixef(m)$g)
  
  # Check how close we are to convergence 
  div <- max(c(beta.mat[s,]-beta.mat[s+1,], alpha.mat[s,]-alpha.mat[s+1,]))
  print(div)
  
}

plot(beta.mat[!is.na(beta.mat)])





# Scraps from here onwards!
# trying Gs ---------------------------------------------------------------


df.obs$e <- feols(y ~ D, df.obs)$resid

df$kk <- kmeans(df.obs$e, centers = 2)$cluster - 1
feols(y ~ -1 + D + kk + D*kk, df, vcov = "hetero")
feols(G ~ kk, df)


plot(kk)

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
