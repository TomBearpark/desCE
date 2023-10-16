library(MASS)
library(tidyverse)
library(fixest)
library(broom)
set.seed(1)
theme_set(theme_bw())

# funcs -------------------------------------------------------------------

run_sim <- function(sim.i, 
                    Ni, 
                    Nt, 
                    NN, 
                    beta, 
                    trend.mu, 
                    trend.Sigma
                    ){
  
  # Draw random shocks
  noise   <- mvrnorm(n = NN, mu = c(0, 0), Sigma = diag(1, nrow = 2))
  x.noise <- noise[,1]
  y.noise <- noise[,2]
  
  # Draw trends
  trends   <- mvrnorm(n = Ni, mu = trend.mu, Sigma = trend.Sigma)
  x.trends <- trends[,1]
  y.trends <- trends[,2]
  
  # Create variables
  tt <- matrix(1:Nt)
  
  df <- 
    expand_grid(i = 1:Ni, t = 1:Nt) %>% 
    mutate(x.trend = as.vector(tt %*% x.trends), 
           y.trend = as.vector(tt %*% y.trends), 
           x.noise = x.noise, 
           y.noise = y.noise, 
           x = x.trend + x.noise,
           y = beta * x + y.trend + y.noise
             )
  
  m0 <- feols(y ~ x, data = df, vcov = "hetero")
  m1 <- feols(y ~ x | i[t], data = df, vcov = "hetero")
  
  # Output the estimates
  bind_rows(
    mutate(tidy(m0), model = "no trend"), 
    mutate(tidy(m1), model = "trend")
  ) %>% filter(term == "x")
  
}

sim.stats <- function(sim.output, beta){
  sim.output %>% 
    mutate(dev = estimate - beta) %>% group_by(model) %>% 
    summarize(mean(dev), mean(abs(dev)), median(dev), sqrt(mean(dev^2)))
}

sim.plots <- function(sim.output){
  sim.output %>% 
    pivot_longer(cols = c(estimate, std.error)) %>% 
    ggplot() + geom_density(aes(x = value, color = model)) + 
    facet_wrap(~name, scales = 'free')
}


# run simulation ----------------------------------------------------------
Nt   <- 40
Ni   <- 100
NN   <- Ni * Nt

beta <- .1
Nsim <- 200

# Case 1: no trends
sim1 <- map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0, 0), 
            trend.Sigma = diag(0, nrow = 2))
  }
)
sim.plots(sim1)
sim.stats(sim1, beta = beta)

# Case 2: only x trending
sim2 <- map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0.1, 0), 
            trend.Sigma = matrix(c(1, 0, 0, 0), nrow = 2))
  }
)
sim.plots(sim2)
sim.stats(sim2, beta = beta)

# Case 3: only y trending
sim3 <- map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0, .1), 
            trend.Sigma = matrix(c(0, 0, 0, 1), nrow = 2))
  }
)
sim.plots(sim3)
sim.stats(sim3, beta = beta)

# Case 4: uncorrelated trends, centered on zero
sim4 <- map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0, 0), 
            trend.Sigma = diag(1, nrow = 2))
  }
)

sim.plots(sim4)
sim.stats(sim4, beta = beta)

# Case 5: positive and correlated trends in both
sim5 <- map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0.1, 0.1), 
            trend.Sigma = matrix(c(1, .5, .5, 1), nrow = 2))
  }
)

sim.plots(sim5)
sim.stats(sim5, beta = beta)
