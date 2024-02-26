pacman::p_load(MASS, fixest, tidyverse, broom, furrr)
# plan(multisession, workers = 6)
seed <- 1
set.seed(1)
theme_set(theme_bw())
code <- '~/Documents/GitHub/desCE/'
source(file.path(code, '/utils/cvFuncs.R'))

# funcs -------------------------------------------------------------------

run_sim <- function(sim.i, 
                    Ni, 
                    Nt, 
                    NN, 
                    beta, 
                    trend.mu, 
                    trend.Sigma
                    ){
  
  message(sim.i)
  
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
    mutate(t2 = t^2, t3 = t^3, 
           x.trend = as.vector(tt %*% x.trends), 
           y.trend = as.vector(tt %*% y.trends), 
           x.noise = x.noise, 
           y.noise = y.noise, 
           x       = x.trend + x.noise,
           x2      = x^2, 
           y       = beta[1] * x + y.trend + y.noise
             )
  if(length(beta)==2) df$y <- df$y + beta[2] * df$x2
  
  m0 <- feols(y ~ x + x2 | i + t, data = df, cluster = "i")
  m1 <- feols(y ~ x + x2 | i + t + i[t], data = df, cluster = "i")
  m2 <- feols(y ~ x + x2 | i + t + i[t] + i[t2], data = df, cluster = "i")
  m3 <- feols(y ~ x + x2 | i + t + i[t] + i[t2] + i[t3], data = df, cluster = "i")
  
  bind_rows(mutate(tidy(m0, conf.int = T), k=0), 
            mutate(tidy(m1, conf.int = T), k=1), 
            mutate(tidy(m2, conf.int = T), k=2), 
            mutate(tidy(m3, conf.int = T), k=3))
}

# parameters --------------------------------------------------------------

Nt   <- 50
Ni   <- 160
NN   <- Ni * Nt

beta <- c(.1, 0.001)
Nsim <- 100

# run simulation ----------------------------------------------------------
# Case 1: no trends
sim1 <- map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0, 0), 
            trend.Sigma = diag(0, nrow = 2))
  }
) %>% 
  mutate(k = paste0("k=", k))

sim1 %>% 
  mutate(inCi = ifelse(conf.low < beta &  conf.high > beta, 1, 0)) %>% 
  group_by(term, k) %>% 
  summarise(coverage      = mean(inCi),
            varBetaHat    = var(estimate), 
            mean_abs_bias = mean(abs(estimate-beta)), 
            se            = mean(std.error)) 

sim1 %>% ggplot() + geom_density(aes(x = estimate, color = k)) + 
  facet_wrap(~term, scales = 'free_x')
sim1 %>% ggplot() + geom_density(aes(x = std.error)) + 
  facet_wrap(~k)


# with correlated trends --------------------------------------------------

# Case 2: correlated positive trends
sim2 <- map_dfr(
  1:Nsim, 
  function(sim.i){    
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0.1, 0.1), 
            trend.Sigma = matrix(c(1, .5, .5, 1), nrow = 2))
  }
)
