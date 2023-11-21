pacman::p_load(MASS, fixest, tidyverse, broom, furrr)
plan(multisession, workers = 6)
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
    mutate(x.trend = as.vector(tt %*% x.trends), 
           y.trend = as.vector(tt %*% y.trends), 
           x.noise = x.noise, 
           y.noise = y.noise, 
           x = x.trend + x.noise,
           y = beta * x + y.trend + y.noise
             )
  
  m0 <- feols(y ~ x | i, data = df, vcov = "hetero")
  m1 <- feols(y ~ x | i + i[t], data = df, vcov = "hetero")
  
  # Separate CV on each variable
  FEs <- list(' ~ 1 | i', ' ~ 1 | i + i[t]')
  results.x <- run_cv(df, "x", FEs, id.vars = c('i', 't'), test.prop = .2, K = 50)
  results.y <- run_cv(df, "y", FEs, id.vars = c('i', 't'), test.prop = .2, K = 50)
  x.winner  <- arrange(results.x, rmse)[[1,1]]
  y.winner  <- arrange(results.y, rmse)[[1,1]]
  
  # dml estimate
  m.x <- feols(best.model(results.x), df, combine.quick = FALSE)
  m.y <- feols(best.model(results.y), df, combine.quick = FALSE)
  for(var in c('x', 'y')) {
    # df[paste0('hat.', var)] <- predict(get(paste0("m.", var)), df)
    df[paste0('tilde.', var)] <- residuals(get(paste0("m.", var)))
  }
  m.dml <- feols(tilde.y ~ tilde.x, df)
  
  # Standard CV
  FEs <- list(' ~ x | i', ' ~ x | i + i[t]')
  results.yx <- run_cv(df, "y", FEs, id.vars = c('i', 't'), test.prop = .2, K = 50)
  yx.winner  <- arrange(results.yx, rmse)[[1,1]]
  
  # Output the estimates
  bind_rows(
    mutate(broom::tidy(m0),    model = "no trend"), 
    mutate(broom::tidy(m1),    model = "trend"), 
    mutate(broom::tidy(m.dml), model = "dml")
  ) %>% 
    filter(str_detect(term, "x")) %>% 
    mutate(ii = sim.i, 
           x.selected  = x.winner, 
           y.selected  = y.winner, 
           yx.selected = yx.winner)
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
get_cv_winners <- function(sim.df){
  map_dfr(str_subset(names(sim.df), 'selected'), 
          \(x) mutate(tally(group_by(sim.df, get(x))), type = x)) %>% 
    relocate(type)
}


# run simulation ----------------------------------------------------------
Nt   <- 40
Ni   <- 100
NN   <- Ni * Nt

beta <- .1
Nsim <- 50

# Case 1: no trends
sim1 <- future_map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0, 0), 
            trend.Sigma = diag(0, nrow = 2))
  }, 
  .options = furrr_options(seed = seed)
)
sim.plots(sim1)
sim.stats(sim1, beta = beta)
get_cv_winners(sim1)

# Case 2: only x trending
sim2 <- future_map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0.1, 0), 
            trend.Sigma = matrix(c(1, 0, 0, 0), nrow = 2))
  }, 
  .options = furrr_options(seed = seed)
)
sim.plots(sim2)
sim.stats(sim2, beta = beta)
get_cv_winners(sim2)

# Case 3: only y trending
sim3 <- future_map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0, .1), 
            trend.Sigma = matrix(c(0, 0, 0, 1), nrow = 2))
  }, 
  .options = furrr_options(seed = seed)
)
sim.plots(sim3)
sim.stats(sim3, beta = beta)
get_cv_winners(sim3)

# Case 4: uncorrelated trends, centered on zero
sim4 <- future_map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0, 0), 
            trend.Sigma = diag(1, nrow = 2))
  }, 
  .options = furrr_options(seed = seed)
)

sim.plots(sim4)
sim.stats(sim4, beta = beta)
get_cv_winners(sim4)

# Case 5: positive and correlated trends in both
sim5 <- future_map_dfr(
  1:Nsim, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, NN = NN, beta = beta, 
            trend.mu = c(0.1, 0.1), 
            trend.Sigma = matrix(c(1, .5, .5, 1), nrow = 2))
  }, 
  .options = furrr_options(seed = seed)
)

sim.plots(sim5)
sim.stats(sim5, beta = beta)
get_cv_winners(sim5)


