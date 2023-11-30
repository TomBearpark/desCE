# Simulation to show why covariate selection on outcomes model doesn't work

pacman::p_load(MASS, fixest, tidyverse, broom, furrr)
plan(multisession, workers = 6)
theme_set(theme_bw())
if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
  code <- '~/Documents/GitHub/desCE/'
}else{
  stop("error")
}
seed <- 123
set.seed(seed)
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")

source(file.path(code, '/utils/cvFuncs.R'))

# calibrate magnitudes ----------------------------------------------------
df   <- read_csv(paste0(dir.data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(., !is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1 = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1, time1 = time, 
         time2   = time^2, time3 = time^3, time4 = time^4) %>% 
  rename(y = growthWDI, 
         x1 = temp1, x2 = temp2, 
         w1 = precip1, w2 = precip2) %>% 
  select(y, x1, x2, w1, w2, iso, time1, time2)

reg0 <- feols(data = df, 
              y ~ x1 + i(iso, ref = "AFG") + i(iso, time1, ref = "AFG"), 
              panel.id = c('time1', 'iso'), 
              se = "hetero")

regY <- feols(data = df, 
              y ~ i(iso, ref = "AFG") + i(iso, time1, ref = "AFG"), 
              panel.id = c('time1', 'iso'), 
              se = "hetero")

regX <- feols(data = df, 
              x1 ~ i(iso, ref = "AFG") + i(iso, time1, ref = "AFG"), 
              panel.id = c('time1', 'iso'), 
              se = "hetero")

trends <- bind_rows(
  mutate(broom::tidy(regX), var = "x1"),
  mutate(broom::tidy(regY), var = "y"), 
  mutate(broom::tidy(reg0), var = "yx")) %>% 
  filter(str_detect(term, ":time1"))

trend.stats <- trends %>% 
  group_by(var) %>% 
  summarize(v = var(estimate), mu = mean(estimate), p = mean(p.value))

trends.xy <- left_join(
  select(filter(trends, var == 'x1'), x = estimate, term), 
  select(filter(trends, var == 'y'),  y = estimate, term), 
) 
covs <- trends.xy %>% 
  summarize(cov.trends = cov(x, y))

beta <- coef(reg0)['x1']

# funcs -------------------------------------------------------------------

run_sim <- function(sim.i, 
                    Ni, 
                    Nt, 
                    beta, 
                    noise.Sigma,
                    
                    trend.mu=NULL, 
                    trend.Sigma=NULL,
                    trends.df=NULL,
                    
                    K = 10
                    ){
  NN <- Ni*Nt
  message(sim.i)
  
  # Draw random shocks
  noise   <- mvrnorm(n = NN, mu = c(0, 0), Sigma = noise.Sigma)
  x.noise <- noise[,1]
  y.noise <- noise[,2]
  
  # Draw trends
  if(is.null(trends.df)){
    trends   <- mvrnorm(n = Ni, mu = trend.mu, Sigma = trend.Sigma)
    x.trends <- trends[,1]
    y.trends <- trends[,2]
  }else{
    trends <- slice_sample(trends.df, n = Ni)
    x.trends <- trends$x
    y.trends <- trends$y
  }

  # Create variables
  tt <- matrix(1:Nt)
  
  df <- expand_grid(i = 1:Ni, t = 1:Nt) %>% 
    mutate(x.trend = as.vector(tt %*% x.trends), 
           y.trend = as.vector(tt %*% y.trends), 
           x.noise = x.noise, 
           y.noise = y.noise, 
           x = x.trend + x.noise,
           y = beta * x + y.trend + y.noise)
  
  # Estimate the two competing models
  browser()
  m0 <- feols(y ~ x | i, data = df, cluster = "i")
  m1 <- feols(y ~ x + i(i, t) | i, data = df, cluster = "i")
  
  # Separate CV on each variable
  FEs <- list(' ~ 1 | i', ' ~ 1 | i + i[t]')
  results.x <- run_cv(df, "x", FEs, id.vars = c('i', 't'), test.prop = .2, K = K)
  results.y <- run_cv(df, "y", FEs, id.vars = c('i', 't'), test.prop = .2, K = K)
  x.winner  <- arrange(results.x, rmse)[[1,1]]
  y.winner  <- arrange(results.y, rmse)[[1,1]]
  
  # Standard CV
  FEs <- list(' ~ x | i', ' ~ x | i + i[t]')
  results.yx <- run_cv(df, "y", FEs, id.vars = c('i', 't'), test.prop = .2, K = K)
  yx.winner  <- arrange(results.yx, rmse)[[1,1]]
  
  # F-test
  f.winner <- ifelse(wald(m1, keep = "t", print = FALSE)$p < 0.05, "trend", "no trend")
  
  # Output the estimates
  bind_rows(
    mutate(broom::tidy(m0), model = "no trend"), 
    mutate(broom::tidy(m1), model = "trend")
  ) %>% 
    filter(str_detect(term, "x")) %>% 
    mutate(ii = sim.i, 
           x.selected  = x.winner, 
           y.selected  = y.winner, 
           yx.selected = yx.winner, 
           f.selected  = f.winner)
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


# globals -----------------------------------------------------------------

Nt   <- 30
Ni   <- 40
Nsim <- 20

# 1. just show we get bias from not including trends -------------------------

# Case 5: positive and correlated trends in both

sim5 <- future_map_dfr(
    1:20, 
    function(sim.i){
      run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, beta = .1, 
              noise.Sigma = diag(1, nrow = 2), 
              trend.mu = c(0.1, 0.1), 
              trend.Sigma = matrix(c(1, .1, .1, 1), nrow = 2))
    }, 
    .options = furrr_options(seed = 123)
  )
  
sim.plots(sim5)
sim.stats(sim5, beta = beta)
get_cv_winners(sim5)

# 2. Show cv doesn't work in reasonable case ------------------------------

# Make it harder to detect

Noise <- matrix(c(var(residuals(regX)), 
                  0, 0, 
                  # var(residuals(regY))
                  # var(df$y), 
                  .3
                  ), 
                nrow = 2)
# Trends.mu <- c(trend.stats$mu[trend.stats$var == "x1"],
#                trend.stats$mu[trend.stats$var == "y"])
# 
# Trends.sigma <- matrix(
#   c(trend.stats$v[trend.stats$var == "x1"],
#     covs$cov.trends, covs$cov.trends, 
#     trend.stats$v[trend.stats$var == "y"]), nrow = 2)

# beta <- coef(reg0)['x1']
beta <- .01
sim6 <- future_map_dfr(
  1:100, 
  function(sim.i){
    run_sim(sim.i = sim.i, Ni = Ni, Nt = Nt, beta = beta, 
            noise.Sigma = Noise, 
            trend.mu = NULL, 
            trend.Sigma = NULL, 
            trends.df = trends.xy
            )
  }, 
  .options = furrr_options(seed = seed)
)

plot.df <- sim6 %>% 
  mutate(across(contains("selected"), 
                ~ifelse(str_detect(.x, "t"), "trend", "no trend"))) %>% 
  mutate(d.selected = ifelse(x.selected == "trend" | y.selected == "trend", 
                             "trend", "no trend"))

# sim.plots(sim6) + 
#   geom_vline(xintercept = beta)
# sim.stats(sim6, beta = beta)
get_cv_winners(sim6)

map_dfr(
  c('yx.selected', 'd.selected', 'f.selected'), 
  function(v){
    plot.df %>% 
      filter(model == .data[[v]]) %>% 
      mutate(mu = mean(estimate))  %>% 
      select(estimate, selection = all_of(v), mu) %>% 
      mutate(rule = v)
  }
) %>% 
  ggplot() + geom_density(aes(x = estimate))+
  geom_vline(xintercept = beta, color = 'red')  +
  geom_vline(aes(xintercept = mu)) + facet_wrap(~rule)
