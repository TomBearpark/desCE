# TO DO
# Jointly draw income and climate as a 2D N_i length vector

pacman::p_load(fixest, MASS, tidyverse, broom, furrr)
dir.out <- '/Users/tombearpark/Dropbox/dbce/out/sim/'
theme_set(theme_bw())
set.seed(1)
plan(multisession, workers = 5)
options(dplyr.summarise.inform = FALSE)

# fixed params ------------------------------------------------------------

N_i <- 100
N_y <- 20
N   <- N_i * N_y * 360
df.days <- 
  expand_grid(year = 1:N_y, month = 1:12, day = 1:30, i = 1:N_i) %>% 
  mutate(yday = (month-1) * 30 + day) 

# funcs -------------------------------------------------------------------

gen_daily_data <- function(df.days, 
                           seasonality_aplitude, 
                           eps_mean, eps_var, 
                           alpha_mean, alpha_var,
                           beta, 
                           ...
                           ){
  
  N_i <- length(unique(df.days$i))
  N   <- length(df.days$i)
  
  # Draw stochastic parts 
  alpha_draws <- MASS::mvrnorm(n = N_i, mu = alpha_mean, Sigma = alpha_var)
  colnames(alpha_draws) <- c("alpha_clim_i", "alpha_y_i")
  alpha_draws <- as_tibble(alpha_draws) %>% mutate(i = row_number())
  
  eps_draws <- MASS::mvrnorm(n = N, mu = eps_mean, Sigma = eps_var) 
  colnames(eps_draws) <- c("eps_clim_i", "eps_y_i")
  eps_draws <- as_tibble(eps_draws)
  
  # Temperature DGP
  df.daily <- df.days %>% 
    
    # generate temp seasonality sin wave
    mutate(seasonality = 
             seasonality_aplitude * sin(2 * pi * (.data[["yday"]]-1) / 360)
           ) %>% 
    
    # Add unit specific level differences and idiosyncratic shocks
    left_join(alpha_draws, by = "i") %>% bind_cols(eps_draws) %>%
    
    # Add seasonality in temperature
    mutate(temp1 = seasonality + alpha_clim_i + eps_clim_i, 
           temp2 = temp1^2) %>% 
    
    # Add the i specific coefs
    beta() %>%
  
    # Generate outcomes
    mutate(F.T = beta_1 * temp1 + beta_2 * temp2, 
           
           y   = F.T + alpha_y_i + eps_y_i)
  
  return(df.daily)
}

collapse_annual <- function(df.daily){
  df.daily %>% 
    group_by(i, year) %>% 
      summarize(temp1 = sum(temp1), 
                temp2 = sum(temp2), 
                y = sum(y), 
                .groups = "drop") %>% 
    ungroup() 
}

estimate <- function(df){
  
  m <- feols(y ~ temp1 + temp2 | i, df)
  broom::tidy(m, conf.int=TRUE)
  
}

sim <- function(df.days, 
                seasonality_aplitude, 
                eps_mean, eps_var, 
                alpha_mean, alpha_var,
                beta, ...){
  
  df.daily <- do.call(gen_daily_data, as.list(environment()))
  df       <- collapse_annual(df.daily)
  m        <- estimate(df)
  return(m)
}

map_sim <- function(Nsim, 
                    sim.tag, 
                    df.days, 
                    seasonality_aplitude, 
                    eps_mean, eps_var, 
                    alpha_mean, alpha_var,
                    beta, ...){
  
  t <- Sys.time()
  
  out.df <- 
    future_map_dfr(
    1:Nsim, 
    function(ii){
      sim(df.days = df.days, 
          seasonality_aplitude = seasonality_aplitude, 
          eps_mean = eps_mean, eps_var = eps_var, 
          alpha_mean = alpha_mean, alpha_var = alpha_var,
          beta = beta) %>% 
        mutate(i = !!ii)
    }, 
    .options = furrr_options(seed = 1)
  )
  
  out.df <- mutate(out.df, sim.tag = !!sim.tag)
  
  message(Sys.time() - t)
  
  return(out.df)
}

sim.density <- function(df){
  ggplot(df) + 
    geom_density(aes(x = estimate)) + 
    facet_wrap(~term, scales = "free_x") 
  
}

sim.mean <- function(df) df %>% group_by(term) %>% summarize(mean(estimate))

# diagnostic plots ---------------------------------------------------------

args <- list(
  df.days = df.days,
  seasonality_aplitude = 1,
  beta = function(df) df %>% mutate(beta_1 = .1, beta_2 = .05),
  # Shock parameterisations
  eps_mean = c(0, 0),
  eps_var = matrix(c(1, 0, 0, 3), nrow = 2),
  alpha_mean = c(0, 0),
  alpha_var = matrix(c(1, 0, 0, 1), nrow = 2)
)

df.daily <- do.call(gen_daily_data, args)
df.daily %>% 
  filter(i %in% 1:3, year %in% 1) %>% 
  ggplot() + 
  geom_line(aes(x = yday, y = temp1)) + 
  # geom_line(aes(x = yday, y = y), color = "red") + 
  facet_wrap(~i)
ggsave(paste0(dir.out, "temp.png"), height = 6, width = 15)

df <- collapse_annual(df.daily)

df %>% 
  filter(i %in% 1:3) %>% 
  ggplot() + 
  geom_line(aes(x = year, y = temp1)) + 
  geom_line(aes(x = year, y = y), color = "red") + 
  facet_wrap(~i)

ggplot(df) + geom_point(aes(x = temp1, y = y, color = i))

feols(y ~ temp1 + temp2, df)
feols(y ~ temp1 + temp2 | i, df)

# baseline ----------------------------------------------------------------

args$Nsim <- 300
args$sim.tag <- "Baseline"
out.df <- do.call(map_sim, args)
out.df %>% sim.density()
sim.mean(out.df)

# dependence of income and climate ----------------------------------------
args.inc <- list_modify(
  args, 
  alpha_var = matrix(c(1,1,1,1), nrow = 2)
)
out.df.inc <- do.call(map_sim, args.inc)
out.df.inc %>% sim.density()
sim.mean(out.df)


# het beta ----------------------------------------------------------------
beta.het <- function(df) {
  df %>% mutate(beta_1 = rnorm(n = N, mean = .1, sd = 10), 
                beta_2 = rnorm(n = N, mean = .05)) 
}

args.het   <- list_modify(args, 
                          beta = beta.het, sim.tag = "Random Het")
out.df.het <- do.call(map_sim, args.het)
out.df.het %>% sim.density()
sim.mean(out.df.het)

# het beta, dependent on clim ---------------------------------------------

beta.het.climdep <- function(df) {
  N <- dim(df)[1]
  df %>% mutate(beta_1 = rnorm(n = N, mean = .1, sd = 10), 
                beta_2 = rnorm(n = N, mean = alpha_clim_i)) 
}

args.het.climdep <- list_modify(args, beta = beta.het.climdep, sim.tag = "Clim Het")
out.df.het.climdep <- do.call(map_sim, args.het.climdep)

out.df.het.climdep %>% sim.density()

bind_rows(out.df, out.df.het, out.df.het.climdep) %>% 
  sim.density() + facet_wrap(vars(term, sim.tag), scales = "free")

# heterogenous beta, depending on climate of the location-month



