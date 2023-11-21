pacman::p_load(MASS, tidyverse, fixest)
set.seed(123)

sim <- function(N.sample, s2.ey, s2.ex, alpha, beta.y, beta.x, 
                N.sim){
  map_dfr(
    1:N.sim, 
    function(ii){
      
      # Draw random data
      x <- rnorm(N.sample, 0, 1)
      e <- mvrnorm(n = N.sample, mu = c(0, 0), Sigma = diag(c(s2.ey, s2.ex)))
      # Propensity regression
      d <- beta.x * x + e[,2]
      # Outcomes model
      y      <- alpha * d + beta.y * x + e[1,]
      
      reg.df <- tibble(y = y, d = d, x = x)
      
      m.structural <- feols(y ~ d + x, reg.df, vcov = "hetero")
      m.y          <- feols(y ~ x, reg.df, vcov = "hetero")
      m.x          <- feols(d ~ x, reg.df, vcov = "hetero")
      
      f.structural <- wald(m.structural, "x", print = FALSE)
      f.y          <- wald(m.y, "x", print = FALSE)
      f.x          <- wald(m.x, "x", print = FALSE)
      
      # Structural decision rule
      if(f.structural$p < .05){
        keep.structual <- 1 
        reg <- m.structural
      } else{
        keep.structual <- 0
        reg <- feols(y ~ d, reg.df, vcov = "hetero")
      }
      
      # Double decision rule
      if(f.structural$p | f.x$p < .05){
        keep.double <- 1 
        reg.double  <- m.structural
      } else{
        keep.double <- 0
        reg.double <- feols(y ~ d, reg.df, vcov = "hetero")
      }
      
      tibble(b1 = coef(reg)['d'], keep1 = keep.structual, 
             b2 = coef(reg.double)['d'], keep2 = keep.double
             )    
    }
  )
}
reformat_output <- function(df){
  bind_rows(
    mutate(select(df, estimate = b1, decision = keep1), 
           rule = "structural"), 
    mutate(select(df, estimate = b2, decision = keep2), 
           rule = "double")
  ) %>% 
    mutate(decision = as.factor(decision))
}


# try stuff ---------------------------------------------------------------

df <- sim(
  N.sample = 300,
  s2.ey    = 1,
  s2.ex    = 1,
  alpha    = 1,
  beta.y   = 1,
  beta.x   = 1,
  N.sim    = 100
)

# Find one where it diverges. Reduce power for structural test
df2 <- sim(
  N.sample = 300,
  s2.ey    = 1,
  s2.ex    = 1,
  alpha    = 1,
  beta.y   = .1,
  beta.x   = 1,
  N.sim    = 100
)

df2 %>% filter(keep1 != 1)
df2 %>% 
  reformat_output() %>% 
  ggplot() + geom_density(aes(x = estimate, color = decision)) + 
  facet_wrap(~rule)

df2 %>% 
  reformat_output() %>% 
  mutate(bias = estimate - 1) %>% 
  group_by(rule) %>% 
  summarise(rmse = sqrt(mean(bias^2)), 
            mad  = mean(abs(bias)))
