## Some fumbling around with the BSM main regression result

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(fixest, tidyverse, janitor, haven, broom, car,useful)
theme_set(theme_bw())


if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
}else{
  stop("error")
}
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out <- paste0(db, "/out/presentation/")

df   <- read_csv(paste0(data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(., !is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1  = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1,
         time2   = time^2, time3 = time^3, time4 = time^4) %>% 
  rename(y = growthWDI, country = countryname) 


# get parameters for x ----------------------------------------------------

df.sim <- df %>% 
  select(year = time, iso, y, x = temp1, year2 = time2) %>% 
  mutate(y = 100 * y)

N <- nrow(df.sim)
Ni <- length(unique(df.sim$iso))

mx.dgp <- feols(x ~ i(iso) + i(iso, year) + i(iso, l(x)), 
                df.sim, panel.id = c('iso', 'year'))
etable(mx.dgp)


# y~x relationship
feols(y ~ x | iso + year, df.sim)

df.sim$u <- residuals(mx.dgp, na.rm = FALSE)

mu.u <- mean(df.sim$u, na.rm = TRUE)
sd.u <- sd(df.sim$u, na.rm = TRUE)

beta <- .05

# Case 1: True model doesn't have trends
sim.out <- map_dfr(
  
  1:100,
  
  function(ii){
    
    df.sim$u_sim <- rnorm(n = N, mu.u, sd.u)
    df.sim$xsim  <- predict(mx.dgp, df.sim)
    
    # Deal with period 0 impact: TP DP
    df.sim <- df.sim %>% 
      group_by(iso) %>% 
      mutate(xsim = ifelse(is.na(xsim), 
                           dplyr::lead(xsim) - u_sim , 
                           xsim )) %>% 
      ungroup() %>% 
      mutate(ysim = beta * xsim + rnorm(N)) 
    
    m0 <- feols(ysim ~ xsim, df.sim, cluster = 'iso')
    m1 <- feols(ysim ~ xsim | iso[year], df.sim, cluster = 'iso')
    m2 <- feols(ysim ~ xsim | iso[year] + iso[year2], df.sim, cluster = 'iso')
    
    beta <- c(coef(m0)['xsim'], coef(m1)['xsim'], coef(m2)['xsim'])
    se <- c(se(m0)['xsim'], se(m1)['xsim'], se(m2)['xsim'])
    
    tibble(ii = ii, m = 0:2, beta = beta, se = se)
    
  })

sim.out %>% 
  pivot_longer(cols = c('m0','m1')) %>% 
  ggplot() + geom_histogram(aes(value, fill = name))


trends <- 
  broom::tidy(feols(y ~ i(iso, year), df.sim)) %>% 
  filter(str_detect(term, "iso")) %>% 
  mutate(iso = str_remove_all(term, "iso::|:year")) %>% 
  select(iso, trend = estimate)


# simulation --------------------------------------------------------------

beta <- .01




df.sim$eps <- residuals(feols(y ~ i(iso, year), df.sim))

df.sim <-  df.sim %>% 
  left_join(trends) %>% 
  group_by(iso) %>% 
  mutate(xi = mean(x)) %>% 
  ungroup()

ggplot(df.sim) + geom_point(aes(x = xi, y = trend))

df.sim <- df.sim %>% 
  mutate(ysim = beta * x + trend * year + eps)

ggplot(df.sim) + geom_point(aes(x = y, y = ysim))

feols(ysim ~ x, df.sim)
feols(ysim ~ x + iso[year], df.sim)


# simpler -----------------------------------------------------------------
df.sim <- df %>% select(year, iso) %>%
  mutate(t = year - 1970)

Ni <- length(unique(df.sim$iso))
# Nt <- length(unique(df.sim$year))
NN <- length(df.sim$year)
beta <- 1
K <- 100
# Draw x

draws <- map_dfr(1:K, 
                 
                 function(ii){
                   
  df.sim <- df.sim %>% 
    
    group_by(iso) %>% 
    
      mutate(a       = rnorm(n = 1, mean = 0, sd = 1), 
             delta.x = rnorm(n = 1, mean = 0, sd = .1),
             delta.y = rnorm(n = 1, mean = 0, sd = .1)
             
             ) %>% 
    
    ungroup() %>% 
    
    mutate(x     = rnorm(n = NN, mean = a, sd = 1),
           xt    = x + delta.x * t, 
           eps   = rnorm(NN, 0, 1), 
           y     = a + beta * x + eps, 
           y.xt  = a + beta * xt + eps, 
           yt.xt = a + beta * xt  + delta.y * t + eps)
  
  
  m01 <- feols(y ~ x, df.sim, cluster = 'iso') %>% tidy() %>% 
    mutate(est = 'corr')
  m02 <- feols(y ~ x | a, df.sim, cluster = 'iso') %>% tidy() %>% 
    mutate(est = 'fe')
  m03 <- feols(y ~ x | a + a[year], df.sim, cluster = 'iso') %>% tidy() %>% 
    mutate(est = 'fe + trend')
  
  m0 <- bind_rows(m01, m02, m03) %>% 
    filter(term == "x") %>% 
    mutate(model = 'no trends') 
  
  m11 <- feols(y.xt ~ xt, df.sim, cluster = 'iso') %>% tidy() %>% 
     mutate(est = 'corr')
  m12 <- feols(y.xt ~ xt | a, df.sim, cluster = 'iso') %>% tidy() %>% 
    mutate(est = 'fe')
  m13 <- feols(y.xt ~ xt | a[t] + a, df.sim, cluster = 'iso') %>% tidy()  %>% 
    mutate(est = 'fe + trend')
  
  m1 <- bind_rows(m11, m12, m13) %>% 
    filter(term == "xt") %>% 
    mutate(model = 'xt trends')
  
  m21 <- feols(yt.xt ~ xt | a, df.sim, cluster = 'iso') %>% tidy()  %>%
    mutate(est = 'fe')
  m22 <- feols(yt.xt ~ xt | a + a[t],  df.sim, cluster = 'iso') %>% tidy() %>% 
    mutate(est = 'fe + trend')
  
  m2 <- bind_rows(m21, m22) %>% 
    filter(term == "xt") %>% 
    mutate(model = 'xt and yt trends')
  
  bind_rows(m0, m1, m2) %>% 
    mutate(ii = ii)
}
)


draws %>% 
  ggplot() + geom_density(aes(x = estimate, 
                          color = est)) + 
  facet_wrap(~model)
draws %>% 
  filter(model == "no trends") %>% 
  ggplot() + geom_density(aes(x = std.error, color = est))
