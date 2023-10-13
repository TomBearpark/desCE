## Some fumbling around with the BSM main regression result

# set up ------------------------------------------------------------------
if(!require(pacman)) install.packages('pacman')
pacman::p_load(fixest, tidyverse, janitor, broom, marginaleffects, useful)
theme_set(theme_bw())

if(Sys.info()['user'] == "tombearpark"){
  root <- "/Users/tombearpark/Dropbox/"  
}else{
  stop("error")
}
db <- file.path(root, "BP_2023_fesearch")
dir.data <- paste0(db, "/data/BurkeHsiangMiguel2015_Replication/data/")
dir.out <- paste0(db, "/out/presentation/")

df   <- read_csv(paste0(dir.data, '/input/GrowthClimateDataset.csv')) %>% 
  # Clean up, add variables (as in their stata code)
  filter(., !is.na(growthWDI), !is.na(UDel_temp_popweight)) %>% 
  mutate(temp1 = UDel_temp_popweight,
         temp2   = temp1 * temp1, 
         precip1  = UDel_precip_popweight / 1000, 
         precip2 = precip1 * precip1,
         time2   = time^2, time3 = time^3, time4 = time^4,
         year2 = year^2) %>% 
  rename(y = growthWDI, country = countryname) 

df.sim <- df %>% 
  select(time1 = time, iso, y, x = temp1, time2 = time2, year, year2) %>% 
  mutate(y = 100 * y)

N  <- nrow(df.sim)
Ni <- length(unique(df.sim$iso))


# are x and y trended -----------------------------------------------------
x.trends <- feols(x ~ -1 + i(iso) + i(iso, time1), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')

y.trends <- feols(y ~ -1 + i(iso) + i(iso, time1), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')
pdf <- bind_rows(
  mutate(tidy(x.trends), var = "Temperature"),  
  mutate(tidy(y.trends), var = "GDP-Growth")
) %>% 
  filter(str_detect(term, "time1"))

pdf %>% 
  ggplot() + 
  geom_rect(aes(xmin=-1.96, xmax=1.96, ymin=-Inf, ymax=Inf), 
            fill = "pink", 
            alpha = .01) + 
  geom_histogram(aes(x = statistic), alpha = .4) + 
  geom_vline(xintercept = 0, color = 'pink') + 
  facet_wrap(~var) + xlab("Time trend t-stat")

ggsave(paste0(dir.out, "time_trends.png"), height = 3, width = 7)

pdf %>% mutate(sig = 1*(abs(statistic)>1.96)) %>% group_by(var) %>% 
  summarize(mean(sig))






#  are they quadratic trended  --------------------------------------------
x.trends2 <- feols(x ~ -1 + i(iso) + i(iso, time1) + i(iso, time2), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')

y.trends2 <- feols(y ~ -1 + i(iso) + i(iso, time1) + i(iso, time2), 
                  df.sim, panel.id = c('iso', 'time1'), 
                  vcov = 'hetero')
pdf2 <- bind_rows(
  mutate(tidy(x.trends2), var = "Temperature"),  
  mutate(tidy(y.trends2), var = "GDP-Growth")
) %>% 
  filter(str_detect(term, "time")) %>% 
  mutate(term = case_when(str_detect(term, "time1") ~ "linear", 
                          str_detect(term, "time2") ~ "quadratic" )) 
pdf2 %>% 
  ggplot() + 
  geom_rect(aes(xmin=-1.96, xmax=1.96, ymin=-Inf, ymax=Inf), 
            fill = "pink", 
            alpha = .01) + 
  geom_histogram(aes(x = statistic), alpha = .4) + 
  geom_vline(xintercept = 0, color = 'pink') + 
  facet_wrap(vars(term, var)) + xlab("Time trend t-stat")

ggsave(paste0(dir.out, "time2_trends.png"), height = 3, width = 7)

pdf %>% mutate(sig = 1*(abs(statistic)>1.96)) %>% group_by(var, term) %>% 
  summarize(mean(sig))




# eda: how do trends in x and y relate ------------------------------------

df.sim
xt <- feols(x ~ -1 + i(iso) + i(time1, ref = "1") + i(iso, time1), 
            df.sim, panel.id = c('iso', 'time1'), 
            vcov = 'hetero')

yt <- feols(y ~ -1 + i(iso) +  i(year, ref = "1") + i(iso, year), 
            df.sim, panel.id = c('iso', 'year'), 
            vcov = 'hetero')

yt.x <- feols(y ~ -1 + x + i(iso) +  i(year, ref = "1") + i(iso, year), 
            df.sim, panel.id = c('iso', 'year'), 
            vcov = 'hetero')

plot.df <- 
  bind_rows(
    left_join(select(tidy(xt), term, xt = estimate), 
              select(tidy(yt), term, yt = estimate)) %>% 
      mutate(y = "raw"),
    left_join(select(tidy(xt), term, xt = estimate), 
              select(tidy(yt.x), term, yt = estimate)) %>% 
      mutate(y = "x control")
    
  ) %>% 
  filter(!str_detect(term, "Intercept")) %>% 
  mutate(type = case_when(
    str_detect(term, "year") & str_detect(term, "iso") ~ "trend", 
    str_detect(term, "year") & !str_detect(term, "iso") ~ "fe time", 
    !str_detect(term, "year") & str_detect(term, "iso") ~ "fe iso",
    .default = NA))


plot.df %>% 
  ggplot() + 
  geom_point(aes(x = xt, y = yt)) + 
  geom_smooth(aes(x = xt, y = yt), method = 'lm') + 
  facet_wrap(vars(y, type), scales = 'free')



 
# Are there level differences in temperature across counties. OFC
m.fe <- feols(x ~ i(iso), 
                df.sim, panel.id = c('iso', 'year'), 
              vcov = 'hetero')

wald(m.fe)

# After removing level differences, are there differences in trends across countries

trends <- i(df.sim$iso, df.sim$year) %>% as_tibble()
names(trends) <- paste0(names(trends), 'trend')
trends.ff     <- paste0(names(trends), collapse = "+")

df.sim.t <- bind_cols(df.sim, trends)
ff <- as.formula(paste0("x ~ i(iso) + ", trends.ff))
m.fe <- feols(ff, 
              df.sim.t, 
              panel.id = c('iso', 'year'), 
              vcov = 'hetero')

wald(m.fe, keep = "trend")

# After removing trends and level differences, is there common shocks
m.fet <- feols(as.formula(paste0("x ~ i(iso) + ", trends.ff, 
                                 "+ i(year)")), 
              df.sim.t, 
              panel.id = c('iso', 'year'), 
              vcov = 'hetero')
wald(m.fet, "year")


# After removing trends and level differences, is there quadratic trend
trends2 <- i(df.sim$iso, df.sim$year2) %>% as_tibble()

m.fet.ar <- feols(as.formula(paste0("x ~ i(iso) + ", trends.ff, 
                                 "+ i(year) + i(iso, l(x))")), 
               df.sim.t, 
               panel.id = c('iso', 'year'), 
               vcov = 'hetero')
wald(m.fet.ar, "l")



# diagnostics on residuals ------------------------------------------------
df.r <- tibble(df.sim, 
          fe = resid(m.fe), 
          fet = resid(m.fet), 
          fet.ar = resid(m.fet.ar, na.rm = FALSE))


